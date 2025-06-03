#!/usr/bin/env python
import pandas as pd
import numpy as np
import os
import argparse
import sys
import time
import multiprocessing
from tqdm import tqdm
import subprocess
import json
import pybedtools

script_path=sys.argv[0]
Mapit_dir=os.path.dirname(script_path)

def argsPrint(args):
    width=max([15]+[len(str(vars(args)[k])) for k in list(vars(args).keys())[:-1]])
    print(format("//"+"="*13, '<'), "MAPIT-seq Pipeline",format("="*(width-2)+"\\\\", '<'))
    print(format("||", '<2'), format("subcommand", '>20'),": ", format('\033[34m'+sys.argv[1]+'\033[0m', '<%d'%(width+14)), format("||", '<2'))
    for arg in list(vars(args).keys())[:-1]:
        print(format("||", '<2'), format(arg, '>20'),": ", format('\033[34m'+str(getattr(args, arg))+'\033[0m', '<%d'%(width+14)), format("||", '<2'))
    print(format("\\\\"+"="*(20+width+11)+"//", '<'))

def infoPrint(type,description):
    colorType={"ERROR":"31","WARN":"33","INFO":"34"}
    print("\033[%sm[%s]\033[32m\t%s\033[0m\t%s"%(colorType[type],type,time.asctime( time.localtime(time.time())),description))
    if type=="ERROR":
        exit(1)

def fileExistCheck(fileCheckList,confFile=""):
    exit_signal=0
    for file in fileCheckList:
        if os.path.exists(file) and os.path.getsize(file)>0:
            continue
        elif confFile:
            exit_signal=1
            infoPrint("WARN","\033[34m%s\033[0m does not exist. Check \033[34m%s\033[0m or regenerate it."%(file,confFile))
        else:
            exit_signal=1
            infoPrint("WARN","\033[34m%s\033[0m does not exist. Check it."%(file))
    return exit_signal

def loadConfig(genomeVersion):
    confFile=Mapit_dir+'/conf/'+genomeVersion+'.json'
    if os.path.exists(confFile):
        with open(confFile) as f:
            conf = json.load(f)
            if genomeVersion != conf['genome']['version']:
                infoPrint(type="WARN",description="Genome Version %s in %s does not match with file name."%(conf['genome']['version'],confFile))
        fileCheckList=[]
        for var in list(conf['genome'].values())[1:]+list(conf['SNP'].values())+list(conf['annotation'].values()):
            fileCheckList.append(var)
        for suffix in ["amb","ann","bwt","pac","sa"]:
            for var in [conf['genome']['abundantRNA'],conf['genome']['bwaIndex']]:
                fileCheckList.append(var+"."+suffix)
        for suffix in range(1,9):
            fileCheckList.append(conf['genome']['hisat2Index']+"."+str(suffix)+".ht2")
        if fileExistCheck(fileCheckList):
            infoPrint("ERROR","Lack of files.\nAborting")
        infoPrint(type="INFO",description="Mapit-seq config file \033[34m%s\033[0m is loaded and has no error."%(confFile))
        return conf
    else:
        infoPrint(type="ERROR",description="\033[34m%s\033[0m does not exist. Use 'Mapit config' to create. "%(confFile))

def isosplit(path):
    df_class=pd.read_csv(path+"/classification.txt",sep="\t")
    df_tag=df_class[["#read_id","gene_id","isoform_id","assignment_type","assignment_events"]].drop_duplicates()
    df_tag['gene_id']="UG:Z:"+df_tag['gene_id']
    df_tag['isoform_id']="UT:Z:"+df_tag['isoform_id']
    df_tag['assignment_type']="AT:Z:"+df_tag['assignment_type']
    df_tag['assignment_events']="AE:Z:"+df_tag['assignment_events'].apply(lambda x:x.replace(":","_"))
    os.system('''samtools view -@ 20 %s/align.bam | \
            awk -F "\t" 'NF==24 {print > "%s/align.txt" } NF==25 {print > "%s/align_SA.txt" } NF>25 {print > "%s/align_FA.txt" }' '''%(
                path,path,path,path))
    os.system("samtools view -@ 20 {path}/align.bam -H > {path}/align_tag.sam".format(path=path))
    infoPrint("INFO","Chromosome %s - adding tag is done!"%(path.split("/")[-1]))
    df_align=pd.read_csv(path+"/align.txt",sep="\t",header=None)
    pd.merge(df_align,  df_tag,left_on=0,right_on="#read_id").drop(columns=["#read_id"]).to_csv(path+"/align_tag.sam", mode='a', header=False,index=False,sep="\t")
    if os.path.exists(path+"/align_SA.txt"):
        df_alignsa=pd.read_csv(path+"/align_SA.txt",sep="\t",header=None)
        pd.merge(df_alignsa,df_tag,left_on=0,right_on="#read_id").drop(columns=["#read_id"]).to_csv(path+"/align_tag.sam", mode='a', header=False,index=False,sep="\t")
    if os.path.exists(path+"/align_FA.txt"):
        df_alignfa=pd.read_csv(path+"/align_FA.txt",sep="\t",header=None,on_bad_lines ="warn")
        pd.merge(df_alignfa,df_tag,left_on=0,right_on="#read_id").drop(columns=["#read_id"]).to_csv(path+"/align_tag.sam", mode='a', header=False,index=False,sep="\t")
    os.system("samtools view -@ 20 {path}/align_tag.sam -hb | samtools sort -@ 20 -o {path}/align_tag.bam".format(path=path))
    os.system("bamtools split -tag UG -in {path}/align_tag.bam".format(path=path))
    infoPrint("INFO","Chromosome %s - spliting genes is done!"%(path.split("/")[-1]))
    os.system("rm {path}/align*.txt {path}/align_tag.sam".format(path=path))
    bamlist = list(filter(lambda x:x.find("TAG")>=0,os.listdir(path+"/")))
    translist=[]
    for gbam in bamlist:
        gene=gbam[17:-4]
        os.makedirs("{path}/{gene}".format(path=path,gene=gene))
        os.system("bamtools split -stub {path}/{gene}/ -in {path}/{bam} -tag UT".format(gene=gene,bam=gbam,path=path))
        tbamlist=list(filter(lambda x:x.find(".bam")>=0 and x.find(".bai")<0,os.listdir("{path}/{gene}".format(path=path,gene=gene))))
        translist=translist+[path+"/"+gene+"/"+tbam for tbam in tbamlist ]
    return translist

def isoredi(tbam):
    os.system("samtools index "+tbam)
    path=os.path.dirname(tbam)
    tbamfile=os.path.basename(tbam)
    command="python /home/gangx/apps/reditools2.0/src/cineca/reditools.py -f {tbam} -q 0 -bq 20 \
            -r {fasta} -o {prefix}_reditools2.txt".format(tbam=tbam, prefix=path+"/"+tbamfile[8:-4],fasta=genomeFasta,)
    out=open(path+"/"+tbamfile[8:-4]+"_reditools2.log",'w')
    ps=subprocess.Popen(command,shell=True,stdout=out,stderr=out)
    signal=ps.wait()
    out.close()
    
def flatten_list(nested_list):
    return [item for sublist in nested_list for item in (flatten_list(sublist) if isinstance(sublist, list) else [sublist])]
    
def exonredi(path,sample,transcript_path,exon_bed):
    df_redi=pd.read_csv(path+"/"+sample+"/"+transcript_path+"_reditools2.txt",sep="\t")
    base=['A','C','G','T']
    df_redi["Start"]=df_redi["Position"]-1
    pybedtools.set_tempdir(os.path.dirname(path+"/"+sample+"/"+transcript_path))
    redi_bed=pybedtools.BedTool.from_dataframe(df_redi[["Region","Start","Position","Reference","Coverage-q30","MeanQ","BaseCount[A,C,G,T]"]])
    redi_bed_exon=redi_bed.intersect(exon_bed,wa=True,wb=True)
    df_redi_exon=redi_bed_exon.to_dataframe(names=["Chrom","Start","Position","REF","Coverage","MeanQ","BaseCount"]+["chrom","exonstart","exonend","exonid","exonorder","strand"])
    if df_redi_exon['strand'].unique()[0]=="+":
        df_redi_exon_edit=df_redi_exon[df_redi_exon['REF'].isin({"A","C"})].copy()
    elif df_redi_exon['strand'].unique()[0]=="-":
        df_redi_exon_edit=df_redi_exon[df_redi_exon['REF'].isin({"T","G"})].copy()
    else:
        exit(1)
    df_redi_exon_edit['ALT']=df_redi_exon_edit['REF'].replace({"A":"G","C":"T","T":"C","G":"A"})
    df_redi_exon_edit['EDIT']=df_redi_exon_edit.apply(lambda x:int(x['BaseCount'][1:-1].split(", ")[base.index(x["ALT"])]),axis=1)
    df_redi_exon_edit['EDIT%']=df_redi_exon_edit['EDIT']/df_redi_exon_edit['Coverage']
    df_redi_exon_edit=df_redi_exon_edit.set_index(["Chrom","Start","Position","REF","ALT","strand"])[["Coverage","EDIT","EDIT%"]].rename(columns={"Coverage":sample+"_DP","EDIT":sample+"_EDIT","EDIT%":sample+"_EDIT%"})
    pybedtools.cleanup()
    return df_redi_exon_edit

def getreditranscript(path):
    translist=[]
    for chr in os.listdir(path+""):
        bamlist = list(filter(lambda x:x.find("TAG")>=0,os.listdir(path+"/"+chr)))
        for gbam in bamlist:
            gene=gbam[17:-4]
            tlist=list(filter(lambda x:x.find(".txt")>=0,os.listdir("{path}/{gene}".format(path=path+"/"+chr,gene=gene))))
            translist=translist+tlist
    return translist

def mergeredi(transcript_path):
    transcript=transcript_path.split("/")[-1]
    pybedtools.set_tempdir("/home/gangx/data/CQX/pacbio/isomapit2")
    exon_bed=pybedtools.BedTool.from_dataframe(df_tbed[df_tbed['transcriptid']==transcript].iloc[:,:6])
    df_edit=pd.concat([exonredi(path="/home/gangx/data/CQX/pacbio/isomapit2",sample=sample,transcript_path=transcript_path,exon_bed=exon_bed) for sample in sampleList],axis=1)
    df_edit=df_edit[df_edit[[sample+"_EDIT%" for sample in sampleList]].sum(1)>0].copy()
    df_edit['transcript']=transcript
    pybedtools.cleanup()
    return df_edit.reset_index()

def isoeditjob(args):
    conf = loadConfig(args.genomeVersion)
    global genomeFasta
    genomeFasta=conf['genome']['fasta']
    chrlist=pd.read_csv(genomeFasta+".fai",sep="\t",header=None)[0].tolist()
    chrlist=list(set(chrlist)-{"chrM"})
    outpath=args.outpath+"/"+args.dirname
    pathlist=[outpath+"/"+chr for chr in chrlist]
    
    df_p=pd.read_csv(args.readAssignments,sep="\t",header=2)
    #os.system("awk ")
    df_transcript=pd.read_csv("/home/gangx/data/CQX/pacbio/gencode.v40.transcript.tsv",sep="\t",header=None,names=["geneid","genesymbol","genetype","transcriptid","transcript","transcripttype"])
    df_pt = pd.merge(df_p[df_p['assignment_type'].isin({"unique","unique_minor_difference"})],
                     df_transcript,left_on=["gene_id","isoform_id"],right_on=["geneid","transcriptid"])
    df_pt_gb=df_pt.groupby("isoform_id")["#read_id"].count().reset_index()
    df_pt_gb_filter=df_pt_gb[df_pt_gb['#read_id']>=10].copy()
    df_pt_filter=df_pt[df_pt['isoform_id'].isin(set(df_pt_gb_filter['isoform_id']))].copy()
    for chr in set(df_pt['chr']):
        os.makedirs(outpath+"/"+chr)
        df_pt_chr=df_pt_filter[df_pt_filter['chr']==chr].copy()
        df_pt_chr.to_csv(outpath+"/"+chr+"/classification.txt",sep="\t",index=False)
        os.system("samtools view -@ {tn} -hb {bam} {chr} | samtools sort -@ {tn} -o {out}/{chr}/align.bam".format(tn=args.thread,bam=args.inputBam,chr=chr,out=outpath))
        os.system("samtools index -@ {tn} {out}/{chr}/align.bam".format(tn=args.thread,chr=chr,out=outpath))
    
    if args.thread==1:
        res=[]
        for path in pathlist:
            res=res+isosplit(path=path)
        for tbam in res:
            isoredi(tbam)
    elif args.thread>1:
        pool = multiprocessing.Pool(args.thread)
        res = pool.map(isosplit,pathlist)
        tbamlist=flatten_list(res)
        pool.close()
        pool.join()
        pbar = tqdm(total=len(tbamlist),position=0,mininterval=10)
        pbar.set_description('isoredi')
        update = lambda *args: pbar.update()
        pool = multiprocessing.Pool(args.thread)
        for tbam in tbamlist:
            pool.apply_async(isoredi, (tbam,), callback=update)
        pool.close()
        pool.join()

    else:
        infoPrint("ERROR","Threads are not valid.")
        
def mergejob(args):
    global df_tbed
    global sampleList
    df_tbed=pd.read_csv("/home/gangx/data/CQX/pacbio/gencode.v40.transcript_exon.bed",sep="\t",header=None,names=["chrom","start","end","exonid","exonorder","strand","geneid","transcriptid"])
    sampleList=args.controlSamples.split(",")+args.treatSamples.split(",")
    i=0
    for sample in sampleList:
        if i==0:
            tlist=set(getreditranscript("/home/gangx/data/CQX/pacbio/isomapit2/"+sample))
            i=1
        else:
            tlist=tlist&set(getreditranscript("/home/gangx/data/CQX/pacbio/isomapit2/"+sample))
    df_tbed_dedup=df_tbed[["chrom","geneid","transcriptid"]].drop_duplicates()
    df_tlist=df_tbed_dedup[df_tbed_dedup['transcriptid'].isin([t[:-15] for t in tlist])]
    with multiprocessing.Pool(100) as p:
        pybedtools.set_tempdir("/home/gangx/data/CQX/pacbio/isomapit2")
        r = list(tqdm(p.imap(mergeredi, df_tlist.apply(lambda x:"/".join(x),axis=1).tolist()), total=df_tlist.shape[0]))
        pybedtools.cleanup(remove_all=True)
    df_edit=pd.concat(r,axis=0)
    df_edit.to_csv(args.output+"",sep="\t",index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='ISO-MAPIT-seq (Isoform-specific \33[4mM\33[0modification \33[4mA\33[0mdded to RNA binding \33[4mP\33[0mrotein \33[4mI\33[0mnteracting \33[4mT\33[0mranscript \33[4mSeq\33[0muencing) is to identify RBP target isoform-transcripts based on adjacently editing by both hADAR2dd and rAPOBEC1. This pipeline is designed to identify RNA-editing events by hADAR2dd and rAPOBEC1 from ISO-MAPIT-seq data (TGS).\n\
        \nDetail information is available at \33[4mhttps://github.com/WangLabPKU/Mapit-seq\33[0m.\n\
        \nAuthor: Gang Xie  gangx1e@stu.pku.edu.cn\n\
        Version: 1.0  May 22, 2024', formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers()
    parser_isoedit = subparsers.add_parser('isoedit', help='Call Isoform-specific Editing Sites')
    parser_isoedit.add_argument('-v', "--genomeVersion", required=True, type=str, help='Genome Build Version')
    parser_isoedit.add_argument("-i", "--inputBam", required=True, type=str, help='Input Bam file')
    parser_isoedit.add_argument("-r", "--readAssignments", required=True, type=str, help='Input read assignments file, the isoquant output')
    parser_isoedit.add_argument('-o', "--outpath", required=True, type=str, help='Output Path for Mapit Result (default: "./isomapit")',default="./isomapit")
    parser_isoedit.add_argument('-d', "--dirname", type=str, help="Directory name in outpath path, default behavior is to use input BAM filename, without .bam extension.")
    parser_isoedit.add_argument('-t', "--thread", type=int, help='Maximum threads used for computation. (default: 10)',default=10)
    parser_isoedit.set_defaults(func=isoeditjob)

    parser_merge = subparsers.add_parser('merge', help='Merge Isoform-specific Editing Result')
    parser_merge.add_argument('-v', "--genomeVersion", required=True, type=str, help='Genome Build Version')
    parser_merge.add_argument("--controlSamples", required=True, type=str, help='Control Samples')
    parser_merge.add_argument("--treatSamples",   required=True, type=str, help='Treat Samples')
    parser_merge.add_argument('-p', "--isomapitPath", required=True, type=str, help='Output Path for Mapit Result (default: "./isomapit")',default="./isomapit")
    parser_merge.add_argument('-o', "--output", type=str, help="Directory name in outpath path, default behavior is to use input BAM filename, without .bam extension.")
    parser_merge.add_argument('-t', "--thread", type=int, help='Maximum threads used for computation. (default: 10)',default=10)
    parser_merge.set_defaults(func=mergejob)

    args = parser.parse_args()
    argsPrint(args)
    args.func(args)