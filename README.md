# MAPIT-seq Pipeline

Author: **Gang Xie**

Email: **gangx1e@stu.pku.edu.cn**

Date: May 28, 2025

Version: 1.6

----

MAPIT-seq (**M**odification **A**dded to RNA binding **P**rotein **I**nteracting **T**ranscript **Seq**uencing) is to identify RBP target transcripts based on adjacently editing by both hADAR2dd and rAPOBEC1. This pipeline is designed to identify RNA-editing events by hADAR2dd and rAPOBEC1 from MAPIT-seq data, and then find the find RBP binding regions with high confidence and *de novo* discover RBP binding motifs.

## Install

All softwares in our pipeline are available in conda. We recommend you to install these softwares by [conda](https://docs.conda.io/en/latest/miniconda.html). Also, you can use directly only if tools it depends on (see `env.yml`) have been installed. 

```
git clone https://github.com/WangLabPKU/MAPIT-seq
cd MAPIT-seq
conda env create -f env.yml
conda activate Mapit-seq
chmod +x Mapit
```

## Install MAPIT dependencies 

### REDItools2
```
cd ..
git clone https://github.com/BioinfoUNIBA/REDItools2 
cd REDItools2
pip install -r requirements.txt

```

### FLARE
```
cd ..
git clone https://github.com/YeoLab/FLARE
conda install snakemake
```


## Configuration

Reference genome, genomic annotations and known SNP sites should be provided to MAPIT-seq pipeline within `GenomeVersion.json` in the `conf` directory, e.g. `conf/GRCh38.json` and `conf/GRCm38.json`. Examples of how to get these annotations have been provided as follow.

### Download reference sequence and annotation 

- Abundant RNA (rRNA, tRNA and mtRNA) sequence (provided in `ref` fold) and create the index by BWA-MEM.

- Reference genome downloaded from [Gencode](https://www.gencodegenes.org/) and [UCSC](https://genome.ucsc.edu/index.html)

Just take "human(GRCh38)" as an example.

```
cd "your_ref_path"
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh38.p13.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.chr_patch_hapl_scaff.annotation.gff3.gz

# RepeatMasker
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz  # for mouse, just replace hg38 to mm10/mm39
gzip -d *


# ERCC spike-in
wget https://assets.thermofisher.cn/TFS-Assets/LSG/manuals/ERCC92.zip
unzip ERCC92.zip
```

We recommend extracting only sequences of autosomes, allosomes, and mitochondrial genomes (sequence name with the prefix "chr") in fasta and GTF/GFF.

### Download known variants annotation: [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/), [1000Genome](https://www.internationalgenome.org/), [EVS](http://evs.gs.washington.edu/EVS/) and [EVA](https://www.ebi.ac.uk/eva)


#### dbSNP

```
mkdir "$your_ref_path"/"$genomeVersion"_SNP
cd "$your_ref_path"/"$genomeVersion"_SNP

# human (hg38/GRCh38)
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/All_20180418.vcf.gz # download data in VCF/GATK fold; column 1 starts with "chr"

# mouse (mm10/GRCm38)
wget https://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz
```

#### Exome Variant Server or European Variation Archive

```
# human (hg38/GRCh38)
wget http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
tar -xvf ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz
for i in {1..22} X Y; do
awk '{if(substr($1, 1, 1) == "#"){print $0}else if((length($4) == 1) && (length($5) == 1)) {gsub("MT","M");{if($1 ~ "chr") print $0; else print "chr"$0 }}}' ESP6500SI-V2-SSA137.GRCh38-liftover.chr${i}.snps_indels.vcf | gzip > EVS_split_chr/chr${i}.gz
done

# mouse (mm10/GRCm38)
wget http://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_3/by_species/mus_musculus/GRCm38.p4/GCA_000001635.6_current_ids.vcf.gz
mkdir EVA_split_chr
zcat GCA_000001635.6_current_ids.vcf.gz | awk -v dir_SNP=EVA_split_chr '{if(substr($1, 1, 1) == "#"){print $0 > "EVA_header"}else if((length($4) == 1) && (length($5) == 1)) {gsub("MT","M");{if($1 ~ "chr") print $0 > dir_SNP"/"$1; else print "chr"$0 > dir_SNP"/chr"$1 }}}'
gzip EVA_split_chr/chr*

# mouse (mm39/GRCm39)
## not recommend, because of lack of mm39/GRCm39 data in dbSNP NCBI. But you could transfer from mm10/GRCm38 data by liftOver.
wget http://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_3/by_species/mus_musculus/GRCm39/GCA_000001635.9_current_ids.vcf.gz
```

#### 1000 Genomes Project

```
# human  
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz
mkdir 1000genomes_split_chr
zcat ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz | awk -v dir_SNP=1000genomes_split_chr '{if(substr($1, 1, 1) == "#"){print $0 > "1000genomes_header"}else if((length($4) == 1) && (length($5) == 1)) {gsub("MT","M");{if($1 ~ "chr") print $0 > dir_SNP"/"$1; else print "chr"$0 > dir_SNP"/chr"$1 }}}'
gzip 1000genomes_split_chr/chr*
```

There is no any known SNP data for mouse in 1000 Genomes Project, so you could create void file in `void_split_chr` fold.

```
# mouse
chroms=($(grep '>' $genome_fasta |sed 's/>//' | awk '{print $1}' | grep 'chr' ))
for chr in ${chroms[@]}; do touch void_split_chr/${chr}; done
gzip void_split_chr/chr*
```

**Attention! All main chromosomes' names in "XXX_split_chr" fold should be started with "chr". If not, please add the prefix "chr".**

#### Create configuration file

After genome sequence `.fasta` file, genome anntation `.gtf` and `.gff3` files, and RepeatMasker annotation downloaded and known SNP data splitted by chromosomes, you could run the command below to complete configuration. These are including creating sequence index, extract annotation of gene elements and creating dbSNP `.vcf` file splitted by chromosomes. 

```
Mapit config --genomeVersion GRCh38 \
             --genomeFasta "full_path"/GRCh38.p13.genome.fa \
             --ERCC "full_path"/"ERCC.fa" \
             --species human \
             --outdir "full_path" \
             --genomeAnno "full_path"/gencode.v40.chr_patch_hapl_scaff.annotation.gff3 \
             --rmsk "full_path"/rmsk.txt \
             --dbSNP "full_path"/GRCh38_SNP/All_20180418.vcf.gz \
             --1000genomesDir "full_path"/GRCh38_SNP/1000genomes_split_chr \
             --EVSEVADir "full_path"/GRCh38_SNP/EVS_split_chr \
             --Reditools "full_path"/Directory_of_RediTools2.0_software \
             --FLARE "full_path"/Directory_of_FLARE_software
```

- Options:

  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `-s {human,mouse}, --species {human,mouse}`
                        Species

  `-f GENOMEFASTA, --genomeFasta GENOMEFASTA`
                        Genome Sequence Fasta File

  `-E ERCC, --ERCC ERCC`  ERCC Spikein Fasta File

  `-a GENOMEANNO, --genomeAnno GENOMEANNO`
                        Genome Annotation File GFF3 Format

  `-r RMSK, --rmsk RMSK`  RepeatMask Annotation Downloaded from UCSC

  `-o OUTPATH, --outpath OUTPATH`
                        Annotation File Output Path for Mapit

  `--dbSNP DBSNP`         VCF File of NCBI dbSNP

  `--1000Genomes 1000GENOMES`
                        Directory of VCF Files of 1000Genomes Splited by Chromosomes

  `--EVSEVA EVSEVA`       Directory of VCF Files of EVS/EVA Splited by Chromosomes

  `--hisat2Index HISAT2INDEX`
                        Path to "hisat2-build" Index, if You've Built Beforehand
  `--Reditools REDITOOLS`
                        Directory of RediTools2.0 software
  `--FLARE FLARE`         Directory of FLARE software
  `--overwrite`           Overwrite Config File or not.


    



## Usage

### Quick start

### Mapping

Mapit-seq pipeline can break down into five main steps, while abundant RNA filtering and reads mapping are integrated into one part during operation:

```
Mapit mapping -v GRCh38 --fq R1.fq.gz  --fq2 R2.fq.gz  --rna-strandness FR -n EN -r 1  -o Mapit_result -t 40
```

- Options:

  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `--fq FQ`               RNA-seq Paired-end Fastq R1 File or Single-end Fastq File

  `--fq2 FQ2`             RNA-seq Paired-end Fastq R2 File

  `--rna-strandness {F,R,FR,RF}`
                        The strand-specific information used for HISAT2 mapping. For single-end reads, use F or R. For paired-end reads, use either FR or RF. Detailed descriptions of
                        this option was available in HISAT2 manual (https://ccb.jhu.edu/software/hisat2/manual.shtml)

  `-n SAMPLENAME, --sampleName SAMPLENAME`
                        Directory name in output directory and the first prefix of file name for the mapping results

  `-r REPLICATE, --replicate REPLICATE`
                        The 2nd prefix of file name for the mapping results

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")

  `-t THREAD, --thread THREAD`
                        Maximum threads used for computation. (default: 10)

  `--ERCC`                When ERCC spikein used, add this parameter. (default: False)


This create 4 sub-directories in the output directory. And some mapping `*.bam` files will create.

### Fine tune

```
Mapit finetuning -v GRCh38 -n EN -r 1 -o "output_path"/Mapit_result
```

- Options:
        
  `-h, --help`            show this help message and exit
  
  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `-n SAMPLENAME, --sampleName SAMPLENAME`
                        Directory name in output directory and the first prefix of file name for the mapping results

  `-r REPLICATE, --replicate REPLICATE`
                        The 2nd prefix of file name for the mapping results

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")


### Call editing events

This step allows for the simultaneous analysis of multiple samples. 

```
Mapit callediting -v GRCh38 --sampleList EN,EM,NN,NM -o "output_path"/Mapit_result --prefix RBP -e Both [-t THREAD]
```

- Options:

  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `--sampleList SAMPLELIST`
                        Directory name for samples inputing in output directory. Use "," to combine.

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")

  `--prefix PREFIX`       Output Directory prefix in Mapit Result Dir

  `-e {ADAR,APOBEC,Both}, --enzyme {ADAR,APOBEC,Both}`
                        RNA-editing enzymes used. ADAR,APOBEC: MAPIT-seq (default: Both); ADAR: HyperTRIBE or TRIBE; APOBEC: STAMP

  `-t THREAD, --thread THREAD`
                        Maximum threads used for computation. (default: 10)


### Call RBP targets

The step will output the "Differential editing analysis" results in the "${outpath}/6-Edit_calling/${prefix}${treatSampleName}/" directory.

```
Mapit calltargets -v GRCh38 -i "output_path"/Mapit_result/6-Edit_calling/RBP/RBP_Edit_GE_RPE_DATA.tsv --treatName EM --controlName EN -o "output_path"/Mapit_result -p 0.05 
```

- Options:
    
  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `-i INPUTEDIT, --inputEdit INPUTEDIT`
                        Input RNA Edit Table File

  `-l {transcript,gene,Both}, --level {transcript,gene,Both}`
                        Perform Differential Editing Analysis in transcript or gene(including intron) level

  `--treatName TREATNAME`
                        Directory name for treatment samples in output directory

  `--controlName CONTROLNAME`
                        Directory name for control samples in output directory

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")

  `-b BINSIZE, --binSize BINSIZE`
                        The length of bins that are split from transcripts. (default: 50)

  `-c COVERAGE, --coverage COVERAGE`
                        Minimun coverage for valid editing sites. (default: 10)

  `--supply-zero`         Supply zero to non-editing bins. (default: False)

  `-p PVALUE, --pvalue PVALUE`
                        Maximum p-value for Poisson filter of editing sites. (default: 0.1)

  `--dropTreatRep DROPTREATREP`
                        Replicate id for treatment samples in output directory. (e.g. 1,4,5)

  `--dropControlRep DROPCONTROLREP`
                        Replicate id for control samples in output directory. (e.g. 1,4,5)

  `--singleCell`          When a replicate represents a cell, add this parameter. (default: False)

### Identifying high-confidence editing clusters and RBP binding motifs

Prepare files for SAILOR workflow.

```
Mapit prepare -v $genomeVersion -n ${samplename} -r ${replicate} -o ${outpath}
```

- Options:

  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `-n SAMPLENAME, --sampleName SAMPLENAME`
                        Directory name in output directory and the first prefix of file name for the mapping results

  `-r REPLICATE, --replicate REPLICATE`
                        The 2nd prefix of file name for the mapping results

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")

Run SAILOR workflow

```
Mapit SAILOR -v $genomeVersion -n ${samplename} -r ${replicate} -c ${coverage} -o ${outpath} -t ${threads}
```

- Options:

  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `-n SAMPLENAME, --sampleName SAMPLENAME`
                        Directory name in output directory and the first prefix of file name for the mapping results

  `-r REPLICATE, --replicate REPLICATE`
                        The 2nd prefix of file name for the mapping results

  `-c COVERAGE, --coverage COVERAGE`
                        Minimun coverage for valid editing sites. (default: 10)

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")

  `-t THREAD, --thread THREAD`
                        Maximum threads used for computation. (default: 10)

Run FLARE workflow to identify edit clusters

```
Mapit FLARE -v $genomeVersion -e {AG,CT} -n ${samplename} -r ${replicate} --regions ${regions} -o ${outpath} -t ${threads}
```

- Options:

  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `-n SAMPLENAME, --sampleName SAMPLENAME`
                        Directory name in output directory and the first prefix of file name for the mapping results

  `-r REPLICATE, --replicate REPLICATE`
                        The 2nd prefix of file name for the mapping results

  `--regions REGIONS`     FLARE configuration directory of files for regions of the genome

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")

  `-t THREAD, --thread THREAD`
                        Maximum threads used for computation. (default: 10)
                        
  `-e {AG,CT}, --edittype {AG,CT}`
                        Editing types

### Identify MAPIT high-confidence edit clusters

```
Mapit hc_cluster -v $genomeVersion -n ${samplename} -o ${outpath} -s ${sloplength}
```

- Options:
  
  `-h, --help`            show this help message and exit

  `-v GENOMEVERSION, --genomeVersion GENOMEVERSION`
                        Genome Build Version

  `-n SAMPLENAME, --sampleName SAMPLENAME`
                        Directory name in output directory and the first prefix of file name for the mapping results

  `-o OUTPATH, --outpath OUTPATH`
                        Output Path for Mapit Result (default: "./Mapit_result")

  `-s SLOPLENGTH, --sloplength SLOPLENGTH`
                        Length of High-confidence clusters expanded for up- and down-stream sides

                        
## License

Copyright (C) 2025 WangLabPKU.