#!/usr/bin/bash
set -e
MYDIR=`dirname $0`

Creating_conf(){
    local genome_build_version=$1
    local organism=$2
    local annotation_path=$3
    local genome_fasta=$4
    local gencode_gff3=$5
    local hisat2_index=$6
    local dbSNP_all=$7
    local rmsk=$8
    local ERCC=$9
    
    local hisat2_index_s=true
    local bwa_index=true

    test -d ${annotation_path}/AbundantRNA || mkdir -p ${annotation_path}/AbundantRNA
    test -f ${annotation_path}/AbundantRNA/${organism}_rRNA_tRNA_mtRNA.fa || cp ${MYDIR}/../data/${organism}_rRNA_tRNA_mtRNA.fa ${annotation_path}/AbundantRNA/
    abundantRNA_index_bwa_mem=${annotation_path}/AbundantRNA/${organism}_rRNA_tRNA_mtRNA.fa
    test -f ${abundantRNA_index_bwa_mem}.sa || bwa index ${abundantRNA_index_bwa_mem}
    
    for((i=1;i<=8;i++));do if [ ! -f ${hisat2_index}.${i}.ht2 ];then hisat2_index_s="false"; fi; done
    for suffix in amb ann bwt pac sa ;do if [ ! -f ${genome_fasta}.${suffix} ];then bwa_index="false"; fi; done

    if [ "$bwa_index" == "false" ]; then 
    bwa index ${genome_fasta} & 
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCreating BWA index..." > /dev/tty
    fi
    if [ "$hisat2_index_s" == "false" ]; then 
    hisat2-build ${genome_fasta} ${hisat2_index} & 
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCreating HISAT2 index..." > /dev/tty
    fi
    if [ ! -f ${genome_fasta%.*}.dict ]; then 
    gatk CreateSequenceDictionary -R ${genome_fasta} & 
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCreating GATK index..." > /dev/tty
    fi
    if [ ! -f ${genome_fasta}.fai ]; then 
    samtools faidx ${genome_fasta} & 
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCreating FAI index..." > /dev/tty
    fi
    if [ ! -f ${gencode_gff3%.*}_all_spsites.bed ]; then 
    gffread ${gencode_gff3} -T | hisat2_extract_splice_sites.py - > ${gencode_gff3%.*}_all_spsites.bed
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCreating HISAT2 file of splice sites index..." > /dev/tty
    fi
    hisat2-build ${ERCC} ${ERCC}

    [ -e ${annotation_path}/fd ] || mkfifo ${annotation_path}/fd
    exec 3<>${annotation_path}/fd
    rm -rf ${annotation_path}/fd
    for ((i=1;i<=10;i++)); do echo >&3; done

    test -f ${gencode_gff3%.*}_exon.bed || awk -v OFS="\t" '{split($9,a,";");split(a[10],d,"=");split(a[3],e,"=");split(a[5],b,"=");split(a[7],c,"=");if($3=="exon") print $1,$4-1,$5,d[2]":"e[2],c[2]":"b[2],$7 }'  ${gencode_gff3} | bedtools sort -i - | uniq >  ${gencode_gff3%.*}_exon.bed
    if [[ ! -f ${gencode_gff3%.*}_exon_mRNA_merge.bed ]]; then
        geneList=($(awk '{if($5 ~ "protein_coding:protein_coding") {split($4,a,":");print a[2] }}' ${gencode_gff3%.*}_exon.bed | sort | uniq ))
        mkdir -p ${annotation_path}/tmp
        for gene in ${geneList[@]}; do 
        read -u3
        {
            awk -v gene=$gene '{if(($4 ~ gene)&&($5 ~ "protein_coding:protein_coding")) print $0}' ${gencode_gff3%.*}_exon.bed | bedtools merge -i - -c 4,5,6 -o distinct -delim ";"> ${annotation_path}/tmp/${gene}.bed
            if [ "`cut -f 6 ${annotation_path}/tmp/${gene}.bed | uniq `" == "-" ]; then 
            awk '{print $0"\t"$3-$2}' ${annotation_path}/tmp/${gene}.bed | sort -n -k 2 -r | awk '{sum+=$7; print $0"\t"sum}' > ${annotation_path}/tmp/${gene}_cumsum.bed
            else
            awk '{print $0"\t"$3-$2}' ${annotation_path}/tmp/${gene}.bed | sort -n -k 2 | awk '{sum+=$7; print $0"\t"sum}' > ${annotation_path}/tmp/${gene}_cumsum.bed
            fi
            rm ${annotation_path}/tmp/${gene}.bed 
            echo >&3
            }&
        done
        wait
        find ${annotation_path}/tmp/ -name *_cumsum* | xargs cat | bedtools sort -i  > ${gencode_gff3%.*}_exon_mRNA_merge.bed
        rm -rf ${annotation_path}/tmp
    fi
    if [[ ! -f ${gencode_gff3%.*}_exon_ncRNA_merge.bed ]]; then
        geneList=($(grep -v "protein_coding" ${gencode_gff3%.*}_exon.bed  | awk '{split($4,a,":");print a[2]}'  | sort | uniq ))
        mkdir -p ${annotation_path}/tmp
        for gene in ${geneList[@]}; do
        read -u3;     
        {
            awk -v gene=$gene '{if(($4 ~ gene)&&($5 !~ "intron" )&&($5 !~ "protein_coding" )) print $0}' ${gencode_gff3%.*}_exon.bed  | bedtools merge -i - -c 4,5,6 -o distinct -delim ";" > ${annotation_path}/tmp/${gene}.bed
            if [ "`cut -f 6 ${annotation_path}/tmp/${gene}.bed | uniq `" == "-" ]; then 
            awk '{print $0"\t"$3-$2}' ${annotation_path}/tmp/${gene}.bed | sort -n -k 2 -r | awk '{sum+=$7; print $0"\t"sum}' > ${annotation_path}/tmp/${gene}_cumsum.bed
            else
            awk '{print $0"\t"$3-$2}' ${annotation_path}/tmp/${gene}.bed | sort -n -k 2 | awk '{sum+=$7; print $0"\t"sum}' > ${annotation_path}/tmp/${gene}_cumsum.bed
            fi
            rm ${annotation_path}/tmp/${gene}.bed 
            echo >&3
        }& 
        done
        wait
        find ${annotation_path}/tmp/ -name *_cumsum* | xargs cat | bedtools sort -i > ${gencode_gff3%.*}_exon_ncRNA_merge.bed
        rm -rf ${annotation_path}/tmp
    fi
    if [[ ! -f ${gencode_gff3%.*}_exon_merge_seq.bed ]]; then
        cat ${gencode_gff3%.*}_exon_mRNA_merge.bed ${gencode_gff3%.*}_exon_ncRNA_merge.bed > ${gencode_gff3%.*}_exon_merge.bed
        bedtools getfasta -bed ${gencode_gff3%.*}_exon_merge.bed -fi ${genome_fasta} -s -bedOut > ${gencode_gff3%.*}_exon_merge_seq.bed
        rm ${gencode_gff3%.*}_exon_mRNA_merge.bed ${gencode_gff3%.*}_exon_ncRNA_merge.bed
    fi 
    # dbSNP
    
    test -d ${annotation_path}/${genome_build_version}_SNP/dbSNP_split_chr || mkdir -p ${annotation_path}/${genome_build_version}_SNP/dbSNP_split_chr
    test -f ${annotation_path}/${genome_build_version}_SNP/dbSNP_all_SNP.vcf || zcat ${dbSNP_all} | awk '{if(substr($1, 1, 1) == "#"){print $0}else if((length($4) == 1) && (length($5) == 1)) {gsub("MT","M");{if($1 ~ "chr") print $0; else print "chr"$0 }}}' > ${annotation_path}/${genome_build_version}_SNP/dbSNP_all_SNP.vcf
    test -f ${annotation_path}/${genome_build_version}_SNP/dbSNP_all_SNP.vcf.idx || gatk IndexFeatureFile -I ${annotation_path}/${genome_build_version}_SNP/dbSNP_all_SNP.vcf
    
    chroms=($(cut -f 1 ${gencode_gff3%.*}_all_spsites.bed | grep chr | uniq))
    flag=0
    for chr in ${chroms[@]}; do test -f ${annotation_path}/${genome_build_version}_SNP/dbSNP_split_chr/${chr}.gz || flag=1;done
    if [ ${flag} == 0 ]; then
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tSplited dbSNP files exist." > /dev/tty
    else
    {
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCreating dbSNP files..." > /dev/tty
    awk -v dir_SNP=${annotation_path}/${genome_build_version}_SNP/dbSNP_split_chr '{if(substr($1, 1, 1) == "#"){print $0 > dir_SNP"/../dbSNP_header"}else {gsub("MT","M"); print $0 > dir_SNP"/"$1}}' ${annotation_path}/${genome_build_version}_SNP/dbSNP_all_SNP.vcf
    pigz ${annotation_path}/${genome_build_version}_SNP/dbSNP_split_chr/chr*
    }
    fi

    test -f ${annotation_path}/UCSC_RepeatMask_All_repetitive.bed  || awk '{print $6"\t"$7"\t"$8"\t"$11"\t"$12"\t"$10}' ${rmsk} |  awk -v nA=FAM_FRAM_FLAM_A_FLAM_C -v OFS="\t" '{gsub("?","",$5);if(($5=="SINE")&&(index(nA,$4)>0 || $4 ~ "Alu")) print $1,$2,$3,$4,"SINE_Alu",$6; else print $0}' >  ${annotation_path}/UCSC_RepeatMask_All_repetitive.bed 

    wait
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tALL done! Checking..." > /dev/tty
}

rRNA_deplete_HISAT2_BWA_mapping(){
    local outname=$1
    local repID=$2
    local layout=$3
    local fq1_path=$4
    local fq2_path=$5
    local output_path=$6
    local thread=$7
    local stranded=$8
    local abundantRNA_index_bwa_mem=$9
    local annotation_splice_sites=${10}
    local genome_index_hisat2=${11}
    local genome_index_bwa_mem=${12}
    local genome_fasta=${13}
    local ERCC_spikein=${14}
    local ERCC_index_hisat2=${15}

    local rRNA_map="${output_path}/0-Remove_rRNA/${outname}"
    local HISAT_map="${output_path}/1-HISAT_map/${outname}"
    local bwa_map="${output_path}/2-BWA_map/${outname}"
    local combine_path="${output_path}/3-Combine_bam/${outname}"

    #echo outname":"$outname,repID":"$repID,layout":"$layout,fq1_path":"$fq1_path,fq2_path":"$fq2_path,output_path":"$output_path,thread":"$thread,stranded":"$stranded
    if [[ -z "$MYDIR" ]] ; then exit 1;fi
    #echo $abundantRNA_index_bwa_mem,$annotation_splice_sites,$genome_index_hisat2,$genome_index_bwa_mem
    test -d ${rRNA_map} || mkdir -p ${rRNA_map}
    test -d ${HISAT_map} || mkdir -p ${HISAT_map}
    test -d ${bwa_map} || mkdir -p ${bwa_map}
    test -d ${combine_path} ||mkdir -p ${combine_path}
    local rep_name=${outname}_${repID}

    if [ "$layout" == "paired" ];then
        ## remove AbundantRNA(rRNA tRNA mtRNA)
        bwa  mem -t ${thread} ${abundantRNA_index_bwa_mem}  ${fq1_path}  ${fq2_path} | samtools view -@ ${thread} -bh -f 4 - | samtools sort -@ ${thread} -n - -o ${rRNA_map}/${rep_name}-rRNA_unmapped_sort.bam
        samtools fastq -@ ${thread} -1 ${rRNA_map}/${rep_name}_R1.fastq.gz -2 ${rRNA_map}/${rep_name}_R2.fastq.gz -s ${rRNA_map}/${rep_name}_singleton.fq ${rRNA_map}/${rep_name}-rRNA_unmapped_sort.bam
        rm ${rRNA_map}/${rep_name}-rRNA_unmapped_sort.bam
        fq1=${rRNA_map}/${rep_name}_R1.fastq.gz
        fq2=${rRNA_map}/${rep_name}_R2.fastq.gz
        echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tremove abundantRNA done! " > /dev/tty

        ### 1. HISAT2 2 mismatches mapping
        hisat2  --rna-strandness ${stranded} --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${annotation_splice_sites} --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t \
                --rg-id ${rep_name} --rg SM:${rep_name} --rg SM:${rep_name} --rg PL:ILLUMINA --rg PU:NovaSeq6000 \
                -p ${thread} -x   ${genome_index_hisat2} -1  $fq1  -2 $fq2  --un-conc-gz ${HISAT_map}/${rep_name}_un_conc_%.fastq.gz -S ${HISAT_map}/${rep_name}_HISAT2_mapped.sam &> ${HISAT_map}/${rep_name}_HISAT2_mapped.summary
        samtools view -@ ${thread} -h -F 4 ${HISAT_map}/${rep_name}_HISAT2_mapped.sam | \
            awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag"  }}' | \
            awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) {print $0 > "'${HISAT_map}/${rep_name}_unique_mismatch2.sam'"}
                    else{print $0 > "'${HISAT_map}/${rep_name}_unique_mismatch_over2.sam'"} } else {print $0 " not have XM tag"} } }' > ${HISAT_map}/${rep_name}_HISAT2.header
        samtools view -@ ${thread} -h -f 4 ${HISAT_map}/${rep_name}_HISAT2_mapped.sam | cat - ${HISAT_map}/${rep_name}_unique_mismatch_over2.sam | samtools view -@ ${thread} -bh  | samtools sort -@ ${thread} -n  -o ${HISAT_map}/${rep_name}_HISAT2_mapped-unmapped_sorted.bam
        samtools fastq -@ ${thread} -1 ${HISAT_map}/${rep_name}_unmapped_1.fastq.gz -2 ${HISAT_map}/${rep_name}_unmapped_2.fastq.gz -s ${HISAT_map}/${rep_name}_unmapped_singleton.fq  ${HISAT_map}/${rep_name}_HISAT2_mapped-unmapped_sorted.bam
        echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tHISAT2 mapping done! " > /dev/tty

        ### 2. ERCC / BWA 0.04 per bp mismatches mapping
        test -d $bwa_map||mkdir -p $bwa_map
        if [ "$ERCC_spikein" == "True" ];then
        hisat2 -q -p ${thread} --no-spliced-alignment --end-to-end --no-mixed --secondary --rdg 10000,10000 --rfg 10000,10000 \
              -x ${ERCC_index_hisat2} -1  ${HISAT_map}/${rep_name}_unmapped_1.fastq.gz  -2 ${HISAT_map}/${rep_name}_unmapped_2.fastq.gz \
              --un-conc-gz ${bwa_map}/${rep_name}_ERCC_un_conc_%.fastq.gz -S ${bwa_map}/${rep_name}_ERCC.sam &> ${HISAT_map}/${rep_name}_ERCC.summary
        echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tMapped to ERCC done! " > /dev/tty
        samtools view -hb -F 4 -@ 20 ${bwa_map}/${rep_name}_ERCC.sam | samtools sort -@ 20 -o ${bwa_map}/${rep_name}_ERCC.bam
        rm ${bwa_map}/${rep_name}_ERCC.sam
        samtools view ${bwa_map}/${rep_name}_ERCC.bam  | cut -f 3 | uniq -c > ${bwa_map}/${rep_name}_ERCC.count.txt
        bwa mem -t ${thread}  -A 1 -B 4 -R "@RG\tID:${rep_name}\tSM:${rep_name}\tSM:${rep_name}\tPL:ILLUMINA\tPU:NovaSeq6000" \
                ${genome_index_bwa_mem}  ${bwa_map}/${rep_name}_ERCC_un_conc_1.fastq.gz  ${bwa_map}/${rep_name}_ERCC_un_conc_2.fastq.gz > ${bwa_map}/${rep_name}_bwa_mapped.sam 
        else
        bwa mem -t ${thread}  -A 1 -B 4 -R "@RG\tID:${rep_name}\tSM:${rep_name}\tSM:${rep_name}\tPL:ILLUMINA\tPU:NovaSeq6000" \
                ${genome_index_bwa_mem}  ${HISAT_map}/${rep_name}_unmapped_1.fastq.gz  ${HISAT_map}/${rep_name}_unmapped_2.fastq.gz > ${bwa_map}/${rep_name}_bwa_mapped.sam 
        fi
    elif [ "$layout" == "single" ];then
	    ## remove AbundantRNA(rRNA tRNA mtRNA)
	    bwa  mem -t ${thread} ${abundantRNA_index_bwa_mem}  $fq1_path | samtools view -@ ${thread} -bh -f 4 - | samtools sort -@ ${thread} -n - -o ${rRNA_map}/${rep_name}-rRNA_unmapped_sort.bam
        samtools fastq -@ ${thread}  ${rRNA_map}/${rep_name}-rRNA_unmapped_sort.bam | pigz - > ${rRNA_map}/${rep_name}.fastq.gz
        rm ${rRNA_map}/${rep_name}-rRNA_unmapped_sort.bam
        fq0=${rRNA_map}/${rep_name}.fastq.gz
        echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tremove abundantRNA done! " > /dev/tty

        ### 1. HISAT2 2 mismatches mapping 
        hisat2 --rna-strandness ${stranded} --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${annotation_splice_sites} --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t \
               --rg-id ${rep_name} --rg SM:${rep_name} --rg SM:${rep_name} --rg PL:ILLUMINA --rg PU:NovaSeq6000 \
               -p ${thread} -x ${genome_index_hisat2} -U $fq0 -S ${HISAT_map}/${rep_name}_HISAT2_mapped.sam &> ${HISAT_map}/${rep_name}_HISAT2_mapped.summary
        samtools view -@ ${thread} -h -F 4 ${HISAT_map}/${rep_name}_HISAT2_mapped.sam | \
            awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag"  }}' | \
            awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) {print $0 >> "'${HISAT_map}/${rep_name}_unique_mismatch2.sam'"}
                    else{print $0 >> "'${HISAT_map}/${rep_name}_unique_mismatch_over2.sam'"} } else {print $0 " not have XM tag"} } }' > ${HISAT_map}/${rep_name}_HISAT2.header
        samtools view -@ ${thread} -h -f 4 ${HISAT_map}/${rep_name}_HISAT2_mapped.sam | cat - ${HISAT_map}/${rep_name}_unique_mismatch_over2.sam | samtools view -@ ${thread} -bh  | samtools sort -@ ${thread} -n  -o ${HISAT_map}/${rep_name}_HISAT2_unmapped_sort.bam
        samtools fastq -@ ${thread}  ${HISAT_map}/${rep_name}_HISAT2_unmapped_sort.bam | pigz -  > ${HISAT_map}/${rep_name}_HISAT2_unmapped.fastq.gz
        echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tHISAT2 mapping done! " > /dev/tty

        ### 2. ERCC / BWA 0.04 per bp mismatches mapping
        if [ "$ERCC_spikein" == "True" ];then
        hisat2 -q -p ${thread} --no-spliced-alignment --end-to-end --no-mixed --secondary --rdg 10000,10000 --rfg 10000,10000 \
              -x ${ERCC_index_hisat2} -U ${HISAT_map}/${rep_name}_HISAT2_unmapped.fastq.gz \
              --un-gz ${bwa_map}/${rep_name}_ERCC_unmapped.fastq.gz -S ${bwa_map}/${rep_name}_ERCC.sam &> ${HISAT_map}/${rep_name}_ERCC.summary
        echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tMapped to ERCC done! " > /dev/tty
        samtools view -hb -F 4 -@ 20 ${bwa_map}/${rep_name}_ERCC.sam | samtools sort -@ 20 -o ${bwa_map}/${rep_name}_ERCC.bam
        rm ${bwa_map}/${rep_name}_ERCC.sam
        samtools view ${bwa_map}/${rep_name}_ERCC.bam | cut -f 3 | uniq -c > ${bwa_map}/${rep_name}_ERCC.count.txt
        bwa mem -t ${thread}  -A 1 -B 4 -R "@RG\tID:${rep_name}\tSM:${rep_name}\tSM:${rep_name}\tPL:ILLUMINA\tPU:NovaSeq6000" \
                ${genome_index_bwa_mem}  ${bwa_map}/${rep_name}_ERCC_unmapped.fastq.gz > ${bwa_map}/${rep_name}_bwa_mapped.sam 
        else
        bwa mem -t ${thread} -R "@RG\tID:${rep_name}\tSM:${rep_name}\tSM:${rep_name}\tPL:ILLUMINA\tPU:NovaSeq6000" \
                ${genome_index_bwa_mem}  ${HISAT_map}/${rep_name}_HISAT2_unmapped.fastq.gz > ${bwa_map}/${rep_name}_bwa_mapped.sam
        fi
    fi

    python ${MYDIR}/bwa_unique_mismatch6.py ${bwa_map}/${rep_name}_bwa_mapped.sam ${bwa_map}/${rep_name}_bwa_unique_mis6_mapq0.sam 
    samtools  view -@ ${thread} -bT ${genome_fasta} ${bwa_map}/${rep_name}_bwa_unique_mis6_mapq0.sam | samtools  sort -@ - -o ${bwa_map}/${rep_name}_unmapped.sort.bam
    samtools  view -@ ${thread} -H ${bwa_map}/${rep_name}_unmapped.sort.bam > ${bwa_map}/${rep_name}_bwa.header

    rm ${bwa_map}/${rep_name}_bwa_mapped.sam ${HISAT_map}/${rep_name}_HISAT2_mapped.sam ${bwa_map}/${rep_name}_bwa_unique_mis6_mapq0.sam
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tBWA mapping done! " > /dev/tty
    sleep 1
    
    cat ${bwa_map}/${rep_name}_bwa.header ${HISAT_map}/${rep_name}_unique_mismatch2.sam | samtools view -@ ${thread} -bT ${genome_fasta} - | samtools sort -@ ${thread} - -o ${combine_path}/${rep_name}_accepted_hits.sort.bam
    samtools merge -@ ${thread} -f ${combine_path}/${rep_name}_combine.bam ${combine_path}/${rep_name}_accepted_hits.sort.bam ${bwa_map}/${rep_name}_unmapped.sort.bam
    samtools index -@ ${thread} ${combine_path}/${rep_name}_combine.bam
    rm ${HISAT_map}/${rep_name}_unique_mismatch2.sam 
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tHISAT2-BWA combining done! " > /dev/tty
}

Filtering_SNV(){
    local input_file=$1
    local output_file=$2
    local SNP_path=$3
    local thread=$4
    local work_path=`dirname $input_file`
    local input_file_basename=`basename $input_file`
    test -d ${work_path}/tmp_deknownSNP/chr  ||mkdir -p ${work_path}/tmp_deknownSNP/chr 
    test -d ${work_path}/tmp_deknownSNP/tmp_result  ||mkdir -p ${work_path}/tmp_deknownSNP/tmp_result
    [ -e ${work_path}/tmp_deknownSNP/fd1 ] || mkfifo ${work_path}/tmp_deknownSNP/fd1
    exec 3<>${work_path}/tmp_deknownSNP/fd1
    rm -rf ${work_path}/tmp_deknownSNP/fd1

    for ((i=1;i<=${thread};i++)); do echo >&3; done

    awk '$0 !~/^#/ && $0 !~/^$/{print  > "'${work_path}'/tmp_deknownSNP/chr/"$1}' ${input_file}
    chr_list=($(ls ${work_path}/tmp_deknownSNP/chr|sort -k1.4n))
    for chr in ${chr_list[@]}; do
    read -u3
    {
        echo $chr| grep -q "_"  &&continue
        if [ -e "${SNP_path}/${chr}.gz" ];then
        awk 'FILENAME==ARGV[1]{array_tmp[$1":"$2]=1}FILENAME==ARGV[2]{if (! array_tmp[$1":"$2] ){print}}'  <(zcat ${SNP_path}/${chr}.gz) ${work_path}/tmp_deknownSNP/chr/${chr} >${work_path}/tmp_deknownSNP/tmp_result/${input_file_basename}_${chr} 
        else
        echo -e "\033[33m[WARN]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\t$chr in $input_file do not exists in $SNP_path" > /dev/tty
        cp ${work_path}/tmp_deknownSNP/chr/${chr} ${work_path}/tmp_deknownSNP/tmp_result/${input_file_basename}_${chr} 
        fi
        echo >&3
    }&
    done
    wait
    cat ${work_path}/tmp_deknownSNP/tmp_result/* >$output_file
    rm -r ${work_path}/tmp_deknownSNP
}

Calling_SNV_filtering_out_knownSNP(){
    local output_path=$1
    local prefix=$2
    local sampleInput=$3
    local thread=$4
    local genome_fasta=$5
    local intervals=$6
    local SNP_dbSNP_divided_by_chromosome=$7
    local SNP_1000Genome_divided_by_chromosome=$8
    local SNP_EVS_divided_by_chromosome=$9

    local var_call_path="${output_path}/4-Var_calling"
    local vcf_filter_path="${output_path}/5-Var_filter"
    local edit_call_path="${output_path}/6-Edit_calling"
    sample_name=`basename ${prefix}`
    local vcf_inter_name=HaplotypeCaller_BQdefault_MAPQ0
    local vcf_split_chr_path="${prefix}"
    local vcf_combine_prefix="${prefix}_${vcf_inter_name}"
    local vcf_deSNP_prefix="${vcf_filter_path}/${sample_name}_${vcf_inter_name}_deAllSNP_"

    test -d ${vcf_split_chr_path} || mkdir -p ${vcf_split_chr_path}

    [ -e ${prefix}/fd ] || mkfifo ${prefix}/fd
    exec 3<>${prefix}/fd
    rm -rf ${prefix}/fd

    processing=$((thread/4+1))
    for ((i=1;i<=$processing;i++)); do echo >&3; done

    chroms=($(cut -f 1 $intervals ))
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tStart to HaplotypeCaller splited by chromosomes." > /dev/tty
    for chr in ${chroms[@]}; do 
    read -u3
    {
        gatk HaplotypeCaller -R ${genome_fasta} -I ${sampleInput//,/" -I "} --minimum-mapping-quality 0 --dont-use-soft-clipped-bases true  -stand-call-conf 0 --force-active true --all-site-pls --heterozygosity 1 --sample-ploidy 4 -O ${vcf_split_chr_path}/${chr}_${vcf_inter_name}.vcf --intervals ${chr} 
        echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\t$chr HaplotypeCaller finished." > /dev/tty
        echo >&3
    }&
    done
    wait
    
    local merge_vcfs=""
    for chr in ${chroms[@]}; do merge_vcfs=${merge_vcfs}" -I ${vcf_split_chr_path}/${chr}_${vcf_inter_name}.vcf"; done
    gatk MergeVcfs ${merge_vcfs} -O ${vcf_combine_prefix}.vcf.gz 
    gatk SelectVariants -select-type SNP -R ${genome_fasta} -V ${vcf_combine_prefix}.vcf.gz -O ${vcf_combine_prefix}_SNP.vcf
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tSNVs were selected from all variants." > /dev/tty

    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tStart to filter out known SNPs." > /dev/tty
    test -d ${vcf_filter_path} || mkdir -p ${vcf_filter_path}
    Filtering_SNV ${vcf_combine_prefix}_SNP.vcf ${vcf_deSNP_prefix}dbSNP.vcf ${SNP_dbSNP_divided_by_chromosome} ${thread}
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tKnown SNPs from dbSNP were filtered out." > /dev/tty
    Filtering_SNV ${vcf_deSNP_prefix}dbSNP.vcf ${vcf_deSNP_prefix}dbSNP_1000genomes.vcf ${SNP_1000Genome_divided_by_chromosome} ${thread}
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tKnown SNPs from 1000Genomes were filtered out." > /dev/tty
    Filtering_SNV ${vcf_deSNP_prefix}dbSNP_1000genomes.vcf ${vcf_deSNP_prefix}dbSNP_1000genomes_EVS.vcf ${SNP_EVS_divided_by_chromosome} ${thread}
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tKnown SNPs from EVS or EVA were filtered out." > /dev/tty
    grep "#CHROM" ${vcf_combine_prefix}_SNP.vcf | head -n 1 | cat - ${vcf_deSNP_prefix}dbSNP_1000genomes_EVS.vcf > ${vcf_deSNP_prefix}dbSNP_1000genomes_EVS_Header.vcf
}

GE_Annotate(){
    local input_file=$1
    local gencode_gff3=$2
    
    output_path=`dirname $input_file`
    gencode_gff3_file=`basename $gencode_gff3`

    test -d ${output_path}/Anno || mkdir -p ${output_path}/Anno
    local output_prefix="${output_path}/Anno/${gencode_gff3_file%.*}"

    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCreating annotation files of gene elements." > /dev/tty
    test -f ${output_prefix}_gene.bed || awk -v OFS="\t" '{split($9,a,";");split(a[2],d,"=");split(a[3],b,"=");split(a[4],c,"=");if($3=="gene") print $1,$4-1,$5,d[2]":"c[2],b[2],$7}' ${gencode_gff3} > ${output_prefix}_gene.bed &
    test -f ${output_prefix}_exon.bed || awk -v OFS="\t" '{split($9,a,";");split(a[10],d,"=");split(a[7],b,"=");split(a[3],c,"=");if($3=="exon") print $1,$4-1,$5,d[2]":"c[2],b[2],$7}' ${gencode_gff3} | grep -v retained_intron | bedtools sort -i | uniq > ${output_prefix}_exon.bed &
    test -f ${output_prefix}_CDS.bed  || awk -v OFS="\t" '{split($9,a,";");split(a[10],d,"=");split(a[7],b,"=");split(a[3],c,"=");if($3=="CDS") print $1,$4-1,$5,d[2]":"c[2],b[2],$7}' ${gencode_gff3} | grep -v retained_intron | bedtools sort -i | uniq > ${output_prefix}_CDS.bed &
    test -f ${output_prefix}_5UTR.bed || awk -v OFS="\t" '{split($9,a,";");split(a[10],d,"=");split(a[7],b,"=");split(a[3],c,"=");if($3=="five_prime_UTR") print $1,$4-1,$5,d[2]":"c[2],b[2],$7}' ${gencode_gff3} | grep -v retained_intron | bedtools sort -i | uniq > ${output_prefix}_5UTR.bed &
    test -f ${output_prefix}_3UTR.bed || awk -v OFS="\t" '{split($9,a,";");split(a[10],d,"=");split(a[7],b,"=");split(a[3],c,"=");if($3=="three_prime_UTR") print $1,$4-1,$5,d[2]":"c[2],b[2],$7}' ${gencode_gff3} | grep -v retained_intron | bedtools sort -i | uniq > ${output_prefix}_3UTR.bed &
    test -f ${output_prefix}_stopcodon.bed || awk -v OFS="\t" '{split($9,a,";");split(a[10],d,"=");split(a[3],e,"=");split(a[5],b,"=");split(a[7],c,"=");if($3=="stop_codon") print $1,$4-1,$5,d[2]":"e[2],c[2]":"b[2],$7 }'  ${gencode_gff3} | bedtools sort -i - | uniq > ${output_prefix}_stopcodon.bed &
    test -f ${output_prefix}_startcodon.bed || awk -v OFS="\t" '{split($9,a,";");split(a[10],d,"=");split(a[3],e,"=");split(a[5],b,"=");split(a[7],c,"=");if($3=="start_codon") print $1,$4-1,$5,d[2]":"e[2],c[2]":"b[2],$7 }'  ${gencode_gff3} | bedtools sort -i - | uniq > ${output_prefix}_startcodon.bed &
    wait
    
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tAnnotating edit data with gene elements." > /dev/tty
    bedtools subtract -a ${output_prefix}_gene.bed -b ${output_prefix}_exon.bed -s > ${output_prefix}_intron.bed
    cat ${output_prefix}_5UTR.bed ${output_prefix}_CDS.bed ${output_prefix}_3UTR.bed > ${output_prefix}_mRNA.bed
    awk 'FILENAME==ARGV[1]{array_tmp[$4]=1}FILENAME==ARGV[2]{if (! array_tmp[$4] ){print}}' ${output_prefix}_mRNA.bed ${output_prefix}_exon.bed > ${output_prefix}_noncoding.bed
    rm ${output_prefix}_mRNA.bed
    
    bedtools intersect -a $input_file -b ${output_prefix}_gene.bed -s -wa -wb | cut -f 1-6,10,11 | bedtools sort -i - | uniq > ${input_file%.*}_gene.bed
    bedtools subtract -a $input_file -b ${output_prefix}_gene.bed  | bedtools sort -i - | uniq > ${input_file%.*}_IGR.bed
    bedtools intersect -a ${input_file%.*}_gene.bed -b ${output_prefix}_exon.bed -s -wa -wb | cut -f 1-8,12,13  | bedtools sort -i - | uniq > ${input_file%.*}_exon.bed
    bedtools intersect -a ${input_file%.*}_gene.bed -b ${output_prefix}_intron.bed -s -wa -wb | cut -f 1-8  | bedtools sort -i - | uniq > ${input_file%.*}_intron.bed
    bedtools intersect -a ${input_file%.*}_exon.bed -b ${output_prefix}_CDS.bed -s -wa -wb | awk '{if(($10==$15)&&($9==$14)) print $0}' |  cut -f 1-10 | uniq > ${input_file%.*}_CDS.bed
    bedtools intersect -a ${input_file%.*}_exon.bed -b ${output_prefix}_3UTR.bed -s -wa -wb | awk '{if(($10==$15)&&($9==$14)) print $0}' |  cut -f 1-10 | uniq > ${input_file%.*}_3UTR.bed
    bedtools intersect -a ${input_file%.*}_exon.bed -b ${output_prefix}_5UTR.bed -s -wa -wb | awk '{if(($10==$15)&&($9==$14)) print $0}' |  cut -f 1-10 | uniq > ${input_file%.*}_5UTR.bed
    bedtools intersect -a ${input_file%.*}_exon.bed -b ${output_prefix}_noncoding.bed -s -wa -wb | awk '{if(($10==$15)&&($9==$14)) print $0}' |  cut -f 1-10 | uniq > ${input_file%.*}_noncoding.bed
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tAnnotating done." > /dev/tty
}

Pre_for_SAILOR(){
    local Mapit_result=$1
    local Sample_name=$2
    local Replicate=$3
    local REFERENCE=$4
    local SIZE_FILE=$5
    local NUM_CORES=$6
    local FLARE_PATH=$7
    local Reditools_PATH=$8

    # 0_MarkDuplicates
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tMarking duplicates." > /dev/tty
    combine_bam=${Mapit_result}/3-Combine_bam/${Sample_name}/${Sample_name}_${Replicate}_combine.bam
    dedupped_bam=${Mapit_result}/3-Combine_bam/${Sample_name}/${Sample_name}_${Replicate}_combine_dedupped.bam
    ( test -f ${dedupped_bam} || test -f ${combine_bam} )|| (echo -e "\033[31m[ERROR]\033[0m\t\033[32m`date '+%a %b %d %H:%M:S %Y'`\033[0m\t"${combine_bam}" not exits. Run mapping first" > /dev/tty && exit 1)
    test -f ${dedupped_bam} || picard MarkDuplicates -I ${combine_bam} -O ${dedupped_bam} \
                                                     --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
                                                     -M ${Mapit_result}/3-Combine_bam/${Sample_name}/${Sample_name}_${Replicate}_MarkDuplicates_output.metrics \
                                                     --REMOVE_DUPLICATES true 

    # 1_split_strands3.1-Split_strand
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tSpliting into positive/negative strand bam." > /dev/tty
    BAM_PREFIX=${Mapit_result}/3.1-Split_strand/${Sample_name}
    VAR_PREFIX=${Mapit_result}/4.1-Redi_sailor/${Sample_name}
    test -d ${BAM_PREFIX} || mkdir -p ${BAM_PREFIX}
    test -d ${VAR_PREFIX} || mkdir -p ${VAR_PREFIX}
    #samtools view -@ $NUM_CORES -F 1024 -hb ${dedupped_bam} \
    #    -o ${BAM_PREFIX}/${Sample_name}_${Replicate}.bam 
    ln -s ${dedupped_bam} ${BAM_PREFIX}/${Sample_name}_${Replicate}.bam
    samtools index -@ $NUM_CORES ${BAM_PREFIX}/${Sample_name}_${Replicate}.bam 
    python ${FLARE_PATH}/workflow_sailor/scripts/split_strands.py --reverse-strand \
        -i ${BAM_PREFIX}/${Sample_name}_${Replicate}.bam  \
        -f ${BAM_PREFIX}/${Sample_name}_${Replicate}.fwd.bam \
        -r ${BAM_PREFIX}/${Sample_name}_${Replicate}.rev.bam
    # 2_index_reads
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tExtracting coverage of positive/negative strand bams following Reditools2.0 workflow." > /dev/tty
    samtools index -@ $NUM_CORES ${BAM_PREFIX}/${Sample_name}_${Replicate}.fwd.bam
    samtools index -@ $NUM_CORES ${BAM_PREFIX}/${Sample_name}_${Replicate}.rev.bam

    test -d ${VAR_PREFIX}/${Sample_name}_${Replicate}.fwd.coverage/ || mkdir -p ${VAR_PREFIX}/${Sample_name}_${Replicate}.fwd.coverage/
    ${Reditools_PATH}/extract_coverage.sh ${BAM_PREFIX}/${Sample_name}_${Replicate}.fwd.bam ${VAR_PREFIX}/${Sample_name}_${Replicate}.fwd.coverage/ $SIZE_FILE
    test -d ${VAR_PREFIX}/${Sample_name}_${Replicate}.rev.coverage/ || mkdir -p ${VAR_PREFIX}/${Sample_name}_${Replicate}.rev.coverage/
    ${Reditools_PATH}/extract_coverage.sh ${BAM_PREFIX}/${Sample_name}_${Replicate}.rev.bam ${VAR_PREFIX}/${Sample_name}_${Replicate}.rev.coverage/ $SIZE_FILE

    # 3_bw (8 in raw sailor)
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tConverting positive/negative strand coverage into bigWig." > /dev/tty
    test -d ${Mapit_result}/3.1-Split_strand/chrom.sizes || cut -f 1,2 $SIZE_FILE > ${Mapit_result}/3.1-Split_strand/chrom.sizes
    bamCoverage -p $NUM_CORES -b ${BAM_PREFIX}/${Sample_name}_${Replicate}.fwd.bam --binSize 1 --normalizeUsing None -o ${BAM_PREFIX}/${Sample_name}_${Replicate}.fwd.sorted.bw 
    bamCoverage -p $NUM_CORES -b ${BAM_PREFIX}/${Sample_name}_${Replicate}.rev.bam --binSize 1 --normalizeUsing None -o ${BAM_PREFIX}/${Sample_name}_${Replicate}.rev.sorted.bw 
}

SAILOR_for_MAPIT(){
    local Mapit_result=$1
    local Sample_name=$2
    local Replicate=$3
    local REFERENCE=$4
    local SIZE_FILE=$5
    local NUM_CORES=$6
    local SNP_dbSNP_divided_by_chromosome=$7
    local SNP_1000Genome_divided_by_chromosome=$8
    local SNP_EVS_divided_by_chromosome=$9
    local FLARE_PATH=${10}
    local Reditools_PATH=${11}

    BAM_PREFIX=${Mapit_result}/3.1-Split_strand/${Sample_name}
    VAR_PREFIX=${Mapit_result}/4.1-Redi_sailor/${Sample_name}
    ## 4.1 reditools 
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tRNA editing detecting." > /dev/tty
    for STRAND in fwd rev
    do
    SOURCE_BAM_FILE=${BAM_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.bam
    OUTPUT_FILE=${VAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}_parallel_table.txt.gz
    TEMP_DIR=${VAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.tmp/
    COVERAGE_DIR=${VAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.coverage/
    COVERAGE_FILE=${COVERAGE_DIR}/${Sample_name}_${Replicate}.${STRAND}.cov
    test -d ${COVERAGE_DIR} || mkdir -p ${COVERAGE_DIR}
    test -d ${TEMP_DIR} || mkdir -p ${TEMP_DIR}
    
    #${Reditools_PATH}/extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
    mpirun -np $NUM_CORES ${Reditools_PATH}/src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR -q 10 -bq 20 
    ${Reditools_PATH}/merge.sh $TEMP_DIR $OUTPUT_FILE $NUM_CORES
    done
    
    ## 4.2 select MAPIT editing
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tMAPIT editing filtering." > /dev/tty
    zcat ${VAR_PREFIX}/${Sample_name}_${Replicate}.fwd_parallel_table.txt.gz | awk '/AG|CT/{if($5>=10) print }' > ${VAR_PREFIX}/${Sample_name}_${Replicate}.fwd_parallel_table_dp10.txt
    zcat ${VAR_PREFIX}/${Sample_name}_${Replicate}.rev_parallel_table.txt.gz | awk '/GA|TC/{if($5>=10) print }' > ${VAR_PREFIX}/${Sample_name}_${Replicate}.rev_parallel_table_dp10.txt
    
    # 5 drop SNPs & format transforming and combine
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tRemoving SNPs in MAPIT editing." > /dev/tty
    FVAR_PREFIX=${Mapit_result}/5.1-Var_filter/${Sample_name}
    EDIT_PREFIX=${Mapit_result}/6.1-Edit_bedgraphs/${Sample_name}
    EDIT_PREFIX2=${Mapit_result}/6.2-Edit_bigwig/${Sample_name}
    test -d ${EDIT_PREFIX2} || mkdir -p ${EDIT_PREFIX2}
    test -d ${FVAR_PREFIX} || mkdir -p ${FVAR_PREFIX}
    test -d ${EDIT_PREFIX} || mkdir -p ${EDIT_PREFIX}
    for STRAND in fwd rev
    do
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\t"${Sample_name}_${Replicate}.${STRAND}"\tStart to filter out known SNPs."
    Filtering_SNV ${VAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}_parallel_table_dp10.txt ${FVAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.dbSNP.txt ${SNP_dbSNP_divided_by_chromosome} ${NUM_CORES}
    Filtering_SNV ${FVAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.dbSNP.txt ${FVAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.dbSNP_1000genomes.txt ${SNP_1000Genome_divided_by_chromosome} ${NUM_CORES}
    Filtering_SNV ${FVAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.dbSNP_1000genomes.txt ${FVAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.dbSNP_1000genomes_EVS.txt ${SNP_EVS_divided_by_chromosome} ${NUM_CORES}
    python ${MYDIR}/rank_edits_based_on_reditools.py -i ${FVAR_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.dbSNP_1000genomes_EVS.txt -o ${EDIT_PREFIX}/${Sample_name}_${Replicate}.${STRAND}.snpfiltered.ranked.bed
    done
    
    grep "C>T" ${EDIT_PREFIX}/${Sample_name}_${Replicate}.fwd.snpfiltered.ranked.bed > ${EDIT_PREFIX}/${Sample_name}_${Replicate}.fwd.CT.snpfiltered.ranked.bed
    grep "G>A" ${EDIT_PREFIX}/${Sample_name}_${Replicate}.rev.snpfiltered.ranked.bed > ${EDIT_PREFIX}/${Sample_name}_${Replicate}.rev.CT.snpfiltered.ranked.bed
    grep "A>G" ${EDIT_PREFIX}/${Sample_name}_${Replicate}.fwd.snpfiltered.ranked.bed > ${EDIT_PREFIX}/${Sample_name}_${Replicate}.fwd.AG.snpfiltered.ranked.bed
    grep "T>C" ${EDIT_PREFIX}/${Sample_name}_${Replicate}.rev.snpfiltered.ranked.bed > ${EDIT_PREFIX}/${Sample_name}_${Replicate}.rev.AG.snpfiltered.ranked.bed
    
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tCombining and reformating." > /dev/tty
    for EDIT in CT AG
    do
    python ${FLARE_PATH}/workflow_sailor/scripts/combine_and_reformat.py --fwd ${EDIT_PREFIX}/${Sample_name}_${Replicate}.fwd.${EDIT}.snpfiltered.ranked.bed \
                                                                         --rev ${EDIT_PREFIX}/${Sample_name}_${Replicate}.rev.${EDIT}.snpfiltered.ranked.bed \
                                                                         --output ${EDIT_PREFIX}/${Sample_name}_${Replicate}.combined.${EDIT}.snpfiltered.ranked.bed
    python ${FLARE_PATH}/workflow_sailor/scripts/edit_fraction_bedgraph.py ${EDIT_PREFIX}/${Sample_name}_${Replicate}.combined.${EDIT}.snpfiltered.ranked.bed \
                                                                           ${EDIT_PREFIX}/${Sample_name}_${Replicate}.${EDIT}.edit_fraction.bedgraph
    awk -v OFS="\t" '{split($5,a,","); if(($4>0.9)&&(a[1]/a[2]<0.7)) print}' ${EDIT_PREFIX}/${Sample_name}_${Replicate}.combined.${EDIT}.snpfiltered.ranked.bed > ${EDIT_PREFIX}/${Sample_name}_${Replicate}.combined.${EDIT}.snpfiltered.ranked.score0.9.bed
    awk -v OFS="\t" '{split($5,a,","); if(($4>0.5)&&(a[1]/a[2]<0.7)) print}' ${EDIT_PREFIX}/${Sample_name}_${Replicate}.combined.${EDIT}.snpfiltered.ranked.bed > ${EDIT_PREFIX}/${Sample_name}_${Replicate}.combined.${EDIT}.snpfiltered.ranked.score0.5.bed
    bedtools sort -i ${EDIT_PREFIX}/${Sample_name}_${Replicate}.${EDIT}.edit_fraction.bedgraph > ${EDIT_PREFIX2}/${Sample_name}_${Replicate}.${EDIT}.edit_fraction.sorted.bedgraph
    bedGraphToBigWig ${EDIT_PREFIX2}/${Sample_name}_${Replicate}.${EDIT}.edit_fraction.sorted.bedgraph ${Mapit_result}/3.1-Split_strand/chrom.sizes ${EDIT_PREFIX2}/${Sample_name}_${Replicate}.${EDIT}.edit_fraction.sorted.bw
    done
    echo -e "\033[34m[INFO]\033[0m\t\033[32m`date '+%a %b %d %H:%M:%S %Y'`\033[0m\tDone." > /dev/tty
}

FLARE_for_MAPIT(){
    local Mapit_result=$1
    local Sample_name=$2
    local Replicate=$3
    local REFERENCE=$4
    local FLARE_REGIONS=$5
    local EDIT_TYPE=$6
    local NUM_CORES=$7
    local FLARE_PATH=$8
    
    test -d ${Mapit_result}/7-FLARE_for_MAPIT/${Sample_name}/ || mkdir -p ${Mapit_result}/7-FLARE_for_MAPIT/${Sample_name}/
    json_file=${Mapit_result}/7-FLARE_for_MAPIT/${Sample_name}/${Sample_name}_${Replicate}_${EDIT_TYPE}.json
    if [ $EDIT_TYPE == "CT" ];then
        score=0.5
    elif [ $EDIT_TYPE == "AG" ];then
        score=0.9
    fi
    
    test -f ${json_file} || awk -v sample=${Sample_name} -v rep=${Replicate} -v path=${Mapit_result} -v genome=${REFERENCE} -v regions=${FLARE_REGIONS} -v edittype=${EDIT_TYPE} -v X=${score} \
        '{gsub(/sample/,sample,$0);gsub(/replicate/,rep,$0);gsub(/path/,path,$0);gsub(/edittype/,edittype,$0);gsub(/0.X/,X,$0);gsub(/genome.fa/,genome,$0);gsub(/FLARE_regions/,regions,$0);
         print}' ${MYDIR}/../conf/FLARE_setting.json > ${json_file}
    snakemake --snakefile ${FLARE_PATH}/workflow_FLARE/Snakefile --configfile ${json_file} -j ${NUM_CORES}
}

HC_cluster(){
    local Mapit_result=$1
    local Sample_name=$2
    local slop_length=$3
    FLARE_dir=${Mapit_result}"/7-FLARE_for_MAPIT/"
    hc_dir=${Mapit_result}"/7-FLARE_for_MAPIT/confident_clusters"
    test -d ${hc_dir} || mkdir -p ${hc_dir}

    for EDIT in CT AG
    do
    bedtools intersect -a ${FLARE_dir}/${Sample_name}/${EDIT}/FLARE/${Sample_name}_1_${EDIT}/${Sample_name}_1_${EDIT}_merged_sorted_peaks.fdr_0.1.d_15.merged.bed \
                    -b ${FLARE_dir}/${Sample_name}/${EDIT}/FLARE/${Sample_name}_2_${EDIT}/${Sample_name}_2_${EDIT}_merged_sorted_peaks.fdr_0.1.d_15.merged.bed \
                    -s -wa -u > ${hc_dir}/${Sample_name}_1_${EDIT}.bed
    bedtools intersect -b ${FLARE_dir}/${Sample_name}/${EDIT}/FLARE/${Sample_name}_1_${EDIT}/${Sample_name}_1_${EDIT}_merged_sorted_peaks.fdr_0.1.d_15.merged.bed \
                    -a ${FLARE_dir}/${Sample_name}/${EDIT}/FLARE/${Sample_name}_2_${EDIT}/${Sample_name}_2_${EDIT}_merged_sorted_peaks.fdr_0.1.d_15.merged.bed \
                    -s -wa -u > ${hc_dir}/${Sample_name}_2_${EDIT}.bed
    cat ${hc_dir}/${Sample_name}_1_${EDIT}.bed ${hc_dir}/${Sample_name}_2_${EDIT}.bed \
    | bedtools sort -i -  | bedtools merge -i - -s -c 4,5,6,7 -o collapse,collapse,distinct,collapse > ${hc_dir}/${Sample_name}.${EDIT}.bed
    done

    bedtools slop -i ${hc_dir}/${Sample_name}.AG.bed -g ${Mapit_result}/3.1-Split_strand/chrom.sizes -b 100 | bedtools intersect -a ${hc_dir}/${Sample_name}.CT.bed -b - -s -wa -u > ${hc_dir}/${Sample_name}.highconfidence.bed
    #cat ${hc_dir}/${Sample_name}.AG.bed ${hc_dir}/${Sample_name}.CT.bed | bedtools sort -i - | bedtools merge -i - -s -c 4,5,6,7 -o collapse,collapse,distinct,collapse > ${hc_dir}/${Sample_name}.merge.bed
    awk -v OFS="\t" '{center=int(($2+$3)/2);print $1,center-1,center,$4,$5,$6}' ${hc_dir}/${Sample_name}.highconfidence.bed | bedtools slop -i - -g ${Mapit_result}/3.1-Split_strand/chrom.sizes -b ${slop_length} > ${hc_dir}/${Sample_name}.highconfidence.slop${slop_length}.bed
}

"$@"