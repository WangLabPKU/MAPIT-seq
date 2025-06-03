# 1.Mapping
STAR  --genomeDir ref/star --readFilesCommand zcat \
 --readFilesIn /home/gangx/data/CQX/240612/fastq/G/G_R2.fq.gz,/home/gangx/data/CQX/240612/fastq2/G/G_R2.fq.gz /home/gangx/data/CQX/240612/fastq/G/G_R1.fq.gz,/home/gangx/data/CQX/240612/fastq2/G/G_R1.fq.gz \
 --soloCBwhitelist ./3M-february-2018.txt \
 --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
 --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
 --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 \
 --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
 --soloMultiMappers EM --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.04 --outFilterType BySJout \
 --outSAMattributes CB UB --outSAMattrRGline ID:G3BP1 SM:G3BP1 PL:Illumina PU:10x  \
 --outFileNamePrefix align/STARsolo.G. --outSAMtype BAM SortedByCoordinate --runThreadN 60

STAR  --genomeDir ref/star --readFilesCommand zcat \
 --readFilesIn /home/gangx/data/CQX/240612/fastq/I/I_R2.fq.gz,/home/gangx/data/CQX/240612/fastq2/I/I_R2.fq.gz /home/gangx/data/CQX/240612/fastq/I/I_R1.fq.gz,/home/gangx/data/CQX/240612/fastq2/I/I_R1.fq.gz \
 --soloCBwhitelist ./3M-february-2018.txt \
 --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
 --soloUMIdedup 1MM_CR --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
 --soloCellFilter EmptyDrops_CR --outFilterScoreMin 30 \
 --soloFeatures Gene GeneFull --soloOutFileNames output/ features.tsv barcodes.tsv matrix.mtx \
 --soloMultiMappers EM --outReadsUnmapped Fastx --outFilterMismatchNoverLmax 0.04 --outFilterType BySJout \
 --outSAMattributes CB UB --outSAMattrRGline ID:IgG SM:IgG PL:Illumina PU:10x \
 --outFileNamePrefix align/STARsolo.I. --outSAMtype BAM SortedByCoordinate --runThreadN 80

samtools index -@ 20 align/STARsolo.G.Aligned.sortedByCoord.out.bam
samtools index -@ 20 align/STARsolo.I.Aligned.sortedByCoord.out.bam

# 2. Cell quality control
# See 10x-scMAPIT-seq_scanpy.ipynb

# 3. BAM filtering cells
# See 10x-scMAPIT-seq_filterbam&scRNAediting.ipynb
samtools index -@ 20 fileredbam/G.bam 
samtools index -@ 20 fileredbam/I.bam

# 4. Drop PCR duplicates
umi_tools dedup -I fileredbam/G.bam -S dedup/G.bam --extract-umi-method=tag --umi-tag=UB --cell-tag=CB --per-cell 
umi_tools dedup -I fileredbam/I.bam -S dedup/I.bam --extract-umi-method=tag --umi-tag=UB --cell-tag=CB --per-cell 
samtools index -@ 20 dedup/G.bam 
samtools index -@ 20 dedup/I.bam 

# 5. BAM spliting to single-cell BAM
bamtools split -stub splitbam/G/ -in dedup/G.bam -tag CB
bamtools split -stub splitbam/I/ -in dedup/I.bam -tag CB

# 6. Perform single cell RNA editing 
# See 10x-scMAPIT-seq_filterbam&scRNAediting.ipynb

# 7. Cell cycle analysis
# See 10x-scMAPIT-seq_Seurat.ipynb

# 8. Cell cycle specific editing analysis
samtools merge -@ 20 -o ccbam_seurat/G_Early-G1.bam  `awk '$2=="Early-G1"  && $3=="Anti-G3BP1"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/G_Late-G1.bam   `awk '$2=="Late-G1"   && $3=="Anti-G3BP1"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/G_Early-S.bam   `awk '$2=="Early-S"   && $3=="Anti-G3BP1"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/G_Late-S.bam    `awk '$2=="Late-S"    && $3=="Anti-G3BP1"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/G_Early-G2M.bam `awk '$2=="Early-G2M" && $3=="Anti-G3BP1"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/G_Late-G2M.bam  `awk '$2=="Late-G2M"  && $3=="Anti-G3BP1"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/I_Early-G1.bam  `awk '$2=="Early-G1"  && $3=="IgG"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/I_Late-G1.bam   `awk '$2=="Late-G1"   && $3=="IgG"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/I_Early-S.bam   `awk '$2=="Early-S"   && $3=="IgG"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/I_Late-S.bam    `awk '$2=="Late-S"    && $3=="IgG"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/I_Early-G2M.bam `awk '$2=="Early-G2M" && $3=="IgG"{print $4}' barcode_ccbam_seurat.tsv`
samtools merge -@ 20 -o ccbam_seurat/I_Late-G2M.bam  `awk '$2=="Late-G2M"  && $3=="IgG"{print $4}' barcode_ccbam_seurat.tsv`

for phase in Early-G1 Late-G1 Early-S Late-S Early-G2M Late-G2M 
do
SOURCE_BAM_FILE="/home/gangx/data/CQX/240612/ccbam_seurat/G_"${phase}".bam"
REFERENCE="/home/gangx/data/Reference/GRCh38/GRCh38.p13.genome.pri.fa"
SIZE_FILE="/home/gangx/data/Reference/GRCh38/GRCh38.p13.genome.pri.fa.fai"

NUM_CORES=40
OUTPUT_FILE="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/G_"${phase}"_parallel_table.txt.gz"
TEMP_DIR="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/G_"${phase}"_temp/"
COVERAGE_FILE="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/coverage_G_"${phase}"/G_"${phase}".cov"
COVERAGE_DIR="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/coverage_G_"${phase}"/"

samtools index -@ 20  $SOURCE_BAM_FILE
/home/gangx/apps/reditools2.0/extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
mpirun -np $NUM_CORES /home/gangx/apps/reditools2.0/src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR
/home/gangx/apps/reditools2.0/merge.sh $TEMP_DIR $OUTPUT_FILE $NUM_CORES

zcat $OUTPUT_FILE | awk -F "\t" -v OFS="\t" '($5 >= 10) && (length($8) > 1){print $1,$2-1,$2,$5,$7,$9,$3,$8}' > cc_redi/G_${phase}_edit.bed
awk -F "\t" -v OFS="\t" '/AG|CT/{print $1,$2,$3,$4,$5,"+",$6,$7,$8}/TC|GA/{print $1,$2,$3,$4,$5,"-",$6,$7,$8}' cc_redi/G_${phase}_edit.bed > cc_redi/G_${phase}_edit2.bed
done

for phase in Early-G1 Late-G1 Early-S Late-S Early-G2M Late-G2M 
do
SOURCE_BAM_FILE="/home/gangx/data/CQX/240612/ccbam_seurat/I_"${phase}".bam"
REFERENCE="/home/gangx/data/Reference/GRCh38/GRCh38.p13.genome.pri.fa"
SIZE_FILE="/home/gangx/data/Reference/GRCh38/GRCh38.p13.genome.pri.fa.fai"

NUM_CORES=40
OUTPUT_FILE="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/I_"${phase}"_parallel_table.txt.gz"
TEMP_DIR="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/I_"${phase}"_temp/"
COVERAGE_FILE="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/coverage_I_"${phase}"/I_"${phase}".cov"
COVERAGE_DIR="/home/gangx/data/CQX/240612/ccbam_seurat/cc_redi/coverage_I_"${phase}"/"

samtools index -@ 20  $SOURCE_BAM_FILE
/home/gangx/apps/reditools2.0/extract_coverage.sh $SOURCE_BAM_FILE $COVERAGE_DIR $SIZE_FILE
mpirun -np $NUM_CORES /home/gangx/apps/reditools2.0/src/cineca/parallel_reditools.py -f $SOURCE_BAM_FILE -o $OUTPUT_FILE -r $REFERENCE -t $TEMP_DIR -Z $SIZE_FILE -G $COVERAGE_FILE -D $COVERAGE_DIR
/home/gangx/apps/reditools2.0/merge.sh $TEMP_DIR $OUTPUT_FILE $NUM_CORES

zcat $OUTPUT_FILE | awk -F "\t" -v OFS="\t" '($5 >= 10) && (length($8) > 1){print $1,$2-1,$2,$5,$7,$9,$3,$8}' > cc_redi/I_${phase}_edit.bed
awk -F "\t" -v OFS="\t" '/AG|CT/{print $1,$2,$3,$4,$5,"+",$6,$7,$8}/TC|GA/{print $1,$2,$3,$4,$5,"-",$6,$7,$8}' cc_redi/I_${phase}_edit.bed > cc_redi/I_${phase}_edit2.bed
done

# 9. MAPIT-seq data call targets 
# See MAPIT-seq `calltargets` part.