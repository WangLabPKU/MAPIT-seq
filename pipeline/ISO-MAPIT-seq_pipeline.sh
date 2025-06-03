# 0. Set up
genomeFasta=/home/gangx/data/Reference/GRCh38/GRCh38.p13.genome.pri.fa
genomeIndex=/home/gangx/data/Reference/GRCh38/GRCh38.p13.genome.pri.mmi
geneAnnoGtf=/home/gangx/data/Reference/GRCh38/gencode.v40.pri.annotation.sorted.gtf
geneAnnoDB=/home/gangx/data/Reference/GRCh38/gencode.v40.pri.annotation.sorted.db
juncBed=/home/gangx/data/Reference/GRCh38/gencode.v40.pri.annotation.junc.bed

pbmm2 index $genomeFasta ${genomeIndex} --preset CCS -j 20
pigeon prepare $geneAnnoGtf $genomeFasta
paftools.js gff2bed $geneAnnoGtf > $juncBed

# 1. Refine
isoseq refine -j 20 --require-polya raw/DG/DG.hifi.bam IsoSeq_v2_primers_12.fasta raw/DG/DG.hifi.flnc.bam
isoseq refine -j 20 --require-polya raw/DI/DI.hifi.bam IsoSeq_v2_primers_12.fasta raw/DI/DI.hifi.flnc.bam
isoseq refine -j 20 --require-polya raw/FG/FG.hifi.bam     IsoSeq_v2_primers_12.fasta raw/FG/FG.hifi.flnc.bam
isoseq refine -j 20 --require-polya raw/FI/FI.hifi.bam     IsoSeq_v2_primers_12.fasta raw/FI/FI.hifi.flnc.bam

# 2. Covert format
bamToFastq -i raw/DG/DG.hifi.flnc.bam -fq /dev/stdout | pigz - > fastq/DG.hifi.flnc.fq.gz
bamToFastq -i raw/DI/DI.hifi.flnc.bam -fq /dev/stdout | pigz - > fastq/DI.hifi.flnc.fq.gz
bamToFastq -i raw/FG/FG.hifi.flnc.bam     -fq /dev/stdout | pigz - > fastq/FG.hifi.flnc.fq.gz
bamToFastq -i raw/FI/FI.hifi.flnc.bam     -fq /dev/stdout | pigz - > fastq/FI.hifi.flnc.fq.gz

# 3. Mapping
minimap2 $genomeFasta fastq/DG.hifi.flnc.fq.gz -ax splice:hq -u f --secondary=no -t 40 -R "@RG\\tID:DG\\tSM:DG" --MD --junc-bed ${juncBed} | samtools view -@ 40 -hb - -o minimap2_align/DG.bam
minimap2 $genomeFasta fastq/DI.hifi.flnc.fq.gz -ax splice:hq -u f --secondary=no -t 40 -R "@RG\\tID:DI\\tSM:DI" --MD --junc-bed ${juncBed} | samtools view -@ 40 -hb - -o minimap2_align/DI.bam
minimap2 $genomeFasta fastq/FG.hifi.flnc.fq.gz -ax splice:hq -u f --secondary=no -t 40 -R "@RG\\tID:FG\\tSM:FG" --MD --junc-bed ${juncBed} | samtools view -@ 40 -hb - -o minimap2_align/FG.bam
minimap2 $genomeFasta fastq/FI.hifi.flnc.fq.gz -ax splice:hq -u f --secondary=no -t 40 -R "@RG\\tID:FI\\tSM:FI" --MD --junc-bed ${juncBed} | samtools view -@ 40 -hb - -o minimap2_align/FI.bam

# 4. Sorting & drop chrM
samtools view -@ 20 -h minimap2_align/DG.bam |  awk '/@/ || ($3 != "chrM") {print $0}' | samtools sort -@ 20 - -o minimap2_align/DG.sorted.bam
samtools view -@ 20 -h minimap2_align/DG.bam | awk '/@HD/ || /@PG/{print $0}/@SQ/&&/chrM/{print $0}$3 == "chrM" {print $0}' | samtools sort -@ 20 - -o minimap2_align/DG.chrM.bam

samtools view -@ 20 -h minimap2_align/DI.bam |  awk '/@/ || ($3 != "chrM") {print $0}' | samtools sort -@ 20 - -o minimap2_align/DI.sorted.bam
samtools view -@ 20 -h minimap2_align/DI.bam | awk '/@HD/ || /@PG/{print $0}/@SQ/&&/chrM/{print $0}$3 == "chrM" {print $0}' | samtools sort -@ 20 - -o minimap2_align/DI.chrM.bam

samtools view -@ 20 -h minimap2_align/FG.bam |  awk '/@/ || ($3 != "chrM") {print $0}' | samtools sort -@ 20 - -o minimap2_align/FG.sorted.bam
samtools view -@ 20 -h minimap2_align/FG.bam | awk '/@HD/ || /@PG/{print $0}/@SQ/&&/chrM/{print $0}$3 == "chrM" {print $0}' | samtools sort -@ 20 - -o minimap2_align/FG.chrM.bam

samtools view -@ 20 -h minimap2_align/FI.bam |  awk '/@/ || ($3 != "chrM") {print $0}' | samtools sort -@ 20 - -o minimap2_align/FI.sorted.bam
samtools view -@ 20 -h minimap2_align/FI.bam | awk '/@HD/ || /@PG/{print $0}/@SQ/&&/chrM/{print $0}$3 == "chrM" {print $0}' | samtools sort -@ 20 - -o minimap2_align/FI.chrM.bam

# 5. Perform isoquant
isoquant.py -t 20 --reference $genomeFasta --genedb $geneAnnoDB --junc_bed_file $juncBed --complete_genedb \
            --bam minimap2_align/DG.sorted.bam -p DG \
            --data_type pacbio_ccs --transcript_quantification unique_only --gene_quantification unique_splicing_consistent --splice_correction_strategy default_pacbio \
            --normalization_method usable_reads --labels DG -o isoquant &
sleep 10
isoquant.py -t 20 --reference $genomeFasta --genedb $geneAnnoDB --junc_bed_file $juncBed --complete_genedb \
            --bam minimap2_align/DI.sorted.bam -p DI \
            --data_type pacbio_ccs --transcript_quantification unique_only --gene_quantification unique_splicing_consistent --splice_correction_strategy default_pacbio \
            --normalization_method usable_reads --labels DI -o isoquant &
sleep 10
isoquant.py -t 20 --reference $genomeFasta --genedb $geneAnnoDB --junc_bed_file $juncBed --complete_genedb \
            --bam minimap2_align/FG.sorted.bam -p FG_1 \
            --data_type pacbio_ccs --transcript_quantification unique_only --gene_quantification unique_splicing_consistent --splice_correction_strategy default_pacbio \
            --normalization_method usable_reads --labels FG_1 -o isoquant &
sleep 10
isoquant.py -t 20 --reference $genomeFasta --genedb $geneAnnoDB --junc_bed_file $juncBed --complete_genedb \
            --bam minimap2_align/FI.sorted.bam -p FI_1 \
            --data_type pacbio_ccs --transcript_quantification unique_only --gene_quantification unique_splicing_consistent --splice_correction_strategy default_pacbio \
            --normalization_method usable_reads --labels FI_1 -o isoquant &
sleep 10
wait

# 6. Perform isoform-mapit-editing analysis 
~/apps/Mapit-seq/isomapit.py isoedit -v GRCh38 -i minimap2_align/DG.sorted.bam -r isoquant/DG/DG.read_assignments.tsv.gz -o ./isomapit -d DG -t 23
~/apps/Mapit-seq/isomapit.py isoedit -v GRCh38 -i minimap2_align/DI.sorted.bam -r isoquant/DI/DI.read_assignments.tsv.gz -o ./isomapit -d DI -t 23
~/apps/Mapit-seq/isomapit.py isoedit -v GRCh38 -i minimap2_align/FG.sorted.bam -r isoquant/FG/FG.read_assignments.tsv.gz -o ./isomapit -d FG -t 23
~/apps/Mapit-seq/isomapit.py isoedit -v GRCh38 -i minimap2_align/FI.sorted.bam -r isoquant/FI/FI.read_assignments.tsv.gz -o ./isomapit -d FI -t 23
