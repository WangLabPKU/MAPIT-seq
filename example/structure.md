## `Mapit_result/`

### `mapping` and `finetuning` steps

#### `0-Remove_rRNA/`

* **Purpose:** Removal of abundant RNA reads before alignment.
* `G3BP1_1_R1.fastq.gz`, `G3BP1_1_R2.fastq.gz`: Paired-end reads after rRNA removal.
* `G3BP1_1_singleton.fq`: Unpaired reads after filtering.

#### `1-HISAT_map/`

* **Purpose:** Primary transcriptome alignment using HISAT2.
* `*.header`: BAM header file.
* `*_mapped.summary`: Mapping statistics summary.
* `*_mapped-unmapped_sorted.bam`: BAM file containing both mapped and unmapped reads.
* `*_unmapped_*.fastq.gz`, `singleton.fq`: Reads that failed to map with HISAT2.

#### `2-BWA_map/`

* **Purpose:** Secondary alignment of unmapped reads using BWA to verify genome origin.
* `*.header`: BAM header.
* `*_unmapped.sort.bam`: Sorted BAM file of previously unmapped reads aligned to the genome.

#### `3-Combine_bam/`

* **Purpose:** Merge HISAT2 and BWA BAMs and apply GATK pre-processing.
* Includes deduplication, base quality score recalibration (BQSR), and GATK-compatible formatting.
* `*_combine_*.bam`: BAMs at different pre-processing stages.
* `*.log`, `*.metrics`: GATK processing logs and QC metrics.

---

### `callediting` and `calltargets` steps

#### `4-Var_calling/`

* **Purpose:** Variant calling using GATK HaplotypeCaller.
* `*_HaplotypeCaller_*.vcf`: Per-chromosome and merged VCF files with SNVs.
* `*.log`, `*.json`: Process logs and input parameters.

#### `5-Var_filter/`

* **Purpose:** Final VCFs with known SNPs removed.
* `*_deAllSNP_*.vcf`: SNP-filtered VCFs used for editing site validation.

#### `6-Edit_calling/`

* **Purpose:** Final annotated, summarized RNA editing sites and differential editing analysis results.
* `Anno/`: Gene annotation BEDs (CDS, UTRs, exons, etc.).
* `Prefix/`:
  * `*_Edit*.bed`: Annotated RNA editing calls.
  * `*_GE_RPE_DATA.tsv`: Editing site-level expression and statistics.
* `PrefixG3BP1/`: Differential editing analysis results in transcript/genefull level.

---

### `SAILOR` and `FLARE` steps

#### `3.1-Split_strand/`

* **Purpose:** Splitting reads by strand for strand-specific analysis.
* `*.bam`, `*.bai`: Strand-separated BAM and index files.
* `*.sorted.bw`: BigWig files for coverage (forward/reverse strand).

#### `4.1-Redi_sailor/`

* **Purpose:** RNA editing site detection using REDItools/Sailor pipeline.
* `*_parallel_table*.txt`: Intermediate editing detection results.
* `*.coverage/`: Coverage files per chromosome for forward/reverse strands.

#### `5.1-Var_filter/`

* **Purpose:** Annotation of variants with known SNP datasets (dbSNP, 1000 Genomes, EVS).
* `*_dbSNP*.txt`: Filtered results by various SNP databases for forward/reverse reads.

#### `6.1-Edit_bedgraphs/`

* **Purpose:** Generate BED and bedGraph files of editing sites.
* `*.edit_fraction.bedgraph`: RNA editing site fractions.
* `*.snpfiltered.ranked*.bed`: Ranked editing sites filtered from SNPs (A-to-G / C-to-U).

#### `6.2-Edit_bigwig/`

* **Purpose:** Generate BigWig files for genome browser visualization.
* `*.sorted.bw`: Normalized RNA editing signal tracks (A>G, C>T).

#### `7-FLARE_for_MAPIT/`

* **Purpose:** Functional clustering of editing sites using FLARE.
* `confident_clusters/`: High-confidence editing clusters (AG/CT).
* `G3BP1/`: Subdirectories containing intermediate files from FLARE pipeline:

  * `editc_outputs/`: Cluster-level editing scores.
  * `FLARE/`: Final clusters with FDR filtering.
  * `bedgraphs/`: Cluster intensity scores in bedGraph format.
  * `bash_scripts/`, `outs/`: Execution scripts and outputs.