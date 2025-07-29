set -e

# 1. mapping to reference genome
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/N_1_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/N_1_val_2.fq.gz --rna-strandness FR -n N -r 1 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/N_2_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/N_2_val_2.fq.gz --rna-strandness FR -n N -r 2 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YI0_1_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YI0_1_val_2.fq.gz --rna-strandness FR -n YI0 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YI0_2_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YI0_2_val_2.fq.gz --rna-strandness FR -n YI0 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YX_1_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YX_1_val_2.fq.gz --rna-strandness FR -n YX -r 1 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YX_2_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YX_2_val_2.fq.gz --rna-strandness FR -n YX -r 2 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YF0_1_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YF0_1_val_2.fq.gz --rna-strandness FR -n YF0 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YF0_2_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YF0_2_val_2.fq.gz --rna-strandness FR -n YF0 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YI1_1_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YI1_1_val_2.fq.gz --rna-strandness FR -n YI1 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YI1_2_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YI1_2_val_2.fq.gz --rna-strandness FR -n YI1 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YF1_1_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YF1_1_val_2.fq.gz --rna-strandness FR -n YF1 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result -t 40
Mapit mapping -v GRCh38 --fq /home/gangx/data/CQX/230301/clean/YF1_2_val_1.fq.gz --fq2 /home/gangx/data/CQX/230301/clean/YF1_2_val_2.fq.gz --rna-strandness FR -n YF1 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result -t 40

# 2.finetuning alignments
Mapit finetuning -v GRCh38 -n N -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n N -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YI0 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YI0 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YX -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YX -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
wait

Mapit finetuning -v GRCh38 -n YF0 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YF0 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YI1 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YI1 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YF1 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit finetuning -v GRCh38 -n YF1 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
wait

# 3. call RNA editing
Mapit callediting -v GRCh38 --sampleList N,YI0,YX,YF0,YI1,YF1 -o /home/gangx/data/CQX/230301/Mapit_result --prefix FL -e Both -t 50

# 4. call RBP targets
Mapit calltargets -v GRCh38 -i /home/gangx/data/CQX/230301/Mapit_result/6-Edit_calling/FL/FL_Edit_GE_RPE_DATA.tsv -l transcript --treatName YF0 --controlName N   -o /home/gangx/data/CQX/230301/Mapit_result -p 0.05 &
Mapit calltargets -v GRCh38 -i /home/gangx/data/CQX/230301/Mapit_result/6-Edit_calling/FL/FL_Edit_GE_RPE_DATA.tsv -l transcript --treatName YF0 --controlName YX  -o /home/gangx/data/CQX/230301/Mapit_result -p 0.05 &
Mapit calltargets -v GRCh38 -i /home/gangx/data/CQX/230301/Mapit_result/6-Edit_calling/FL/FL_Edit_GE_RPE_DATA.tsv -l transcript --treatName YF0 --controlName YI0 -o /home/gangx/data/CQX/230301/Mapit_result -p 0.05 &
Mapit calltargets -v GRCh38 -i /home/gangx/data/CQX/230301/Mapit_result/6-Edit_calling/FL/FL_Edit_GE_RPE_DATA.tsv -l transcript --treatName YF1 --controlName YI1 -o /home/gangx/data/CQX/230301/Mapit_result -p 0.05 &
wait

# 5. prepare files for SAILOR workflow
Mapit prepare -v GRCh38 -n N -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n N -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YI0 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YI0 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YX -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YX -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
wait

Mapit prepare -v GRCh38 -n YF0 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YF0 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YI1 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YI1 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YF1 -r 1 -o /home/gangx/data/CQX/230301/Mapit_result &
Mapit prepare -v GRCh38 -n YF1 -r 2 -o /home/gangx/data/CQX/230301/Mapit_result &
wait

# 6. run SAILOR and FLARE identidy edit clusters and high-confidence MAPIT edit clusters
regions=/home/gangx/data/Reference/GRCh38/GRCh38_FLARE_regions

for sample in N YI0 YX YF0 YI1 YF1
do
for rep in 1 2
do
Mapit SAILOR -v GRCh38 -n $sample -r $rep -c 10 -o /home/gangx/data/CQX/230301/Mapit_result -t 20
Mapit FLARE -v GRCh38 -e CT -n $sample -r $rep --regions ${regions} -o /home/gangx/data/CQX/230301/Mapit_result -t 20
Mapit FLARE -v GRCh38 -e AG -n $sample -r $rep --regions ${regions} -o /home/gangx/data/CQX/230301/Mapit_result -t 20
done
Mapit hc_cluster -v GRCh38 -n $sample -o /home/gangx/data/CQX/230301/Mapit_result -s 150
done

# 7. de nove call RBP binding motifs
hc_cluster_path=${outpath}/7-FLARE_for_MAPIT/confident_clusters
hc_cluster_bed=${hc_cluster_path}/${samplename}/${samplename}_${slop}.bed
output_motif=${hc_cluster_path}/${samplename}/${samplename}.motif

findMotifsGenome.pl $hc_cluster_bed hg38 $output_motif -noknown -rna -size given -len 5,6,7,8 -p ${threads}
