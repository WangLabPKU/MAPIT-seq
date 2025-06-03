import sys
import pysam
import re
insam=sys.argv[1]
outsam=sys.argv[2]

mapq_threshold=0
n_primary=0
n_not_xa=0
n_not_xa_mapq_gt_20=0
bam_f=insam
out_bam=outsam
bam = pysam.AlignmentFile(bam_f, 'rb')
bam_unique_mis6=pysam.AlignmentFile(out_bam, "w", template=bam)
for read in bam:
    if read.is_secondary:  # not the primary alignment
        continue
    n_primary+=1
    #if read.has_tag('MD'):break
    if read.has_tag('XA'):continue
    n_not_xa+=1
    if read.mapping_quality <= mapq_threshold: continue
    n_not_xa_mapq_gt_20+=1
    md_tag = read.get_tag('MD').split()[0]
    mis=len(re.findall(r'(\d+[A|T|C|G])',md_tag))
    if mis/read.query_alignment_length > 0.04: continue
    bam_unique_mis6.write(read)
bam_unique_mis6.close()
bam.close()
print("n_primary:",n_primary)
print("n_not_xa:",n_not_xa)
print("n_not_xa_mapq_gt_20:",n_not_xa_mapq_gt_20)

