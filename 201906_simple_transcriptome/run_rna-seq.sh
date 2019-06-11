#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe make 12

NDEX="~/db/Zea_mays.AGPv4.dna.toplevel.hisat2/Zea_mays.AGPv4.dna.toplevel.hisat2"
conda activate base

ls *R1.fq.gz | while read i
do
time hisat2 -p 10 --dta -x $INDEX -1 $i -2 ${i%R1*}R2.fq.gz | samtools view -Sb -q20 - | samtools sort - -o ${i%R1*}sorted.uniq.bam
time samtools index ${i%R1*}sorted.uniq.bam &

#cufflinks -p 12 -o ${i%R1*} -u -G /NAS7/home/gaoxiang/db/ann/Zea_mays.B73_RefGen_v4.43.gff3 ${i%R1*}sorted.uniq.bam

## ues stringtie instead of cufflinks, because cufflinks has some bugs
time ~/softwares/stringtie-1.3.6.Linux_x86_64/stringtie ${i%R1*}sorted.uniq.bam -G ~/db/ann/Zea_mays.B73_RefGen_v4.43.gff3 -A ${i%.R1*}/${i%R1*}gene_abund.tab -o ${i%.R1*}/stringtie.gtf -p 12 -e

done

