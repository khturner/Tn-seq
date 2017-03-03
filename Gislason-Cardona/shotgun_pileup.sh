cd ~/cenocepacia-s3
aws s3 cp s3://cenocepacia-tnseq/K56-2-1383382_S1_L001_R1_001.fastq.gz ./
aws s3 cp s3://cenocepacia-tnseq/K56-2-1383382_S1_L001_R2_001.fastq.gz ./

# Doh they're already trimmed
# cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o shotgun.r1.fastq -p shotgun.r2.fastq K56-2-1383382_S1_L001_R1_001.fastq.gz K56-2-1383382_S1_L001_R2_001.fastq.gz 

bowtie2 -p 4 -x k56.fasta -1 K56-2-1383382_S1_L001_R1_001.fastq -2 K56-2-1383382_S1_L001_R2_001.fastq -S shotgun.sam
samtools view -b -S shotgun.sam > shotgun.bam
samtools sort shotgun.bam shotgun.sorted
samtools index shotgun.sorted.bam
samtools mpileup -f k56.fasta shotgun.sorted.bam > shotgun.pileup

