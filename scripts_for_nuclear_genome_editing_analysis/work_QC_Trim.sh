
# ${sample}_R1.fastq.gz and ${sample}_R2.fastq.gz are raw WGS data for one sample


# FastQC:quality control
fastqc ${sample}_R1.fastq.gz 
fastqc ${sample}_R2.fastq.gz

# Trimmomatic PE
trimmomatic PE -phred33 -threads 64 ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz ${sample}_paired_R1.clean.fastq.gz ${sample}_unpaired_R1.clean.fastq.gz ${sample}_paired_R2.clean.fastq.gz ${sample}_unpaired_R2.clean.fastq.gz ILLUMINACLIP:trim_adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
