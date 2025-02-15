

# Raw fastq files are saved in directory 0-0_fastq

# FastQC:quality control
mkdir 1-0_fastqc
fastqc ../0-0_fastq/*.fq.gz -o .

# Trimmomatic PE
mkdir 2-0_trimmomatic
parallel -j 1 --rpl ',, s/.*\/(.*).R1.fq.gz/$1/' 'trimmomatic PE -threads 10 {} {= s/.R1./.R2./ =} ,,.R1.trimmed.fq.gz  ,,.R1.unpaired.fq.gz  ,,.R2.trimmed.fq.gz ,,.R2.unpaired.fq.gz ILLUMINACLIP:/data/bioind1/xingma/Projects/BEhuman/CFBI_test/TruSeq3-PE-2.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &> ,,.log' ::: ../0-0_fastq/*.R1.fq.gz

# reads sampling
mkdir 3-0_seqtk
parallel -j 1 --plus 'seqtk sample -s100 {} 30000000 |pigz -p 10 > {/...}.37M.fq.gz' ::: ../2-0_trimmomatic/*.trimmed.fq.gz
