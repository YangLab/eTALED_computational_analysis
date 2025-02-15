
/data/bioind1/fengh/software/anaconda3/envs/BEIDOU_env/bin/BEIDOU/BEIDOU -1 ${sample}_paired_R1.clean.fastq.gz -2 ${sample}_paired_R2.clean.fastq.gz -o `pwd`/${sample}_beidou_output -n ${sample}_beidou_output -g hg38 -t 8 -c /data/bioind1/gaobaoqing/DNA_Cas9_ChenJia/DNA_Cas9_ChenJia_20240318_WGS/mapping/bwa_men_sh/BEIDOU_config_hg38 -f SNV -p True
