
Time_m_path=/data/rnomics6/fuzhican/project/RNA_editing_Time_flow_fine_tune
ABE_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
ABE_dep_path=${ABE_m_path}/dep_files
user_bin=/picb/rnomics4/rotation/fuzhican/bin
export PATH=/picb/rnomics4/rotation/fuzhican/bin:$PATH

annovar_annotation_summary_table_Para(){  
    ###### Thu May 21 10:08:33 CST 2020  
    ref_genome=$1;summary_file=$2  
    [[ "$ref_genome" == "hg38" ]]&& species=Human  
    [[ "$ref_genome" == "hg19" ]]&& species=Human  
    [[ "$ref_genome" == "rheMac10" ]]&& species=Macaque  
    [[ "$ref_genome" == "mm10" ]]&& species=Mouse  
    [[ "$ref_genome" == "Bter1" ]]&& species=Bee  
    [[ "$ref_genome" == "galGal6" ]]&& species=Chicken  
    [[ "$ref_genome" == "ce11" ]]&& species=CaenorhabditisElegans  
     
    test "$species" == "" && species=Human  

    if [[ "$ref_genome" == "hg38" ]]||[[ "$ref_genome" == "hg19" ]]||[[ "$ref_genome" == "rheMac10" ]]||[[ "$ref_genome" == "mm10" ]];then  
    echo $ref_genome |grep -q Refall && buildver=${ref_genome} || buildver=${ref_genome}Refall3  
    elif [[ "$ref_genome" == "Bter1" ]]||[[ "$ref_genome" == "galGal6" ]]||[[ "$ref_genome" == "ce11" ]];then  
    echo $ref_genome |grep -q Refall && buildver=${ref_genome} || buildver=${ref_genome}Refall1  
    fi  
    ANNOVAR_db=/picb/rnomics4/rotation/fuzhican/software/annovar/${species}Db_ref_all/  
    summary_path=`dirname $summary_file`  
    summary_file_basename=`basename $summary_file|cut -f1 -d.`  

    ANNOVAR_outpath=${summary_path}/ANNOVAR_splicing_threshold_3  
    test -d $ANNOVAR_outpath || mkdir -p $ANNOVAR_outpath  
    avi_input=${ANNOVAR_outpath}/${summary_file_basename}.avinput  
    ANNOVAR_out=${ANNOVAR_outpath}/${summary_file_basename}_ANNOVAR_annatation_${buildver}  
    
    awk 'BEGIN{OFS="\t"}$0 !~/_ER_summary/||$0 !~/^\t/{split($1,a,":");print a[1],a[2],a[2],"A","G","het","50.0","10"}' $summary_file >${avi_input}  

    echo "annotate_variation.pl --geneanno  --exonsort --exonicsplicing --splicing_threshold 3 --dbtype refGene --buildver $buildver --outfile ${ANNOVAR_out} ${avi_input} ${ANNOVAR_db}  "  
    annotate_variation.pl --geneanno  --exonsort --exonicsplicing --splicing_threshold 3 --dbtype refGene --buildver $buildver --outfile ${ANNOVAR_out} ${avi_input} ${ANNOVAR_db}  

}  

Annotate_Genome_Location_all_tissue_Para(){  
    ###### Sun Jul 12 21:12:46 CST 2020  
    method=site  
    ref_genome=$1;input_file=$2  
    output_file=`echo $input_file|cut -d. -f1`.Genome_Location  
    Annotate_Genome_Location $method $input_file $ref_genome $output_file  
}  
Annotate_Genome_Location(){
    ###### Fri Jul 3 14:50:02 CST 2020
    #plus_sh_workscript.sh Annotate_Genome_Location $method $input_file $ref_genome $output_file
    method=$1
    #method=bmc
    dep_path=$ABE_dep_path
    input_file=$2
    ref_genome=$3

    test -z $4 && output_name=${input_file}.Annotate_Genome_Location || output_name=$4 
    bmc_result_file=$input_file
    editing_sites=`dirname $bmc_result_file`
    outname=`basename $editing_sites`
    tag=_recal
    
    avinput_name=${input_file}.avinput
    
    if [ "$method" == "gatk" ];then
        vcf_name=$input_file
        convert2annovar.pl -format vcf4 ${vcf_name} > ${avinput_name}
    elif [ "$method" == "bmc" ];then
        awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),a[1]":"a[2],$3,$16,$NF}' $bmc_result_file >${avinput_name}
    elif [ "$method" == "site" ];then
        awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],"0","0"}' $bmc_result_file >${avinput_name}
    fi
    filter_region_type=Genome_Location_annotation
    out_1baseBed=${input_file}_${filter_region_type}_contain
    annotate_variation.pl ${avinput_name} ${dep_path}/${ref_genome} --buildver ${ref_genome} -bedfile Genome_Location_annotation.bed -dbtype bed -regionanno -out ${out_1baseBed} -colsWanted 7
    
    if [ "$method" == "gatk" ];then
    awk -F"Name=" '{print $2}' ${out_1baseBed}.${ref_genome}_bed |awk 'BEGIN{OFS="\t"}{if ($1 ~/Alu,/){$1="Alu"};print $2":"$3,$5">"$6,$1}' >$output_name 
    elif [ "$method" == "bmc" ];then
    awk -F"Name=" '{print $2}' ${out_1baseBed}.${ref_genome}_bed |awk 'BEGIN{OFS="\t"}{if ($1 ~/Alu,/){$1="Alu"};print $7,$5,$6,$8,$9,$10,$1}' >$output_name 
    elif [ "$method" == "site" ];then
    awk -F"Name=" '{print $2}' ${out_1baseBed}.${ref_genome}_bed |awk 'BEGIN{OFS="\t"}{if ($1 ~/Alu,/){$1="Alu"};print $2":"$3,".",$1}' >$output_name 
    
    fi
    echo $output_name
    
    

}
extract_context_seq_Para(){  
    ref_genome=$1;input_file=$2  
    ref_genome_path=${ABE_m_path}/dep_files/${ref_genome}/mapping_tools_index/${ref_genome}_all.fa  
    output_file=`echo $input_file|cut -d. -f1`.Extract_context_seq  
    internal_bed=`echo $input_file|cut -d. -f1`.bed  
    echo $internal_bed  
    awk 'BEGIN{OFS="\t"}NR!=1{split($1,a,":");if (a[2]-11>=0){ss=a[2]-11}else{ss=0}print a[1],ss,a[2]+10,$1}' $input_file >${internal_bed}  
    echo "${user_bin}/bedtools getfasta -fi $ref_genome_path -bed $internal_bed -name -tab >$output_file"
    ${user_bin}/bedtools getfasta -fi $ref_genome_path -bed $internal_bed -name -tab >$output_file 
    
}  
Annotate_exon_intron_interval(){
    ###### Fri Jul 3 14:50:02 CST 2020
    #plus_sh_workscript.sh Annotate_exon_intron_interval $method $input_file $ref_genome $output_file
    method=$1
    #method=bmc
    dep_path=$ABE_dep_path
    input_file=$2
    ref_genome=$3
    test -z $4 && output_name=${input_file}.Annotate_exon_intron_interval || output_name=$4 
    bmc_result_file=$input_file
    editing_sites=`dirname $bmc_result_file`
    outname=`basename $editing_sites`
    tag=_recal
    
    avinput_name=${input_file}.avinput
    
    if [ "$method" == "gatk" ];then
        vcf_name=$input_file
        convert2annovar.pl -format vcf4 ${vcf_name} > ${avinput_name}
    elif [ "$method" == "bmc" ];then
        awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],toupper($2),toupper($17),a[1]":"a[2],$3,$16,$NF}' $bmc_result_file >${avinput_name}
    elif [ "$method" == "site" ];then
        awk 'BEGIN{OFS="\t"}{split($1,a,":");print a[1],a[2],a[2],"0","0"}' $bmc_result_file >${avinput_name}
    fi
    filter_region_type=Exon_Intron_interval_annotation
    out_1baseBed=${input_file}_${filter_region_type}_contain

    dbpath=${dep_path}/${ref_genome}/bed_file
    [[ "$ref_genome" == "hg38" ]] && dbpath=${dep_path}/${ref_genome}/bed_file/170823_annotation
    ${user_bin}/annotate_variation.pl ${avinput_name} $dbpath  --buildver ${ref_genome} -bedfile ref_all_exon_intron_interval_GeneName.bed -dbtype bed -regionanno -out ${out_1baseBed} -colsWanted 4 ###### Mon Mar 1 23:58:51 CST 2021
    if [ "$method" == "gatk" ];then
    awk -F"Name=" '{print $2}' ${out_1baseBed}.${ref_genome}_bed |awk 'BEGIN{OFS="\t"}{if ($1 ~/Alu,/){$1="Alu"};print $2":"$3,$5">"$6,$1}' >$output_name 
    elif [ "$method" == "bmc" ];then
    awk -F"Name=" '{print $2}' ${out_1baseBed}.${ref_genome}_bed |awk 'BEGIN{OFS="\t"}{if ($1 ~/Alu,/){$1="Alu"};print $7,$5,$6,$8,$9,$10,$1}' >$output_name 
    elif [ "$method" == "site" ];then
    awk -F"Name=" '{print $2}' ${out_1baseBed}.${ref_genome}_bed |awk 'BEGIN{OFS="\t"}{if ($1 ~/Alu,/){$1="Alu"};print $2":"$3,".",$1}' >$output_name 
    fi
    echo $output_name

}
Annotate_Exon_Intron_interval_all_tissue_Para(){  
    ###### Sun Jul 12 21:12:46 CST 2020  
    method=site  
    ref_genome=$1;input_file=$2  
    output_file=`echo $input_file|cut -d. -f1`.Exon_Intron_interval  
    Annotate_exon_intron_interval $method $input_file $ref_genome $output_file  

}  
Annotate_Annotate_common_SNP_all_tissue_Para(){  
    ###### Sun Jul 12 21:12:46 CST 2020  
    method=site  
    ref_genome=$1;input_file=$2  
    output_file=`echo $input_file|cut -d. -f1`.Annotate_common_SNP  
    SNP_vcf=${ABE_m_path}/dep_files/${ref_genome}/NCBI_dbSNP_b151_common_${ref_genome}.vcf.gz  
    Annotate_vcf_core $method $input_file $SNP_vcf $output_file  
    echo $output_file  
}  
Annotate_vcf_core(){  
    method=$1  
    input_file=$2  
    SNP_vcf=$3  
    output_file=$4  
    token=`bcftools view -H $SNP_vcf|head|grep -q "^chr" || echo "chr" `  
    if [ "$method" == "site" ];then  
    awk 'FILENAME==ARGV[1]{array["'$token'"$1":"$2]++}FILENAME==ARGV[2]&&$1~/^chr/{if (!array[$1]){print $1"\t"1}else{print $1"\t"0}}' <(bcftools view -H $SNP_vcf) $input_file >$output_file  
    fi  
}  
"$@"
