#!/usr/bin/env bash
touch $0\.arg 2>/dev/null && args_path=`dirname $0` ||args_path=`pwd`
set -eo pipefail

echo -e `date` dir:$PWD  command:$0 "$*" >>$args_path/`basename $0`\.arg  ###### Sun Dec 8 10:27:47 CST 2019 fzc introduce $args_path
######!!!!!!!!!Warining!!!!!!!!!!############
##USE bash instead of sh to run this shell script!
######!!!!!!!!!Warining!!!!!!!!!!############
###############Usage############
###paired end:
#bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller -1 ${path_of_fastq1} -2 ${path_of_fastq2}  -o ${output_dir} -n ${name} -g ${ref_genome} -t ${maximum_threads}
###single end:
#bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller -s ${path_of_fastq} -o ${output_dir} -n ${name} -g ${ref_genome} -t ${maximum_threads}

#ref_genome="hg19" or "hg38"

#bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller -1 /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_RNA_20191016/backup/fqtrimmed/TET1_ctrl_1st_HEK293FT_R1_trimmed.fq -2 /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_RNA_20191016/backup/fqtrimmed/TET1_ctrl_1st_HEK293FT_R2_trimmed.fq  -o /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/tmp_dir -n TET1_ctrl_1st_HEK293FT -g hg19 -t 10


################################



PATH=/picb/rnomics4/rotation/fuzhican/software/conda2_4_7/bin:/picb/rnomics4/rotation/fuzhican/software/conda2_4_7/condabin:/picb/rnomics4/rotation/fuzhican/software/conda/bin:/picb/rnomics4/rotation/fuzhican/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/picb/intprog/bin:/picb/extprog/bin:/picb/extprog/bin/linux-x86:/picb/extprog/bin/Linux-amd64:/picb/ports/Linux26-x86-64-prefix/usr/bin:/picb/ports/Linux26-x86-64-prefix/bin:/picb/ports/Linux26-x86-64/bin:/picb/ports/Linux26-x86-64/sbin:/picb/ports/Linux26-x86-64-nix/usr/bin

user_bin=/picb/rnomics4/rotation/fuzhican/bin

Time_m_path=/data/rnomics6/fuzhican/project/RNA_editing_Time_flow_fine_tune
ABE_m_path=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target
dep_path=${ABE_m_path}/dep_files
sh_path=/data/rnomics10/yuanguohua/20220925_RNAtBE/RADAR_scripts

picard=/picb/rnomics4/rotation/fuzhican/bin/picard.jar
gatk4=/picb/rnomics4/rotation/fuzhican/bin/gatk4
samtools1=/picb/rnomics4/rotation/fuzhican/software/conda/envs/circ/bin/samtools
seqtk='/picb/rnomics4/rotation/fuzhican/software/conda/envs/clear/bin/seqtk'
bwa=/picb/rnomics4/rotation/fuzhican/software/conda/envs/variant/bin/bwa

fastqc=/picb/rnomics4/rotation/fuzhican/bin/fastqc

init_1(){
    ref_genome=$1
    mapping_tools_index_path=${dep_path}/${ref_genome}/mapping_tools_index
    rDNA_index_bwa_mem=${mapping_tools_index_path}/rDNA.fa
    # rDNA_index_bwa_mem=${mapping_tools_index_path}/RNA45SN4.fa
    ref_genome_path=${mapping_tools_index_path}/${ref_genome}_all.fa

}
RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller(){
    ###### Wed Sep 25 20:00:00 CST 2019
    ####################init#################################################
    ##Usage:
    ###paired end:
    #bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller -1 ${path_of_fastq1} -2 ${path_of_fastq2}  -o ${output_dir} -n ${name} -g ${ref_genome} -t ${maximum_threads}
    ###single end:
    #bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller -s ${path_of_fastq} -o ${output_dir} -n ${name} -g ${ref_genome} -t ${maximum_threads}

    #ref_genome="hg19" or "hg38"

    #bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller -1 /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_RNA_20191016/backup/fqtrimmed/TET1_ctrl_1st_HEK293FT_R1_trimmed.fq -2 /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_RNA_20191016/backup/fqtrimmed/TET1_ctrl_1st_HEK293FT_R2_trimmed.fq  -o /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/tmp_dir -n TET1_ctrl_1st_HEK293FT -g hg19 -t 10
    ####################init END#################################################
    while getopts :1:2:s:o:n:g:t: ARGS  
    do  
    case $ARGS in   
        1)  
            fq1=$OPTARG
            ;;  
        2)  
            fq2=$OPTARG  
            ;;
        s)  
            fq0=$OPTARG
            ;;  
        t)  
            threads=$OPTARG
            ;;  
        o)  
            m0_path=$OPTARG
            ;;  
        n)  
            name=$OPTARG
            ;; 
        g)  
            ref_genome=$OPTARG
            ;;  
        *)  
            echo "Unknown option: $ARGS"
            ;;
        \?)
        echo "Invalid option: -$OPTARG" 
        ;;
    esac
    done
    srr=$name
    init_1 $ref_genome
    tmp_work_path_0=$m0_path
    tmp_a="bb${fq0}bb"
    test "$tmp_a" == "bbbb" && {
        layout=paired
        fq_path=`dirname ${fq1}`
        }|| {
            layout=single
            fq_path=`dirname ${fq0}`
        }
    threads_tmp_a="bb${threads}bb"
    test "$threads_tmp_a" == "bbbb" && threads=10
    m0_path=$tmp_work_path_0
    HISAT2_mapping_path="${tmp_work_path_0}/${ref_genome}/1-1-0HISAT_map"
    BWA_mapping_path="${tmp_work_path_0}/${ref_genome}/2-0-0BWA_map"
    
    fine_tune_path="${tmp_work_path_0}/${ref_genome}/3-0-0Combine_bam"
    ###1. mapping, HISAT2_2mismatch_following_BWA_6mismatch_mapping
    
    HISAT2_2mismatch_following_BWA_6mismatch_mapping $srr $layout $fq_path $HISAT2_mapping_path $BWA_mapping_path $fine_tune_path ${threads} ###### Thu Nov 21 20:28:22 CST 2019 fzc  support set threads number when running mapping
    
    ###2. sam_fine_tune, Markduplicate, BQSR
    bam_file=${fine_tune_path}/${srr}_combine.bam
    sam_fine_tune $bam_file $fine_tune_path

    ###3. GATK_HaplotypeCaller Calling
    bam_file_prefix=${srr}_combine
    bam_file=$fine_tune_path/${bam_file_prefix}_rgadd_dedupped_split_recal.bam

    gatk_Variant_filter $bam_file "_2pass" $m0_path ${srr} $ref_genome
    
}
RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_Samtools_mpileup(){
    #RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_Samtools_mpileup $workspace $fq_path $srr $layout

    ###### Tue Dec 3 14:35:27 CST 2019
    ####################init#################################################
    ##Usage:
    ###paired end:
    #bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_Samtools_mpileup -1 ${path_of_fastq1} -2 ${path_of_fastq2}  -o ${output_dir} -n ${name} -g ${ref_genome} -t ${maximum_threads}
    ###single end:
    #bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_Samtools_mpileup -s ${path_of_fastq} -o ${output_dir} -n ${name} -g ${ref_genome} -t ${maximum_threads}

    #ref_genome="hg19" or "hg38"

    #bash GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_Samtools_mpileup -1 /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_RNA_20191016/backup/fqtrimmed/TET1_ctrl_1st_HEK293FT_R1_trimmed.fq -2 /data/rnomics6/xuew/Human/BE_WangLiJie/SBE_RNA_20191016/backup/fqtrimmed/TET1_ctrl_1st_HEK293FT_R2_trimmed.fq  -o /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/tmp_dir -n TET1_ctrl_1st_HEK293FT -g hg19 -t 10
    ####################init END#################################################
    while getopts :1:2:s:o:n:g:t:p: ARGS  
    do  
    case $ARGS in   
        1)  
            fq1=$OPTARG
            ;;  
        2)  
            fq2=$OPTARG  
            ;;
        s)  
            fq0=$OPTARG
            ;;  
        t)  
            threads=$OPTARG
            ;;  
        o)  
            m0_path=$OPTARG
            ;;  
        n)  
            name=$OPTARG
            ;; 
        g)  
            ref_genome=$OPTARG
            ;;  
        p) ###### Sat Feb 22 15:48:50 CST 2020 fzc; specify run step when fix 
            step=$OPTARG
            ;;  
        *)  
            echo "Unknown option: $ARGS"
            ;;
        \?)
        echo "Invalid option: -$OPTARG" 
        ;;
    esac
    done
    srr=$name
    init_1 $ref_genome
    tmp_work_path_0=$m0_path
    tmp_a="bb${fq0}bb"
    test "$tmp_a" == "bbbb" && {
        layout=paired
        fq_path=`dirname ${fq1}`
        }|| {
            layout=single
            fq_path=`dirname ${fq0}`
        }
    threads_tmp_a="bb${threads}bb"
    test "$threads_tmp_a" == "bbbb" && threads=10
    step_tmp_a="bb${step}bb"
    test "$step_tmp_a" == "bbbb" && step=0123
    m0_path=$tmp_work_path_0
    rRNAdeplete_path="${tmp_work_path_0}/${ref_genome}/0-0-0remove_rRNA"
    HISAT2_mapping_path="${tmp_work_path_0}/${ref_genome}/1-1-0HISAT_map"
    BWA_mapping_path="${tmp_work_path_0}/${ref_genome}/2-0-0BWA_map"
    
    
    fine_tune_path="${tmp_work_path_0}/${ref_genome}/3-0-0Combine_bam"
    
    echo "#####[`date`]${srr} BEGIN..."

    Y_or_N=$(echo "$step"|grep -q "0" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###0. remove rRNA, remove_rRNA_mapped_reads
    echo "#####[`date`]BEGIN remove_rRNA_mapped_reads"
    if [ "$layout" == "paired" ];then
    out_fq1=${rRNAdeplete_path}/${name}_R1.fastq.gz
	out_fq2=${rRNAdeplete_path}/${name}_R2.fastq.gz
    remove_rRNA_mapped_reads $fq1 $fq2 $out_fq1 $out_fq2  ###### Wed Dec 23 17:07:21 CST 2020  
    fq1=$out_fq1
    fq2=$out_fq2
    elif [ "$layout" == "single" ];then
    out_fq0=${rRNAdeplete_path}/${name}.fastq.gz
    remove_rRNA_mapped_reads $fq0 $out_fq0
    fq0=$out_fq0
    fi
    fi



    Y_or_N=$(echo "$step"|grep -q "1" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###1. mapping, HISAT2_2mismatch_following_BWA_6mismatch_mapping
    echo "#####[`date`]BEGIN HISAT2_2mismatch_following_BWA_6mismatch_mapping"
    HISAT2_2mismatch_following_BWA_6mismatch_mapping $srr $layout $HISAT2_mapping_path $BWA_mapping_path $fine_tune_path  ###### Tue Dec 3 14:35:36 CST 2019 fzc  support set threads number when running mapping
    fi


    Y_or_N=$(echo "$step"|grep -q "2" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###2. sam_fine_tune, Markduplicate, BQSR
    echo "######[`date`]BEGIN sam_fine_tune"
    bam_file=${fine_tune_path}/${srr}_combine.bam  ###### Sat Dec 21 06:25:53 CST 2019 fzc; change ${combine_bam} to ${fine_tune_path}
    sam_fine_tune $bam_file $fine_tune_path 
    fi


    Y_or_N=$(echo "$step"|grep -q "3" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###3. Samtools_mpileup Calling
    echo "#####[`date`]BEGIN bmc_Variant_filter"
    bam_file_prefix=${srr}_combine
    bam_file=$fine_tune_path/${bam_file_prefix}_rgadd_dedupped_split_recal.bam
    bmc_Variant_filter $bam_file "_recal" $m0_path ${srr} $ref_genome
    fi

    echo "#####[`date`]${srr} END"
}

RNA_Editing_Calling_Pipeline_STAR_followed_by_GATK_HaplotypeCaller(){
    #RNA_Editing_Calling_Pipeline_STAR_followed_by_GATK_HaplotypeCaller $workspace_2 $fq_path $srr $layout $max_readlength

    ###### Tue Oct 8 18:43:00 CST 2019
    ####################init#################################################
    ##Usage:
    #RNA_Editing_Calling_Pipeline_STAR_followed_by_GATK_HaplotypeCaller $tmp_work_path_0 $fq_path $srr $layout $max_readlength
    #Parameters example:
    #tmp_work_path_0="${ABE_m_path}/SameBam_TwoCaller/HISAT2_6mis_bam"
    #fq_path="${ABE_m_path}/dep_files/fastq"
    #srr="NT-Cas9_rep1"
    #layout="paired"
    ##
    #used for RNA_Editing_Calling_Pipeline_STAR_followed_by_GATK_HaplotypeCaller,\
    #assemble subfunctions:mapping, sam fine tune, caliing and variants filtering.

    ##required subfunctions:
    #STAR_10mis_mapping $srr $layout $fq_path $mapping_path $max_readlength
    #sam_fine_tune $bam_file $fine_tune_path
    #gatk_Variant_filter $bam_file "_2pass" $m0_path ${srr}
    ####################init END#################################################
    tmp_work_path_0=$1
    fq_path=$2
    srr=$3
    layout=$4
    max_readlength=$5


    m0_path=$tmp_work_path_0
    STAR_mapping_path="${tmp_work_path_0}/${ref_genome}/1AND2-0-0STAR_maping_TwoPass"
    mapping_path=$STAR_mapping_path
    
    fine_tune_path="${tmp_work_path_0}/${ref_genome}/3-0-0Sam_fine-tune"
    ###1. mapping, STAR_10mis_mapping
    
    STAR_10mis_mapping $srr $layout $fq_path $mapping_path $max_readlength
    
    ###2. sam_fine_tune, Markduplicate, BQSR
    bam_file=${mapping_path}/2-0-0STAR_mapping_SecondPass/${srr}_Aligned.out.bam
    bam_file_basename_prefix=$srr
    sam_fine_tune $bam_file $fine_tune_path $bam_file_basename_prefix

    ###3. GATK_HaplotypeCaller Calling
    bam_file_prefix=${bam_file_basename_prefix}
    bam_file=$fine_tune_path/${bam_file_prefix}_rgadd_dedupped_split_recal.bam
    gatk_Variant_filter $bam_file "_2pass" $m0_path ${srr}

}
STAR_10mis_mapping(){
    #STAR_10mis_mapping $srr $layout $fq_path $mapping_path $max_readlength
    ###### Tue Oct 8 19:12:45 CST 2019
    #STAR mapping paprameters from GATK RNA-seq calling pipeline(https://software.broadinstitute.org/gatk/documentation/article.php?id=3891) and "Transcriptome-wide off-target RNA editing induced by CRISPR-guided DNA base editors" GEO SOFT formatted family file
    ###### Tue Oct 8 18:55:08 CST 2019
    srr=$1
    layout=$2
    fq_path=$3
    mapping_path=$4
    max_readlength=$5

    STAR=/picb/rnomics4/rotation/fuzhican/bin/STAR   ##STAR_2.6.0c
    ref_genome=${ref_genome}
    star_index1=${dep_path}/${ref_genome}/star_index_${ref_genome}
    k_parameter=1
    test -e ${fq_path}/${srr}_1.fastq.gz && {
        fq_suffix_1=_1.fastq.gz
        fq_suffix_2=_2.fastq.gz
    }
    test -e ${fq_path}/${srr}_R1.fastq.gz && {
        fq_suffix_1=_R1.fastq.gz
        fq_suffix_2=_R2.fastq.gz
    }
    
    
    star_map1=${mapping_path}/1-0-0STAR_mapping_FirstPass
    star_map2=${mapping_path}/2-0-0STAR_mapping_SecondPass
    test -d $star_map1 ||mkdir -p $star_map1
    test -d $star_map2 ||mkdir -p $star_map2

    ###1. First pass mapping, in order to produce splice junction information

    test -d ${star_map1}/${srr}_STARtmp && rm -r ${star_map1}/${srr}_STARtmp
    $STAR --genomeDir $star_index1 --readFilesIn ${fq_path}/${srr}${fq_suffix_1} ${fq_path}/${srr}${fq_suffix_2} --outTmpDir ${star_map1}/${srr}_STARtmp --readFilesCommand zcat --runThreadN 12 --outFileNamePrefix ${star_map1}/${srr}_ --outSAMtype BAM Unsorted #--outFilterScoreMinOverLread 0 


    ###2. Second pass mapping, mapping with splice junction information produced by the first pass

    star_index2=${star_map2}/star_index_2pass_${srr}
    test -d ${star_index2}||mkdir ${star_index2}
    test -d ${star_map2}/${srr}_STARtmp && rm -r ${star_map2}/${srr}_STARtmp
    ## a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass
    $STAR --runMode genomeGenerate --genomeDir ${star_index2} --genomeFastaFiles /picb/rnomics1/database/Human/${ref_genome}/genome/${ref_genome}_all.fa --sjdbFileChrStartEnd ${star_map1}/${srr}_SJ.out.tab --outTmpDir ${star_map2}/${srr}_STARtmp --sjdbGTFfile ${dep_path}/${ref_genome}/${ref_genome}_annotation/ref_all.gtf --sjdbOverhang `expr $max_readlength - 1` --runThreadN 12 

    ##The resulting index is then used to produce the final alignments as follows
    test -d ${star_map2}/${srr}_STARtmp && rm -r ${star_map2}/${srr}_STARtmp
    $STAR --genomeDir ${star_index2} --readFilesIn ${fq_path}/${srr}${fq_suffix_1} ${fq_path}/${srr}${fq_suffix_2} --outTmpDir ${star_map2}/${srr}_STARtmp --readFilesCommand zcat  --outFilterMultimapNmax $k_parameter --runThreadN 12 --outFileNamePrefix ${star_map2}/${srr}_ --outSAMtype BAM Unsorted 
    ##Filnal mapped bam was deposited at ${star_map2}

}
remove_rRNA_mapped_reads(){
    outname=$name
    rRNAdeplete_path=$rRNAdeplete_path
    thread=$threads
    test -d $rRNAdeplete_path ||mkdir -p $rRNAdeplete_path
    if [ "$layout" == "paired" ];then
    in_fq1=$1
    in_fq2=$2
    out_fq1=$3
    out_fq2=$4
    ## remove rRNA
	# bwa  mem -t ${thread} ${rDNA_index_bwa_mem}  $in_fq1  $in_fq2 > ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam 
	#2>${rRNAdeplete_path}/log_BWA_rRNA_`date +%Y_%m_%d`.log
	# ${samtools1} view --threads ${thread} -bh -f 4 ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam > ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
	# ${samtools1} sort --threads ${thread} -n ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam  -o ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
	# ${samtools1} fastq --threads ${thread} -1 ${out_fq1} -2 ${out_fq2} -s ${rRNAdeplete_path}/${outname}_singleton.fq ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
	# rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
	# rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
    echo "bwa mem -t ${thread} ${rDNA_index_bwa_mem}  $in_fq1  $in_fq2|${samtools1} view --threads ${thread} -bh  -f 4 |${samtools1} sort --threads ${thread} -n  -o ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam 2>${rRNAdeplete_path}/log_BWA_derRNA_`date +%Y_%m_%d`.log"
    bwa mem -t ${thread} ${rDNA_index_bwa_mem}  $in_fq1  $in_fq2|${samtools1} view --threads ${thread} -bh  -f 4 |${samtools1} sort --threads ${thread} -n  -o ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam 2>${rRNAdeplete_path}/log_BWA_derRNA_`date +%Y_%m_%d`.log
    ${samtools1} fastq --threads ${thread} -1 ${out_fq1} -2 ${out_fq2} -s /dev/null ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
   	#fastqc_out=`dirname ${out_fq1}`/fastqc_out
    #test -d $fastqc_out ||mkdir  $fastqc_out
    #$fastqc -o ${fastqc_out} -d ~/tmp -q ${out_fq1} ${out_fq2}  

    # ls *_bwa_mapped_rRNA.sam|awk -F "_bwa_mapped_rRNA.sam" '{print $1}'|pll -j 2 "${samtools1} view --threads ${thread} -bh  -f 4 {}_bwa_mapped_rRNA.sam |${samtools1} sort --threads ${thread} -n  -o ${rRNAdeplete_path}/{}-rRNA_unmapped_sort.bam "&

	
    elif [ "$layout" == "single" ];then
    in_fq0=$1
    out_fq0=$2
    # bwa  mem -t ${thread} ${rDNA_index_bwa_mem}  $in_fq0  > ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam 
	# #2>${rRNAdeplete_path}/log_BWA_rRNA_`date +%Y_%m_%d`.log
    #     ${samtools1} view --threads ${thread} -bh -f 4 ${rRNAdeplete_path}/${outname}_bwa_mapped_rRNA.sam > ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
    #     ${samtools1} sort --threads ${thread} -n ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam  -o ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam
    #     ${samtools1} fastq --threads ${thread}  ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam | pigz -p 3 > $out_fq0
    #     rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped.bam
    #     rm ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam

        bwa mem -t ${thread} ${rDNA_index_bwa_mem}  $in_fq0|${samtools1} view --threads ${thread} -bh  -f 4 |${samtools1} sort --threads ${thread} -n  -o ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam 2>${rRNAdeplete_path}/log_BWA_derRNA_`date +%Y_%m_%d`.log
        ${samtools1} fastq --threads ${thread} -s ${rRNAdeplete_path}/${outname}_singleton.fq ${rRNAdeplete_path}/${outname}-rRNA_unmapped_sort.bam | pigz -p 3 > $out_fq0

        #fastqc_out=`dirname ${out_fq0}`
        #$fastqc -o ${fastqc_out} -d ~/tmp -q ${out_fq0} 
	fi

    
}
HISAT2_2mismatch_following_BWA_6mismatch_mapping(){
    #HISAT2_2mismatch_following_BWA_6mismatch_mapping $srr $layout $fq_path $HISAT2_mapping_path $BWA_mapping_path $combine_bam
    ###### Wed Sep 25 19:04:01 CST 2019
    #used for HISAT2 and BWA two pass mapping.
    #First, HISAT2 2 mismatches
    #HISAT2 unmapped reads are extracted and send to BWA 6 mismatches mapping.
    #for paired ends reads, mismatches are all count at read level not fragment level, for example, HISAT2 2 mismatches means read1's mismatches maximum count is 2,read2's mismatches maximum count is 2, read1 mismatches plus read2 mismatches can up to 4.

    srr=$1
    layout=$2
    HISAT2_mapping_path=$3
    BWA_mapping_path=$4
    combine_bam=$5

    test -d $HISAT2_mapping_path || mkdir -p $HISAT2_mapping_path
    test -d $BWA_mapping_path || mkdir -p $BWA_mapping_path
    test -d ${combine_bam} ||mkdir -p ${combine_bam}
    HISAT_map=$HISAT2_mapping_path
    bwa_map=${BWA_mapping_path}

    if [ "$layout" == "paired" ];then

        ### 1. HISAT2 2 mismatches mapping
        ##################need confirm: --rna-strandness RF update:###### Wed Sep 25 19:56:57 CST 2019 confirmed, 0.12 VS 0.88
        #hisat2  --rna-strandness RF --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${dep_path}/${ref_genome}/${ref_genome}_annotation/ref_all_spsites.txt --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p 10 -x /picb/rnomics1/database/Human/${ref_genome}/genome/${ref_genome}_all -1  ${fq_path}/${srr}${fq_suffix_1}  -2 ${fq_path}/${srr}${fq_suffix_2}  --un-conc-gz ${HISAT_map}/${srr}_un_conc_%.fastq.gz -S ${HISAT_map}/${srr}_HISAT2_mapped.sam  2>${HISAT_map}/log_HISAT2_2mismatch_${srr}_`date +%Y_%m_%d`.log 
        /picb/rnomics4/rotation/fuzhican/bin/hisat2  --rna-strandness RF --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${dep_path}/${ref_genome}/${ref_genome}_annotation/ref_all_spsites.txt --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p ${threads} -x ${mapping_tools_index_path}/${ref_genome}_all -1  $fq1  -2 $fq2  --un-conc-gz ${HISAT_map}/${srr}_un_conc_%.fastq.gz -S ${HISAT_map}/${srr}_HISAT2_mapped.sam  2>${HISAT_map}/log_HISAT2_2mismatch_${srr}_`date +%Y_%m_%d`.log ###### Wed Nov 20 08:06:39 CST 2019 fzc
        
        $samtools1 view -h -F 4 ${HISAT_map}/${srr}_HISAT2_mapped.sam|awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) print $0 } else print $0 " not have XM tag"}}'|awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag"  }}' >${HISAT_map}/${srr}_unique_mismatch2.sam &
        
        $samtools1 view -bS -f 4 -o ${HISAT_map}/${srr}_HISAT2_unmapped.bam ${HISAT_map}/${srr}_HISAT2_mapped.sam

        $samtools1 view -f 4 -S ${HISAT_map}/${srr}_HISAT2_mapped.sam |awk 'BEGIN{FS="\t"}{print $1}'|sort|uniq >${HISAT_map}/${srr}_hisat2_unmap.readid 

        $seqtk subseq ${HISAT_map}/${srr}_un_conc_1.fastq.gz ${HISAT_map}/${srr}_hisat2_unmap.readid |gzip > ${HISAT_map}/${srr}_unmapped_1.fastq.gz  &
        $seqtk subseq ${HISAT_map}/${srr}_un_conc_2.fastq.gz ${HISAT_map}/${srr}_hisat2_unmap.readid |gzip > ${HISAT_map}/${srr}_unmapped_2.fastq.gz  
        wait

        ### 2. BWA 6 mismatches mapping
        
        test -d $bwa_map||mkdir -p $bwa_map
        bwa mem -t ${threads}  -A 1 -B 4  ${ref_genome_path}  ${HISAT_map}/${srr}_unmapped_1.fastq.gz  ${HISAT_map}/${srr}_unmapped_2.fastq.gz > ${bwa_map}/${srr}_bwa_mapped.sam 2>${bwa_map}/log_BWA_6mismatch_`date +%Y_%m_%d`.log
    elif [ "$layout" == "single" ];then

        test -d ${HISAT_map}||mkdir -p ${HISAT_map}
        #hisat2 --secondary --no-temp-splicesite --known-splicesite-infile ${dep_path}/${ref_genome}/${ref_genome}_annotation/ref_all_spsites.txt --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p 10 -x /picb/rnomics1/database/Human/${ref_genome}/genome/${ref_genome}_all -U ${fq_path}/${srr}.fastq.gz -S ${HISAT_map}/${srr}_HISAT2_mapped.sam 2>${log_path}/${ref_genome}/bmc/HISAT2/log_hisat2_${srr}.log 

        /picb/rnomics4/rotation/fuzhican/bin/hisat2 --secondary --no-temp-splicesite --known-splicesite-infile ${dep_path}/${ref_genome}/${ref_genome}_annotation/ref_all_spsites.txt --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p ${threads} -x ${mapping_tools_index_path}/${ref_genome}_all -U $fq0 -S ${HISAT_map}/${srr}_HISAT2_mapped.sam 2>${HISAT_map}/log_HISAT2_2mismatch_${srr}_`date +%Y_%m_%d`.log ###### Tue Dec 17 19:58:15 CST 2019 fzc; change log path

        $samtools1 view -h -F 4 ${HISAT_map}/${srr}_HISAT2_mapped.sam|awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) print $0 } else print $0 " not have XM tag" }}'|awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag" }}' >${HISAT_map}/${srr}_unique_mismatch2.sam &
        $samtools1 view -bS -f 4 -o ${HISAT_map}/${srr}_HISAT2_unmapped.bam ${HISAT_map}/${srr}_HISAT2_mapped.sam 
        /picb/rnomics4/rotation/fuzhican/software/conda/envs/circ/bin/bamToFastq -i  ${HISAT_map}/${srr}_HISAT2_unmapped.bam  -fq /dev/stdout | gzip >${HISAT_map}/${srr}_HISAT2_unmapped.fastq.gz
        ### 2. BWA 6 mismatches mapping
        $bwa mem -t ${threads} ${ref_genome_path}  ${HISAT_map}/${srr}_HISAT2_unmapped.fastq.gz > ${bwa_map}/${srr}_bwa_mapped.sam
    else
        echo "Unknown layout:$layout, please specify single or paired." 

    fi
    python ${sh_path}/bwa_unique_mismatch6.py ${bwa_map}/${srr}_bwa_mapped.sam ${bwa_map}/${srr}_bwa_unique_mis6_mapq0.sam 
    samtools_log_file=${bwa_map}/samtools_${srr}.log 
    ${samtools1} view -bT ${ref_genome_path} -o ${bwa_map}/${srr}_unmapped.nbam ${bwa_map}/${srr}_bwa_unique_mis6_mapq0.sam 2>>$samtools_log_file
    ${samtools1} sort ${bwa_map}/${srr}_unmapped.nbam -o ${bwa_map}/${srr}_unmapped.sort.bam 2>>$samtools_log_file
    ${samtools1} view -H ${bwa_map}/${srr}_unmapped.sort.bam > ${bwa_map}/${srr}_bwa.header 2>>$samtools_log_file 
    
    wait
    cat ${bwa_map}/${srr}_bwa.header ${HISAT_map}/${srr}_unique_mismatch2.sam > ${combine_bam}/${srr}_accepted_hits.nsam 2>>$samtools_log_file
    ${samtools1} view -bT ${ref_genome_path} -o  ${combine_bam}/${srr}_accepted_hits.nbam ${combine_bam}/${srr}_accepted_hits.nsam 2>>$samtools_log_file
    ${samtools1} sort ${combine_bam}/${srr}_accepted_hits.nbam -o ${combine_bam}/${srr}_accepted_hits.sort.bam 2>>$samtools_log_file

    ${samtools1} merge -f ${combine_bam}/${srr}_combine.bam ${combine_bam}/${srr}_accepted_hits.sort.bam ${bwa_map}/${srr}_unmapped.sort.bam 2>>$samtools_log_file
    ${samtools1} flagstat -@ $threads ${combine_bam}/${srr}_combine.bam >${combine_bam}/${srr}_combine_flagstat.log 
    ${samtools1} index ${combine_bam}/${srr}_combine.bam 2>>$samtools_log_file
    need_rm_internal_file_list=(${combine_bam}/${srr}_accepted_hits.nsam ${combine_bam}/${srr}_accepted_hits.nbam ${combine_bam}/${srr}_accepted_hits.sort.bam ${bwa_map}/${srr}_bwa.header ${bwa_map}/${srr}_unmapped.nbam ${bwa_map}/${srr}_unmapped.nbam ${bwa_map}/${srr}_unmapped.sort.bam ${bwa_map}/${srr}_bwa_unique_mis6_mapq0.sam ${bwa_map}/${srr}_bwa_mapped.sam ${HISAT_map}/${srr}_HISAT2_mapped.sam ${HISAT_map}/${srr}_unique_mismatch2.sam ${HISAT_map}/${srr}_unmapped_1.fastq.gz ${HISAT_map}/${srr}_unmapped_2.fastq.gz ${HISAT_map}/${srr}_HISAT2_unmapped.fastq.gz  ${HISAT_map}/${srr}_hisat2_unmap.readid ${HISAT_map}/${srr}_un_conc_2.fastq.gz ${HISAT_map}/${srr}_un_conc_1.fastq.gz) #${HISAT_map}/${srr}_HISAT2_unmapped.bam

    echo "#####[`date`]${srr} satrt rm internal files"
    for need_rm_internal_file in ${need_rm_internal_file_list[@]}
    do
        test -e $need_rm_internal_file && rm $need_rm_internal_file
    done
    echo "#####[`date`]${srr} end rm internal files"
}
sam_fine_tune(){
    #sam_fine_tune $bam_file $fine_tune_path $bam_file_basename_prefix
    bam_file=$1
    fine_tune_path=$2
    bam_file_basename_prefix=$3  ###### Tue Oct 8 19:27:03 CST 2019 fzc added
    test -d $fine_tune_path || mkdir -p $fine_tune_path
    if [ "$bam_file_basename_prefix" == "" ];then
        bam_file_basename_prefix=`basename $bam_file|awk -F".bam" '{print $1}'`
    fi
    bam_file_prefix=${fine_tune_path}/${bam_file_basename_prefix}

    picard=/picb/rnomics4/rotation/fuzhican/bin/picard.jar
    gatk4=/picb/rnomics4/rotation/fuzhican/bin/gatk4
    #knownSNP_for_BQSR=${dep_path}/${ref_genome}/NCBI_dbSNP_b151_all_${ref_genome}.vcf
    knownSNP_for_BQSR=${dep_path}/${ref_genome}/NCBI_dbSNP_all_${ref_genome}.vcf ###### Wed Feb 12 22:05:18 CST 2020 fzc start
    test -e ${knownSNP_for_BQSR} ||{
        knownSNP_for_BQSR=${dep_path}/${ref_genome}/NCBI_dbSNP_all_${ref_genome}.vcf.gz
    }
    test -e ${knownSNP_for_BQSR} ||{
        knownSNP_for_BQSR=${dep_path}/${ref_genome}/NCBI_dbSNP_b151_all_${ref_genome}.vcf
    } 
    test -e ${knownSNP_for_BQSR} ||{
         echo "check knownSNP_for_BQSR, ${dep_path}/${ref_genome}/NCBI_dbSNP_all_${ref_genome}.vcf, ${dep_path}/${ref_genome}/NCBI_dbSNP_all_${ref_genome}.vcf.gz and ${dep_path}/${ref_genome}/NCBI_dbSNP_b151_all_${ref_genome}.vcf non-exist!"
         exit
          } ###### Wed Feb 12 22:05:18 CST 2020 fzc add 
    java -jar ${picard} AddOrReplaceReadGroups I=${bam_file} O=${bam_file_prefix}_rgadd.bam SO=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20 

    java -jar $picard MarkDuplicates I=${bam_file_prefix}_rgadd.bam O=${bam_file_prefix}_rgadd_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${bam_file_prefix}_rgadd_MarkDuplicates_output.metrics 

    ${gatk4} SplitNCigarReads -R ${ref_genome_path} -I ${bam_file_prefix}_rgadd_dedupped.bam -O ${bam_file_prefix}_rgadd_dedupped_split.bam

    $gatk4  BaseRecalibrator -R ${ref_genome_path} -I ${bam_file_prefix}_rgadd_dedupped_split.bam --known-sites ${knownSNP_for_BQSR}  -O ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk4.grv 
    $gatk4 ApplyBQSR  -R ${ref_genome_path} -I ${bam_file_prefix}_rgadd_dedupped_split.bam  --bqsr-recal-file ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk4.grv -O ${bam_file_prefix}_rgadd_dedupped_split_recal.bam

    ${samtools1} flagstat -@ $threads ${bam_file_prefix}_rgadd_dedupped_split_recal.bam >${bam_file_prefix}_rgadd_dedupped_split_recal_flagstat.log & 

    tmp_remove_file_list=(${bam_file_prefix}_rgadd.bam ${bam_file_prefix}_rgadd_dedupped.bam ${bam_file_prefix}_rgadd_dedupped.bai ${bam_file_prefix}_rgadd_dedupped_split.bam ${bam_file_prefix}_rgadd_dedupped_split.bai ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk4.grv ${bam_file_prefix}_rgadd_MarkDuplicates_output.metrics)

    for file1 in ${tmp_remove_file_list[@]}
    do
    test -e $file1 && rm $file1
    done
}
gatk_Variant_filter_hg19(){
    ###### Sun Nov 10 16:42:24 CST 2019
    #gatk_Variant_filter $bam_file $tag $m0_path
    #tag="_2pass"
    #m0_path=${ABE_m_path}
    local from_bam=$1
    local tag=$2
    local m0_path=$3
    local srr=$4
    
    m1_path=${m0_path}/${ref_genome}/gatk
    local var_call=${m1_path}/4-0-0var_calling
    test -d ${var_call}||mkdir -p ${var_call}

    #--base-quality-score-threshold default 18
    $gatk4 HaplotypeCaller -R ${ref_genome_path} -I ${from_bam}  --minimum-mapping-quality 0 --dont-use-soft-clipped-bases true  -stand-call-conf 0 -O ${var_call}/${srr}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0.vcf  ###### Sun Nov 10 16:37:06 CST 2019 remove "--min-base-quality-score 20"
    #default --min-base-quality-score 10

    echo -e "[`date`]BEGIN filter..."
    local vcf_filter_path=${m1_path}/5-0-0vcf_filter/${srr}
    test -d ${vcf_filter_path}||mkdir -p ${vcf_filter_path}
    #extract SNP
    $gatk4 SelectVariants -select-type SNP -R ${ref_genome_path} -V ${var_call}/${srr}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0.vcf -O ${vcf_filter_path}/${srr}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP.vcf
    inter_name="_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP"
    
    #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_hg19.sh run_2013_NatMethods_filter_SNP_flexible gatk $srr ${m0_path} ${inter_name}
    #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_hg19.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr ${m0_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam
    bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m gatk -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
    bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m gatk -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" -t "" -b $from_bam -g ${ref_genome}





    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible gatk $srr ${m0_path} ${inter_name}
    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr ${m0_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam
}
gatk_Variant_filter(){
    #gatk_Variant_filter $bam_file $tag $m0_path
    #tag="_2pass"
    #m0_path=${ABE_m_path}
    local from_bam=$1
    local tag=$2
    local m0_path=$3
    local srr=$4
    local ref_genome_tmp=$5
    a_tmp=aa${ref_genome_tmp}aa
    test "$a_tmp" != "aaaa" && ref_genome=$ref_genome_tmp
    
    m1_path=${m0_path}/${ref_genome}/gatk
    local var_call=${m1_path}/4-0-0var_calling
    test -d ${var_call}||mkdir -p ${var_call}

    #--base-quality-score-threshold default 18
    $gatk4 HaplotypeCaller -R ${ref_genome_path} -I ${from_bam} --minimum-mapping-quality 0 --dont-use-soft-clipped-bases true  -stand-call-conf 0 -O ${var_call}/${srr}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0.vcf ###### Sun Nov 10 16:37:06 CST 2019 remove "--min-base-quality-score 20"

    echo -e "[`date`]BEGIN filter..."
    local vcf_filter_path=${m1_path}/5-0-0vcf_filter/${srr}
    test -d ${vcf_filter_path}||mkdir -p ${vcf_filter_path}
    #extract SNP
    $gatk4 SelectVariants -select-type SNP -R ${ref_genome_path} -V ${var_call}/${srr}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0.vcf -O ${vcf_filter_path}/${srr}_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP.vcf
    inter_name="_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP"
    
    bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m gatk -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
    bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m gatk -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" -t "" -b $from_bam -g ${ref_genome}

    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_SNP_flexible gatk $srr ${m0_path} ${inter_name}
    #bash ${sh_path}/support_flow_19_7_8.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr ${m0_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam
}
bmc_Variant_filter(){
    #bmc_Variant_filter $bam_file $tag $m0_path $srr
    #Used for: Calling and filtering, Calling mutation use samtools mpileup, filtering same as gatk_Variant_filter.
    #tag="_recal"
    #m0_path=${ABE_m_path}
    local from_bam=$1
    local tag=$2
    local m0_path=$3
    local srr=$4
    local ref_genome_tmp=$5
    a_tmp=aa${ref_genome_tmp}aa
    test "$a_tmp" != "aaaa" && ref_genome=$ref_genome_tmp
    #####All editing sites (BQ20; overhang6; ES95)
    m1_path=${m0_path}/${ref_genome}/bmc
    local editing_sites=${m1_path}/4-0-0Editing_sites/${srr}
    test -d ${editing_sites} ||mkdir -p ${editing_sites}
    #$samtools1 mpileup -Q 0 -d 10000000  -OI -f ${ref_genome_path}  $from_bam |gzip >${editing_sites}/${srr}_mpileup_direct_out.gz & ###### Sun Feb 23 11:11:58 CST 2020

    perl ${sh_path}/npileupBam_sszhu.pl.backup -i $from_bam -s ${ref_genome_path} -depth 10000000 -minBQ 20 -o 6 -HPB 0 -eSignal 0.95 -v 0 --cRatio 0 >${editing_sites}/${srr}.BQ20o6ES95v0${tag}

    perl -ane 'print "$F[0]:$F[1]\t",join("\t",@F[2..$#F]),"\n"' ${editing_sites}/${srr}.BQ20o6ES95v0${tag} >${editing_sites}/${srr}.BQ20o6ES95v0.new${tag}
    #pigz -p 10 ${editing_sites}/${srr}.BQ20o6ES95v0${tag} &
    test -e  ${editing_sites}/${srr}.BQ20o6ES95v0${tag} && rm  ${editing_sites}/${srr}.BQ20o6ES95v0${tag}


    python ${sh_path}/npileup_to_ES_variant_allvariants.py ${editing_sites}/${srr}.BQ20o6ES95v0.new${tag} 0.95 2 >${editing_sites}/${srr}.BQ20o6ES95v2.allvariants${tag}
    rm ${editing_sites}/${srr}.BQ20o6ES95v0.new${tag} &
    if [ "$ref_genome" == "hg38" ] || [ "$ref_genome" == "hg19" ];then
        inter_name=".BQ20o6ES95v2.allvariants"

        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" -t "" -b $from_bam -g ${ref_genome}
        #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_SangerInstitute" -t "" -b $from_bam -g ${ref_genome}

    elif [ "$ref_genome" == "mm10" ];then
        inter_name=".BQ20o6ES95v2.allvariants"
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_SangerInstitute" -t "" -b $from_bam -g ${ref_genome}
    elif [ "$ref_genome" == "rheMac10" ];then
        inter_name=".BQ20o6ES95v2.allvariants"
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP" -t "" -b $from_bam -g ${ref_genome}
    elif [ "$ref_genome" == "ce11" ];then
        inter_name=".BQ20o6ES95v2.allvariants"
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP" -t "" -b $from_bam -g ${ref_genome}
    else
        echo "unsupport genome_build_version: $ref_genome"
        exit 1
    fi
    wait

}
patch_begin_at_Blat(){
    #patch_begin_at_Blat $tmp_work_path_0 $srr
    ###### Thu Sep 26 17:25:11 CST 2019
    #useed for continue gatk_Variant_filter at Blat point
    #
    tmp_work_path_0=$1
    srr=$2
    #20190906_1_HaplotypeCaller_Variants_2pass_BQ20_MAPQ0_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_v2_non-Alu_v3_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand_Nonrepetitive.vcf
    inter_name="_HaplotypeCaller_Variants_2pass_BQ20_MAPQ0_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_v2_non-Alu_v3_deSimpleRepeat_intronic4bp_deHomopolymer"
    method="gatk"
    m0_path=$tmp_work_path_0
    m1_path=${m0_path}/${ref_genome}/gatk
    vcf_filter_path=${m1_path}/5-0-0vcf_filter/${srr}
    out_final_part_name=${inter_name}"_blat" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "prepare_for_blat_2013_Natmethods_version $input_result $bam_file $output_result $method"
    bash ${sh_path}/bmc_or_gatk_filter_spurious_sites.sh prepare_for_blat_2013_Natmethods_version $input_result $bam_file $output_result $method

    ################################strandness################################

    inter_name=$out_final_part_name 
    out_final_part_name=${inter_name}"_annotationStrand" 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    output_result="${tmp_work_path}/${srr}${out_final_part_name}${suffix}"
    echo "distinguish_plus_minus_main_annotation_version $method $input_result $output_result"
    bash ${sh_path}/bmc_or_gatk_filter_spurious_sites.sh distinguish_plus_minus_main_annotation_version $method $input_result $output_result

    ################################split_non-alu_to_nonAlu_nonRepeat#################################

    #_HaplotypeCaller_Variants_2pass_SNP_deAllSNP_dbSNP_b151_1000genomes_EVS_hits10_non-Alu_deSimpleRepeat_intronic4bp_deHomopolymer_blat_annotationStrand.
    inter_name=$out_final_part_name 
    input_result="${tmp_work_path}/${srr}${inter_name}${suffix}"
    #echo $input_result
    echo "core_of_split_non_alu_to_nonAlu_nonRepeat $method $input_result #output name is auto determined by input "
    bash ${sh_path}/bmc_or_gatk_filter_spurious_sites.sh core_of_split_non_alu_to_nonAlu_nonRepeat $method $input_result #output name is auto determined by input 

}
patch_begin_at_filter(){
    #patch_begin_at_filter $tmp_work_path_0 $srr
    ###### Thu Sep 26 18:23:58 CST 2019
    ###### Wed Nov 20 20:09:08 CST 2019
    #patch_begin_at_filter -o ${output_dir} -n TET1_ctrl_1st_HEK293FT -g ${ref_genome}
    #output_dir=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/results_after_19_9_25/
    #bash /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/sh/GATK_RNA_seq_HISAT2_BWA_20_1_15.sh patch_begin_at_filter -o ${output_dir} -n NT-BE3_rep1
    #output_dir=""
    #output_dir=/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/tmp_dir
    #name=TET1_ctrl_1st_HEK293FT
    #bash /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/sh/GATK_RNA_seq_HISAT2_BWA_20_1_15.sh patch_begin_at_filter -o ${output_dir} -n ${name}

    
    while getopts :o:n:g: ARGS  
    do  
    case $ARGS in   
        o)  
            m0_path=$OPTARG
            ;;  
        n)  
            name=$OPTARG
            ;; 
        g)  
            ref_genome=$OPTARG
            ;; 
        *)  
            echo "Unknown option: $ARGS"
            ;;
        \?)
        echo "Invalid option: -$OPTARG" 
        ;;
    esac
    done
    tmp_work_path_0=$m0_path
    srr=$name
    tag="_2pass"
    bam_file_prefix=${srr}_combine
    fine_tune_path="${tmp_work_path_0}/${ref_genome}/3-0-0Combine_bam"
    #bam_file="${tmp_work_path_0}/${ref_genome}/bmc/3-0-0Combine_bam/${bam_file_prefix}_readgroup_sort_dedupped_recal_gatk4.bam"
    #bam_file=${tmp_work_path_0}/THP1_0hr_pAminus_tmp/3-0-0Combine_bam/THP1_0hr_pAminus_combine_rgadd_dedupped_split_recal.bam
    method=bmc
    #inter_name="_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP"
    #inter_name="_HaplotypeCaller_Variants${tag}_BQdefault_MAPQ0_SNP"
    #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr ${m0_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam

    #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_hg19.sh run_2013_NatMethods_filter_SNP_flexible gatk $srr ${m0_path} ${inter_name}
    
    #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_hg19.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible gatk $srr ${m0_path} "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" "" $from_bam
    
    #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_hg19.sh run_2013_NatMethods_filter_SNP_flexible -m gatk -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
    #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m gatk -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" -t "" -b $from_bam -g ${ref_genome}

    if [ "$method" == "bmc" ];then
    {
        bam_file="${tmp_work_path_0}/${ref_genome}/3-0-0Combine_bam/${bam_file_prefix}_rgadd_dedupped_split_recal.bam"
        from_bam=$bam_file
        if [ "$ref_genome" == "hg38" ] || [ "$ref_genome" == "hg19" ];then
            inter_name=".BQ20o6ES95v2.allvariants"
            #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
            bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" -t "" -b $from_bam -g ${ref_genome}
            #bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}" -t "" -b $from_bam -g ${ref_genome}
        elif [ "$ref_genome" == "mm10" ];then
            inter_name=".BQ20o6ES95v2.allvariants"
            bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
            bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_SangerInstitute" -t "" -b $from_bam -g ${ref_genome}
        elif [ "$ref_genome" == "rheMac10" ];then
            inter_name=".BQ20o6ES95v2.allvariants"
            bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
            bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP" -t "" -b $from_bam -g ${ref_genome}
        elif [ "$ref_genome" == "ce11" ];then
            inter_name=".BQ20o6ES95v2.allvariants"
            bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_SNP_flexible -m bmc -n $srr -o ${m0_path} -i ${inter_name} -g ${ref_genome}
            bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m bmc -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP" -t "" -b $from_bam -g ${ref_genome}
        fi
    }
    elif [ "$method" == "gatk" ];then
    {
        inter_name="_HaplotypeCaller_Variants_2pass_BQdefault_MAPQ0_SNP"
        bash ${sh_path}/bmc_or_gatk_filter_spurious_sites_v2.sh run_2013_NatMethods_filter_regions_in_bed_core_flexible -m gatk -n $srr -o ${m0_path} -i "${inter_name}_deAllSNP_dbSNP_b151_1000genomes_EVS" -t "" -b $from_bam -g ${ref_genome}
    }
    else
    echo "unsupport method:$method"
    fi

}
run_strandness_Alu(){
    ###### Mon May 25 17:30:17 CST 2020
    m0_path=$1
    srr=$2
    method=$3
    ref_genome=$4
    m_path=$m0_path
    bash /data/rnomics6/fuzhican/project/ABE_transcriptomics_off_target/sh/bmc_or_gatk_filter_spurious_sites_v2.sh strandness_Alu $m0_path $srr $method $ref_genome
}
#bash /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/sh/GATK_RNA_seq_HISAT2_BWA_19_9_25.sh RNA_Editing_Calling_Pipeline_HISAT2_BWA_followed_by_GATK_HaplotypeCaller -1 /picb/rnomics4/xuew/cell/TET1_ctrl_1st_HEK293FT_R1_trimmed.fq -2 /picb/rnomics4/xuew/cell/TET1_ctrl_1st_HEK293FT_R2_trimmed.fq  -o /picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target/tmp_dir -n TET1_ctrl_1st_HEK293FT -g hg19 -t 16

"$@"
