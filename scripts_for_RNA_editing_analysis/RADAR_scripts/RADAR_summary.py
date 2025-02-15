import pandas as pd
import numpy as np
import time,os 
class SNP_editing:
    def __init__(self):
        pass
    def check_dep_files(self,dep_file,function_name,site_list_file,ref_genome=None):
        #annovar_annotation_summary_table_para
        if ref_genome==None:ref_genome=self.ref_genome
        if not ( os.path.exists(dep_file) and os.path.getsize(dep_file)>100):
            print(("check_dep_files",dep_file,function_name))
            cmd="bash %s/RADAR_summary_depfiles_prepare.sh %s %s %s"%(sh_path,function_name,ref_genome,site_list_file)
            print(cmd)
            os.system(cmd)
    def run_editing_ratio_files_summary_total_samples_paraList(self):
        site_range="filtered_SNP"
        site_range="Four_validate_sites" 
        site_range="all" 
        site_range="filtered_spurious" 
        # self.site_range=site_range
        self.wo_filtered_ER="_filtered_ER"
        self.wo_filtered_ER=""
        #self.whether_semi_samples="_semi-samples"
        self.whether_semi_samples=""
        self.at_least_occur_3_times="_occur3" ###### Mon Dec 30 16:36:34 CST 2019
        self.at_least_occur_3_times=""
        # self.whether_distnguish_NA_0="_Distinguish_NA0"
        self.whether_distnguish_NA_0="" ###### Wed Aug 11 20:15:32 CST 2021
        #if site_range== "all":
        self.whether_only_AGTC_CTGA="_only_AGTC_CTGA"
        self.whether_only_AGTC_CTGA=""
        self.whether_only_AGTC="_only_AGTC"
        self.whether_only_AGTC=""
    def run_editing_ratio_files_summary_total_samples_big_file(self,SampleList_SepBySemicolon=None,sample_source_token=None,results_path=None,summary_path=None,sample_name_parser="%;0;0",species="Human",ref_genome="hg38",sample_source="2019_nature",outfile=None):  
        print("species=\"%s\",ref_genome=\"%s\",sample_source=\"%s\""%(species,ref_genome,sample_source))
        print("sample_name_parser:"+str(sample_name_parser))
        #################fine-tune_para#############
        self.ref_genome=ref_genome
        self.sample_source=sample_source
        #################fine-tune_para#############
        self.method="bmc"  
        # SampleSource=species+"_all_samples"
        t1_GeneBody={'exonic;splicing': 'splicing', 'upstream;downstream': 'intergenic', 'downstream': 'intergenic', 'splicing': 'splicing', 'UTR5': 'UTR5', 'ncRNA_exonic': 'exonic', 'intergenic': 'intergenic', 'intronic': 'intronic', 'exonic': 'exonic', 'ncRNA_exonic;splicing': 'splicing', 'UTR3': 'UTR3', 'ncRNA_splicing': 'splicing', 'upstream': 'intergenic', 'UTR5;UTR3': 'UTR3', 'ncRNA_intronic': 'intronic'}
        self.df_Genome_Location=""
        self.df_HostGene=""
        self.df_Extract_context_seq=""
        self.df_Exon_Intron_interval=""
        inter_name=".BQ20o6ES95v2"
        suffixName_part_list_str="%(inter_name)s.allvariants_recal %(inter_name)s.allvariants_annotationStrand_recal _recommended_final_result_All_recal %(inter_name)s.allvariants_deAllSNP_dbSNP_b151_1000genomes_EVS_recal"%{"inter_name":inter_name} ######
        suffixName_part_list_key_str="all all_site_annotationStrand filtered_spurious filtered_SNP"
        sample_list=SampleList_SepBySemicolon.split(";")
        ABE_m_path="/picb/rnomics4/rotation/fuzhican/project/ABE_transcriptomics_off_target"
        if ref_genome=="hg38":
            ###dep files
            RADAR_db_file=os.path.join(ABE_m_path,"dep_files/hg38/RADAR_All_Sites_liftOver_hg38_sort_20200508.sites")
            REDIportal_db_file=os.path.join(ABE_m_path,"dep_files/hg38/REDIportal_All_Sites_liftOver_hg38_sort_20200719.sites")
            HostGene_info_file=os.path.join(summary_path,"Editing_ratio_summary/ANNOVAR_splicing_threshold_3/%s_ANNOVAR_annatation_hg38Refall3.variant_function"%sample_source_token)
            info_filename=os.path.join(summary_path,"Editing_ratio_summary/"+sample_source_token)
            site_list_file=info_filename+".tsv"
            Genome_Location_file=info_filename+".Genome_Location"
            Extract_context_seq_file=info_filename+".Extract_context_seq"
            Exon_Intron_interval_file=info_filename+".Exon_Intron_interval"
            filtered_common_SNP_file=info_filename+".Annotate_common_SNP"
            print(HostGene_info_file,info_filename)
            ###RADAR known sites
            df_RADAR_db_file=pd.read_csv(RADAR_db_file,sep="\t",header=None)
            df_RADAR_db_file.rename(columns={0:"Site"},inplace=True)
            df_RADAR_db_file.set_index("Site",inplace=True)
            self.df_RADAR_db_file=df_RADAR_db_file.copy()
            ###REDIportal known sites
            df_REDIportal_db_file=pd.read_csv(REDIportal_db_file,sep="\t",header=None)
            df_REDIportal_db_file.rename(columns={0:"Site"},inplace=True)
            df_REDIportal_db_file.set_index("Site",inplace=True)
            self.df_REDIportal_db_file=df_REDIportal_db_file.copy()
        elif ref_genome=="rheMac10" or ref_genome=="mm10" :
            ###dep files
            HostGene_info_file=os.path.join(summary_path,"Editing_ratio_summary/ANNOVAR_splicing_threshold_3/%s_ANNOVAR_annatation_%sRefall3.variant_function"%(sample_source_token,ref_genome))
            info_filename=os.path.join(summary_path,"Editing_ratio_summary/"+sample_source_token)
            site_list_file=info_filename+".tsv"
            Genome_Location_file=info_filename+".Genome_Location"
            Extract_context_seq_file=info_filename+".Extract_context_seq"
            Exon_Intron_interval_file=info_filename+".Exon_Intron_interval"
            if ref_genome=="mm10":
                REDIportal_db_file=os.path.join(ABE_m_path,\
                                                "dep_files/mm10/REDIportal_All_Sites_mm10_20200720.sites")
                ###REDIportal known sites
                df_REDIportal_db_file=pd.read_csv(REDIportal_db_file,sep="\t",header=None)
                df_REDIportal_db_file.rename(columns={0:"Site"},inplace=True)
                df_REDIportal_db_file.set_index("Site",inplace=True)
                self.df_REDIportal_db_file=df_REDIportal_db_file.copy()
        elif ref_genome=="galGal6" or ref_genome=="Bter1" or ref_genome=="ce11":
            ###dep files
            HostGene_info_file=os.path.join(summary_path,"Editing_ratio_summary/ANNOVAR_splicing_threshold_3/%s_ANNOVAR_annatation_%sRefall1.variant_function"%(sample_source_token,ref_genome))
            info_filename=os.path.join(summary_path,"Editing_ratio_summary/"+sample_source_token)
            site_list_file=info_filename+".tsv"
            Genome_Location_file=info_filename+".Genome_Location"
            Extract_context_seq_file=info_filename+".Extract_context_seq"
            Exon_Intron_interval_file=info_filename+".Exon_Intron_interval"
        else:
            print("ref_genome(%s) is not supported!"%ref_genome)
        t1_function_dep_file={dep_file:function_name for dep_file,function_name in zip("HostGene_info_file Genome_Location_file Extract_context_seq_file Exon_Intron_interval_file filtered_common_SNP_file".split() ,"annovar_annotation_summary_table_Para Annotate_Genome_Location_all_tissue_Para extract_context_seq_Para Annotate_Exon_Intron_interval_all_tissue_Para Annotate_Annotate_common_SNP_all_tissue_Para".split())}
        for dep_file, function_name in t1_function_dep_file.items():
            if dep_file in locals(): 
                self.check_dep_files(locals()[dep_file],function_name,site_list_file)
        if ref_genome=="hg38":
            ###common SNP
            self.df_filtered_common_SNP=pd.read_csv(filtered_common_SNP_file,sep="\t",header=None)
            self.df_filtered_common_SNP.rename(columns={0:"Site",1:"filtered_common_SNP"},inplace=True)
            self.df_filtered_common_SNP.set_index("Site",inplace=True)

        ###HostGene 
        print(HostGene_info_file)
        df_HostGene=pd.read_csv(HostGene_info_file,sep="\t",header=None)
        df_HostGene=df_HostGene.assign(HostGene_select=df_HostGene.loc[:,1].str.split(r"\(|,|;",n=1,expand=True).loc[:,0],
                        HostGene_count=df_HostGene.loc[:,1].str.replace('\(.+?\)', '', regex=True).str.split(",").apply(lambda x:len(x)),
                        GeneBody_select=df_HostGene.iloc[:,0].map(t1_GeneBody),
                        Site=lambda x:x[2]+":"+x[3].astype(str)
                        )
        df_HostGene=df_HostGene[["Site",1,"HostGene_select","HostGene_count",0,"GeneBody_select"]]
        df_HostGene.rename(columns={1:"HostGene_ori",0:"GeneBody"},inplace=True)
        df_HostGene.set_index("Site",inplace=True)
        self.df_HostGene=df_HostGene.copy()

        ###Genome_Location
        if ref_genome=="hg38" or ref_genome=="rheMac10" or ref_genome=="mm10" :
            df_Genome_Location=pd.read_csv(Genome_Location_file,sep="\t",header=None)
            df_Genome_Location=df_Genome_Location[[0,2]]
            df_Genome_Location.rename(columns={0:"Site",2:"Genome_Location"},inplace=True)
            df_Genome_Location.set_index("Site",inplace=True)
            self.df_Genome_Location=df_Genome_Location.copy()
        ###RADAR known sites
        self.df_Extract_context_seq=pd.read_csv(Extract_context_seq_file,sep="\t",header=None,names=["Site","Context_20bp_seq"],index_col=0)
        self.df_Exon_Intron_interval=pd.read_csv(Exon_Intron_interval_file,sep="\t",header=None,names=["Site","Exon_Intron_interval"],usecols=(0,2),index_col=0)
        self.df_Exon_Intron_interval=self.df_Exon_Intron_interval.assign(Exon_Intron_interval_select=self.df_Exon_Intron_interval.iloc[:,0].str.split(r",",n=1,expand=True).loc[:,0])

        if ref_genome=="hg38":
            suffixName_part_list_str="%(inter_name)s.allvariants_recal %(inter_name)s.allvariants_annotationStrand_recal _recommended_final_result_All_recal %(inter_name)s.allvariants_deAllSNP_dbSNP_b151_1000genomes_EVS_recal"%{"inter_name":inter_name} ######
            suffixName_part_list_key_str="all all_site_annotationStrand filtered_spurious filtered_SNP"
            self.suffixName_part_t={suffixName_part_list_key_str.split()[x]:suffixName_part_list_str.split()[x] for x in range(0,len(suffixName_part_list_str.split()))}
        elif ref_genome=="rheMac10":
            suffixName_part_list_str="%(inter_name)s.allvariants_recal _recommended_final_result_All_recal %(inter_name)s.allvariants_deAllSNP_dbSNP_recal"%{"inter_name":inter_name} ######
            suffixName_part_list_key_str="all filtered_spurious filtered_SNP"
            self.suffixName_part_t={suffixName_part_list_key_str.split()[x]:suffixName_part_list_str.split()[x] for x in range(0,len(suffixName_part_list_str.split()))} 
        elif ref_genome=="mm10":
            suffixName_part_list_str="%(inter_name)s.allvariants_recal _recommended_final_result_All_recal %(inter_name)s.allvariants_deAllSNP_dbSNP_SangerInstitute_recal"%{"inter_name":inter_name} ######
            suffixName_part_list_key_str="all filtered_spurious filtered_SNP"
            self.suffixName_part_t={suffixName_part_list_key_str.split()[x]:suffixName_part_list_str.split()[x] for x in range(0,len(suffixName_part_list_str.split()))} 
        elif ref_genome=="Bter1" or ref_genome=="galGal6" or ref_genome=="ce11":
            suffixName_part_list_str="%(inter_name)s.allvariants_recal"%{"inter_name":inter_name} ######
            suffixName_part_list_key_str="all"
            self.suffixName_part_t={suffixName_part_list_key_str.split()[x]:suffixName_part_list_str.split()[x] for x in range(0,len(suffixName_part_list_str.split()))}
        else:
            print("Please check ref_genome(%s), whose value is not supported!"%ref_genome)
            exit()
        
        
        self.results_path_prefix=results_path
        self.list_df_sample=[self.editing_ratio_files_summary_total_samples_big_file_one_sample(sample_name,sample_name_parser) for sample_name in sample_list]
        df_source_all_results=pd.concat(self.list_df_sample,sort=False)
        if outfile==None:
            date_token=time.strftime("%Y_%m_%d", time.localtime())  
            outfile=os.path.join(summary_path,species+"_"+date_token+"_big_file_sample_source-%s_site_range-%s.tsv.gz"%(sample_source,self.site_range))

        df_source_all_results.to_csv(outfile,sep="\t")
        print(outfile)
       
        
    def editing_ratio_files_summary_total_samples_big_file_one_sample(self,sample_name,sample_name_parser=None):
        print("editing_ratio_files_summary_total_samples_big_file_one_sample,sample_name,sample_name_parser",sample_name,sample_name_parser)
        using_site_range="all"
        if self.method!="bmc":
            print("method(%s) is not support!"%(self.method))
            exit()
        result_path=os.path.join(self.results_path_prefix,sample_name)
        if sample_name_parser!=None:
            sample_name_sep,SampleSource_id,Treat_id=sample_name_parser.split(";")
            print("sample_name,sample_name_sep,SampleSource_id,Treat_id:",sample_name,sample_name_sep,SampleSource_id,Treat_id)
            Treat=sample_name.split(sample_name_sep)[int(Treat_id)]
            SampleSource=sample_name.split(sample_name_sep)[int(SampleSource_id)]
        # print(Treat,SampleSource)
        source_all_results_file=os.path.join(result_path,sample_name+self.suffixName_part_t[using_site_range])
        if not os.path.exists(source_all_results_file):
            return(None)
        #print(source_all_results_file)
        df_source_all_results=pd.read_csv(source_all_results_file,sep="\t",header=None)
        
        df_source_all_results[['mut','ER']] = df_source_all_results.iloc[:,16].str.split(" ",expand=True) 
        # slow solution:
        #df_source_all_results[['mut','ER']] = df_source_all_results.iloc[:,16].apply(lambda x: pd.Series(str(x).split(" "))) 
        df_source_all_results["BS"]=df_source_all_results.loc[:,1]+">"+df_source_all_results.loc[:,"mut"]
        df_source_all_results=df_source_all_results[[0,"BS","ER",15,2]]
        df_source_all_results.rename(columns={0:"Site",15:"HPB",2:"Covered_reads"},inplace=True)
        df_source_all_results=df_source_all_results.assign(Sample_name=sample_name,Treat=Treat,
                                        SampleSource=SampleSource)
        df_source_all_results.set_index("Site",inplace=True)
        for processing_file_suffix_key,processing_file_suffix in self.suffixName_part_t.items():
            if processing_file_suffix_key==using_site_range:continue
            processing_file=os.path.join(result_path,sample_name+processing_file_suffix)
            if processing_file_suffix_key=="all_site_annotationStrand":
                df_one_sample_ann=pd.read_csv(processing_file,sep="\s+",header=None,usecols=(0,1,16))
                df_one_sample_ann.rename(columns={0:"Site"},inplace=True)
                df_one_sample_ann["BS_annotated"]=df_one_sample_ann[1]+">"+df_one_sample_ann[16]
                df_one_sample_ann.set_index("Site",inplace=True)
                df_source_all_results=pd.merge(df_source_all_results,df_one_sample_ann["BS_annotated"],how="left",on="Site",copy=False)
                print(processing_file_suffix_key,df_source_all_results.shape,df_one_sample_ann.shape)
                continue
            #print(processing_file)
            df_processing=pd.read_csv(processing_file,sep="\s+",usecols=[0],header=None,names=["Site"],index_col="Site")
            df_source_all_results[processing_file_suffix_key]=np.where(df_source_all_results.index.isin(df_processing.index),1,0)
        if self.ref_genome=="hg38" or self.ref_genome=="mm10":
            df_source_all_results["In_REDIportal"]=np.where(df_source_all_results.index.isin(self.df_REDIportal_db_file.index),1,0)
            if self.ref_genome == "hg38" :
                df_source_all_results["In_RADAR"]=np.where(df_source_all_results.index.isin(self.df_RADAR_db_file.index),1,0)
        df_source_all_results=pd.merge(df_source_all_results, self.df_Exon_Intron_interval, how='left', on="Site",
                    left_index=False, right_index=False, sort=False,
                    suffixes=("", ''), copy=False)
        if df_source_all_results.index.name != "Site":
            df_source_all_results.set_index("Site",inplace=True)
        if self.ref_genome=="hg38" :
            df_source_all_results=pd.concat([df_source_all_results,\
                                            self.df_filtered_common_SNP,self.df_Genome_Location,self.df_HostGene,self.df_Extract_context_seq],join='inner',axis=1,sort=False)   
        elif  self.ref_genome=="rheMac10"  or  self.ref_genome=="mm10" :
            df_source_all_results=pd.concat([df_source_all_results,\
                                            self.df_Genome_Location,self.df_HostGene,self.df_Extract_context_seq],join='inner',axis=1,sort=False) 
        elif  self.ref_genome=="Bter1" or  self.ref_genome=="galGal6"  or  self.ref_genome=="ce11" :
            df_source_all_results=pd.concat([df_source_all_results,\
                                            self.df_HostGene,self.df_Extract_context_seq],join='inner',axis=1,sort=False) 
        print(sample_name,df_source_all_results.shape)
        return(df_source_all_results)
    def run_editing_ratio_files_summary_total_samples(self,species="Human",ref_genome="hg38",sample_source="2019_nature",SampleList_SepBySemicolon=None,results_path=None,output_file=None):  
        print("species=\"%s\",ref_genome=\"%s\""%(species,ref_genome))
       
        #################fine-tune_para#############
        inter_name=".BQ20o6ES95v2"  
        site_range="all"
        site_range=self.site_range  
        print(("site_range: "+site_range+"\nwo_filtered_ER: "+self.wo_filtered_ER+"\nwhether_semi_samples: "+self.whether_semi_samples+"\nwhether_distnguish_NA_0: "+self.whether_distnguish_NA_0+"\nat_least_occur_3_times: "+self.at_least_occur_3_times)+"\nwhether_only_AGTC_CTGA: "+self.whether_only_AGTC_CTGA+"\nwhether_only_AGTC: "+self.whether_only_AGTC)
        #################fine-tune_para#############
        method="bmc"  
        suffixName_part="%s.allvariants_recal"%(inter_name) 
        
        if not os.path.exists(os.path.split(output_file)[0]):
            os.makedirs(os.path.split(output_file)[0])
        self.editing_ratio_files_summary_df_total_samples(method,sample_source,results_path,suffixName_part,SampleList_SepBySemicolon,output_file)  
        print(output_file)
    def _editing_ratio_files_summary_sub_bmc(self,line):  
        site=line.split()[0]  
        er=line.split()[-1]  
        return site,er  
    def editing_ratio_files_summary_df_total_samples(self,method,tissue,results_path,suffixName_part_str,SampleList_SepBySemicolon,output_file):  
        suffixName_part_list=suffixName_part_str.split()  
        time_flow_list=SampleList_SepBySemicolon.split(';')  
 
        t1={}#sites-time-er  
        for srr in time_flow_list:  
            for suffixName_part in suffixName_part_list:  
                file1_full_pathname=os.path.join(results_path,srr,srr+suffixName_part)  
                if not os.path.exists(file1_full_pathname):  
                    print("Not exits!",file1_full_pathname)  
                    continue  
                with open(file1_full_pathname) as file1 :  
                    for i in file1:
                        if method=="bmc":  
                            site,er=self._editing_ratio_files_summary_sub_bmc(i)  
                        t1.setdefault(site,{})  
                        t1[site][srr]=er  
        tmp_t1_file=output_file+"_t1"
        print(tissue+" dict t1 finished!\nwriting to "+tmp_t1_file)
        open(tmp_t1_file,'w').write(str(t1))
        os.system("pigz -p 5 -f "+tmp_t1_file)
        df = pd.DataFrame(t1,index=time_flow_list).T
        df_isin=df.isin([1,np.nan])
        row_name=[x for i,x in enumerate(df.index) if set(df_isin.iloc[i,]) == {True}]
        df=df.drop(row_name,axis=0)
        if self.whether_distnguish_NA_0=="_Distinguish_NA0": ###### Mon Aug 23 20:49:09 CST 2021, fuzhican
            df=df.fillna(value=self.df_t1_chrPlussite_sample_replaced)
        df.to_csv(output_file,sep="\t")
    def run_editing_ratio_files_summary_total_samples_shell(self,BatchName,RADAR_OutPath,SampleName_List,sample_name_parser):
        site_range="all"
        species="Human";ref_genome="hg38"; 
        results_path=os.path.join(RADAR_OutPath,"%s/bmc/4-0-0Editing_sites"%ref_genome)  
        summary_path=os.path.join(RADAR_OutPath,ref_genome+"/bmc/6-0-0Summary")
        sample_source=BatchName
        print(SampleName_List,":\n",SampleName_List)
        date_token=time.strftime("%Y_%m_%d", time.localtime())  
        sample_source_token="%s_%s_site_range-%s"%(sample_source,date_token,site_range)
        self.run_editing_ratio_files_summary_total_samples_paraList()
        self.site_range=site_range
        editing_ratio_files_summary_file=os.path.join(summary_path,"Editing_ratio_summary",sample_source_token+".tsv") 
        editing_ratio_files_summary_total_samples_big_file=os.path.join(summary_path,species+"_"+date_token+"_big_file_sample_source-%s_site_range-%s.tsv.gz"%(sample_source,self.site_range))
        if not os.path.exists(editing_ratio_files_summary_file):
            self.run_editing_ratio_files_summary_total_samples(species=species,ref_genome=ref_genome,sample_source=sample_source,SampleList_SepBySemicolon=SampleName_List,results_path=results_path,output_file=editing_ratio_files_summary_file)
        if not os.path.exists(editing_ratio_files_summary_total_samples_big_file):
            sample_name_parser="%;0;0"
            self.run_editing_ratio_files_summary_total_samples_big_file(SampleList_SepBySemicolon=SampleName_List,sample_source_token=sample_source_token,results_path=results_path,summary_path=summary_path,sample_name_parser=sample_name_parser,species=species,ref_genome=ref_genome,sample_source=sample_source,outfile=editing_ratio_files_summary_total_samples_big_file)
if __name__=="__main__":
    sh_path="/data/rnomics10/yuanguohua/20220925_RNAtBE/RADAR_scripts"
    Batch_list = ["NT","EGFP","ADAR2-P3","ADAR2-A3D-P3","ADAR2-A3F-P3","ADAR2-CT2","ADAR2-A3D-CT2","ADAR2-A3F-CT2"]
    for BatchName in Batch_list:
        RADAR_OutPath="/data/rnomics10/yuanguohua/20220925_RNAtBE/04_RADAR_output/" + BatchName + "_output"
        SampleList_SepBySemicolon=BatchName 
        sample_name_parser = "%;0;0"
        MY_SNP_editing=SNP_editing()
        MY_SNP_editing.run_editing_ratio_files_summary_total_samples_shell(BatchName,RADAR_OutPath,SampleList_SepBySemicolon,sample_name_parser)
