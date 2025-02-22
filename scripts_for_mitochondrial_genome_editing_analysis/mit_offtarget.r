library(plyr)
library(stringr)

#----------- read CFBI output files
path<-"./CFBI_output_files"
full_path<-list.files(path, pattern = "*.xls", full.names = T)
point_result<-data.frame()
for(i in full_path)
{
  #i=full_path[1]
  sam<-str_extract(basename(i),"\\d+")
  res1<-read.table(i,as.is=T,header=F,fill = T)
  res1[,18]<-sam
  res1[,19]<-"20241204"
  point_result<-rbind(point_result,res1,stringsAsFactors=F,make.row.names=F)
}
point_result[,20]<-paste(point_result[,19],point_result[,18],point_result[,1],sep = "_")
colnames(point_result)<-c("gene","pos","ref","ref_3","all_num","mut_num","mut_freq","A_num","A_freq","C_num","C_freq","G_num","G_freq","N_num","N_freq","T_num","T_freq","sample","date","id")
point_result$id2<-paste(point_result$sample,point_result$pos,sep = ":")



#---------sample target editing information
edit_sam<-read.table("./sample_edit_site.txt",as.is=T,header=T,fill = T)
on_site<-c()
on_site_all<-c()
for(i in 1:nrow(edit_sam))
{
  pos<-unlist(strsplit(edit_sam[i,2],"/",fix=T))
  on_site<-c(on_site,paste(edit_sam[i,1],pos,sep = ":"))
  rg<-as.numeric(unlist(strsplit(edit_sam[i,3],"-",fix=T)))
  range<-rg[1]:rg[2]
  on_site_all<-c(on_site_all,paste(edit_sam[i,1],range,sep = ":"))
}



#---------293FT SNV identification
point_result_293<-point_result
blank_sample_293<-c(4,5,6) # sample 4,5,6 are NT i
point_result_293<-point_result_293[point_result_293$sample%in%blank_sample_293,]
snv_293<-point_result_293[point_result_293$mut_freq>=0.1,]
snv_site_293_rep1<-snv_293[snv_293$sample==blank_sample_293[1],"pos"]
snv_site_293_rep2<-snv_293[snv_293$sample==blank_sample_293[2],"pos"]
snv_site_293_rep3<-snv_293[snv_293$sample==blank_sample_293[3],"pos"]
snv_293_site<-union(union(snv_site_293_rep1,snv_site_293_rep2),snv_site_293_rep3)



#-----------a>g off-target evaluation
ave_off_point_result<-point_result[!point_result$pos%in%snv_293_site,]
ave_off_point_result<-ave_off_point_result[!ave_off_point_result$id2%in%on_site_all,]
ave_off_point_result<-ave_off_point_result[ave_off_point_result$ref%in%c("A","T"),]
ave_off_point_result$ag_tc_num<-0
ave_off_point_result$mut_ag_tc_num<-0

ave_off_point_result[ave_off_point_result$ref=="A","ag_tc_num"]<-ave_off_point_result[ave_off_point_result$ref=="A","A_num"]
ave_off_point_result[ave_off_point_result$ref=="T","ag_tc_num"]<-ave_off_point_result[ave_off_point_result$ref=="T","T_num"]
ave_off_point_result[ave_off_point_result$ref=="A","mut_ag_tc_num"]<-ave_off_point_result[ave_off_point_result$ref=="A","G_num"]
ave_off_point_result[ave_off_point_result$ref=="T","mut_ag_tc_num"]<-ave_off_point_result[ave_off_point_result$ref=="T","C_num"]

ave_off_point_result$mut_freq_ag_tc<-ave_off_point_result$mut_ag_tc_num/ave_off_point_result$all_num
ave_off_point_result<-ave_off_point_result[ave_off_point_result$mut_freq_ag_tc>=0.01,]
ave_off_stat<-ddply(ave_off_point_result,.(sample),summarise,mut_reads=sum(mut_ag_tc_num))
ave_off_stat<-ave_off_stat[order(as.numeric(ave_off_stat[,1])),]
point_result_at<-point_result[point_result$ref%in%c("A","T"),]
point_result_at_stat<-ddply(point_result_at,.(sample),summarise,all_at_reads=sum(all_num))
point_result_at_stat<-point_result_at_stat[order(as.numeric(point_result_at_stat$sample)),]

ave_off_stat$all_at_reads<-point_result_at_stat$all_at_reads
ave_off_stat$ave_off_all<-ave_off_stat$mut_reads/ave_off_stat$all_at_reads

write.table(ave_off_stat,"./ave_off_stat_ag.txt",quote = F,sep = "\t",col.names = T,row.names = F)



#-----------c>t off_target
ave_off_point_result_ct<-point_result[!point_result$pos%in%snv_293_site,]
ave_off_point_result_ct<-ave_off_point_result_ct[!ave_off_point_result_ct$id2%in%on_site_all,]
ave_off_point_result_ct<-ave_off_point_result_ct[ave_off_point_result_ct$ref%in%c("C","G"),]
ave_off_point_result_ct$ct_ga_num<-0
ave_off_point_result_ct$mut_ct_ga_num<-0

ave_off_point_result_ct[ave_off_point_result_ct$ref=="C","ct_ga_num"]<-ave_off_point_result_ct[ave_off_point_result_ct$ref=="C","C_num"]
ave_off_point_result_ct[ave_off_point_result_ct$ref=="G","ct_ga_num"]<-ave_off_point_result_ct[ave_off_point_result_ct$ref=="G","G_num"]
ave_off_point_result_ct[ave_off_point_result_ct$ref=="C","mut_ct_ga_num"]<-ave_off_point_result_ct[ave_off_point_result_ct$ref=="C","T_num"]
ave_off_point_result_ct[ave_off_point_result_ct$ref=="G","mut_ct_ga_num"]<-ave_off_point_result_ct[ave_off_point_result_ct$ref=="G","A_num"]

ave_off_point_result_ct$mut_freq_ct_ga<-ave_off_point_result_ct$mut_ct_ga_num/ave_off_point_result_ct$all_num
ave_off_point_result_ct<-ave_off_point_result_ct[ave_off_point_result_ct$mut_freq_ct_ga>=0.01,]
ave_off_stat_ct<-ddply(ave_off_point_result_ct,.(sample),summarise,mut_reads=sum(mut_ct_ga_num))
ave_off_stat_ct<-ave_off_stat_ct[order(as.numeric(ave_off_stat_ct[,1])),]
point_result_cg<-point_result[point_result$ref%in%c("C","G"),]
point_result_cg_stat<-ddply(point_result_cg,.(sample),summarise,all_at_reads=sum(all_num))
point_result_cg_stat<-point_result_cg_stat[order(as.numeric(point_result_cg_stat$sample)),]


ave_off_stat_ct$all_at_reads<-point_result_cg_stat$all_at_reads
ave_off_stat_ct$ave_off_all<-ave_off_stat_ct$mut_reads/ave_off_stat_ct$all_at_reads

write.table(ave_off_stat_ct,"./ave_off_stat_ct.txt",quote = F,sep = "\t",col.names = T,row.names = F)








