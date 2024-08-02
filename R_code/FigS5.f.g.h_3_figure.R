### input data1:COAD data--------------------------------
tcga_brca1<-read.table("~/fshare_file_move/others_project_money/taoyuwei/cellphonedb/other_tumor_type_score/TCGA-COAD.htseq_fpkm.tsv",header = T,sep="\t")
load("~/fshare_file_move/others_project_money/taoyuwei/cellphonedb/other_tumor_type_score/gene_name_convert.rdata")
rownames(tcga_brca1)<-tcga_brca1$Ensembl_ID
tcga_brca1<-tcga_brca1[match(intersect(sub("(^.*)\\..*$","\\1",rownames(tcga_brca1)),gene_name_convert$ENSEMBL),sub("(^.*)\\..*$","\\1",rownames(tcga_brca1))),]
symbol<-gene_name_convert$SYMBOL[match(sub("(^.*)\\..*$","\\1",rownames(tcga_brca1)),gene_name_convert$ENSEMBL)]
tcga_brca1<-tcga_brca1[match(unique(symbol),symbol),]
rownames(tcga_brca1)<-unique(symbol)
tcga_brca1<-tcga_brca1[,2:dim(tcga_brca1)[2]]

#### input dataset2 TNBC:------------------------------
tcga_brca2<-read.table("~/fshare_file_move/others_project_money/taoyuwei/tcga/FUSCCTNBC_Expression_RNAseqFPKM.csv",header = T,sep=",")
rownames(tcga_brca2)<-tcga_brca2$gene
tcga_brca2<-tcga_brca2[,2:dim(tcga_brca2)[2]]


#### input dataset3:SKCM:--------------------------------------
tcga_brca3<-read.table("~/fshare_file_move/others_project_money/taoyuwei/cellphonedb/other_tumor_type_score/TCGA-SKCM.htseq_fpkm.tsv",header = T,sep="\t")
load("~/fshare_file_move/others_project_money/taoyuwei/cellphonedb/other_tumor_type_score/gene_name_convert.rdata")
rownames(tcga_brca3)<-tcga_brca3$Ensembl_ID
tcga_brca3<-tcga_brca3[match(intersect(sub("(^.*)\\..*$","\\1",rownames(tcga_brca3)),gene_name_convert$ENSEMBL),sub("(^.*)\\..*$","\\1",rownames(tcga_brca3))),]
symbol<-gene_name_convert$SYMBOL[match(sub("(^.*)\\..*$","\\1",rownames(tcga_brca3)),gene_name_convert$ENSEMBL)]
tcga_brca3<-tcga_brca3[match(unique(symbol),symbol),]
rownames(tcga_brca3)<-unique(symbol)
tcga_brca3<-tcga_brca3[,2:dim(tcga_brca3)[2]]


####Plot PDL1 expression:3 dataset----------------------------------
use_mat<-as.matrix(tcga_brca1)
use_mat<-as.matrix(tcga_brca2)
use_mat<-as.matrix(tcga_brca3)

####Plot PDL1 expression:3 dataset----------------------------------
use_mat<-as.matrix(tcga_brca3)

low_group<-colnames(use_mat)[(use_mat[grep("^MALT1$",rownames(use_mat)),]) <= (quantile(use_mat[grep("^MALT1$",rownames(use_mat)),])[2])]
high_group<-colnames(use_mat)[(use_mat[grep("^MALT1$",rownames(use_mat)),]) >= (quantile(use_mat[grep("^MALT1$",rownames(use_mat)),])[4])]
gene_plot_data<-use_mat[intersect("CD274",rownames(use_mat)),
                        intersect(c(high_group,low_group),colnames(use_mat))]

cd274_score<-gene_plot_data
pdf("Cd274_colon.pdf",width=5,height=6)
data_plot<-data.frame("group"=c(rep("high_group",length(high_group)),rep("low_group",length(low_group))),"cd274_score"=cd274_score)
p_value_1<-wilcox.test(data_plot$cd274_score[data_plot$group=="high_group"],data_plot$cd274_score[data_plot$group=="low_group"])$p.value
plot<-data.frame("group"=c(rep("high_group",length(high_group)),rep("low_group",length(low_group))),"cd274_score"=cd274_score)
plot$group<-factor(plot$group,levels = c("low_group","high_group"))
ggpubr::ggboxplot(plot, 
                  x="group", y="cd274_score",title =paste("p_value",p_value_1,sep=":"),
                  color = "group", add = "jitter")+ 
  ggplot2::scale_color_manual(values = c("#8b8c8d","#cc0000"))+
  ggpubr::rotate_x_text(angle = 45)+ 
  #geom_hline(yintercept = mean(mydata$proportion), linetype=2)+# Add horizontal line at base mean 
  ggpubr::stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = "Tumor")   # Pairwise comparison against all                                                                                                                 

dev.off()
write.csv(use_mat,file="COAD_Cd274_expression_MALT1group.csv",quote = F)
write.csv(use_mat,file="TNBC_Cd274_expression_MALT1group.csv",quote = F)
write.csv(use_mat,file="SKCM_Cd274_expression_MALT1group.csv",quote = F)

