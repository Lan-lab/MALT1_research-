###input data:------------------
gene133<-read.table("133_gene.txt",header=T)
load("~/other_project/taoyuwei/ICB_result/dat.Kim.RData")
kim_data<-dat[["mRNA.norm3"]]
re<-colnames(kim_data)[dat[["response"]]==1]
non_re<-colnames(kim_data)[dat[["response"]]==0]

rownames(kim_data)<-dat$genes
kim_data2<-kim_data[(!is.na(rownames(kim_data)) & !(rownames(kim_data)=="")),]
kim <- CreateSeuratObject(counts = as.matrix(kim_data2), project = "kim", min.cells = 0, min.features = 1)
kim <- NormalizeData(kim, normalization.method = "LogNormalize", scale.factor = 10000)
#kim <- FindVariableFeatures(kim, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(kim)
kim <- ScaleData(kim, features = all.genes)
kim@assays$RNA$data<-kim@assays$RNA$counts
genelist<-list()
genelist[[1]]<-intersect(rownames(kim_data),toupper(gene133$Gene))
kim<-AddModuleScore(
  kim,
  features=genelist,
  name = "gene133",
  seed = 1,
  search = FALSE,
  slot = "data"
)
plot<-data.frame("data"=kim$gene1331,
                 "sample"=colnames(kim))
plot$group<-"NA"
plot$group[match(re,plot$sample)]<-"re"
plot$group[match(non_re,plot$sample)]<-"non_re"
plot$group<-factor(plot$group,levels = c("re","non_re"))
p_value_1<-t.test(plot$data[plot$group=="re"],plot$data[plot$group=="non_re"])$p.value

pdf(paste("133gene_kimdataset",".pdf",sep=""),width=3,height=4)
p<-ggpubr::ggboxplot(plot, 
                     x="group", y="data",title =paste("133_gene",p_value_1,sep=":"),
                     color = "group", add = "jitter")+ 
  ggplot2::scale_color_manual(values = c("#cc0000","#8b8c8d"))+
  ggpubr::rotate_x_text(angle = 45)+ 
  #geom_hline(yintercept = mean(mydata$proportion), linetype=2)+# Add horizontal line at base mean 
  ggpubr::stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = "Tumor")# Pairwise comparison against all                                                                                                                 

p
dev.off()
write.csv(plot,"./ICB_result/gene133_gene_score_kimdata.csv",quote = F)
