###input raw counts
gene_count1 <- read.table("sourcedata/bulkrna2_count.csv",header = T,sep=",")
rownames(gene_count1)<-gene_count1$X
gene_count<-gene_count1[,2:dim(gene_count1)[2]]
gene_count<-as.matrix(gene_count)
#rownames(gene_count)<-gene_count1$Geneid
colnames(gene_count)<-c("Control_rep1","Control_rep2","Control_rep3",
                        "Malt1_C461A_rep1","Malt1_C461A_rep2","Malt1_C461A_rep3",
                        "Malt1_V87R_rep1","Malt1_V87R_rep2","Malt1_V87R_rep3",
                        "Malt1_L88D_rep1","Malt1_L88D_rep2","Malt1_L88D_rep3",
                        "Malt1_WT_rep1","Malt1_WT_rep2","Malt1_WT_rep3")
#x<-gene_count[unique(rownames(gene_count)),]

x<-gene_count
#Malt1_C461A_vs_Control:-------------
x<-x[,c(1:6)]
#Malt1_V87R_vs_Control:--------------
x<-x[,c(1:3,7:9)]
#Malt1_L88D_vs_Control:---------------
x<-x[,c(1:3,10:12)]
#Malt1_WT_vs_Control:-----------------
x<-x[,c(1:3,13:15)]

table <- data.frame(name =colnames(x),
                    condition =factor(sub("(^.*)_rep.*$","\\1",colnames(x))))


x<- x[rowSums(x)>0,]

dds <- DESeqDataSetFromMatrix(x, colData=table, design= ~ condition) 
rld <- rlog(dds) 
#plotPCA(rld, intgroup=c("condition"))
dds <- DESeq(dds)
res <- results(dds) 
x<-as.data.frame(res) 

norcounts <- counts(dds, normalized=T)  
y<-cbind(x,norcounts)     
deseq2_result<-y

write.table(file="deseq2.result_Malt1_C461A_vs_Control.csv",deseq2_result,row.names=T,col.names=T,quote=F,sep="\t") 
write.table(file="deseq2.result_Malt1_V87R_vs_Control.csv",deseq2_result,row.names=T,col.names=T,quote=F,sep="\t") 
write.table(file="deseq2.result_Malt1_L88D_vs_Control.csv",deseq2_result,row.names=T,col.names=T,quote=F,sep="\t") 
write.table(file="deseq2.result_Malt1_WT_vs_Control.csv",deseq2_result,row.names=T,col.names=T,quote=F,sep="\t") 


save(deseq2_result,file="deseq2_result_Malt1_C461A_vs_Control.rdata")
save(deseq2_result,file="deseq2_result_Malt1_V87R_vs_Control.rdata")
save(deseq2_result,file="deseq2_result_Malt1_L88D_vs_Control.rdata")
save(deseq2_result,file="deseq2_result_Malt1_WT_vs_Control.rdata")


###Volcano plot:-------------------------------------
threshold<-as.factor((deseq2_result$log2FoldChange > (1.5)& deseq2_result$padj<0.01))
threshold_label<-rep("no_sig",length(rownames(deseq2_result)))
threshold_label[!is.na(deseq2_result$log2FoldChange) & deseq2_result$log2FoldChange > 
                  (1.5) & !is.na(deseq2_result$padj) & deseq2_result$padj<0.01]<-
  rep("pos",length(threshold_label[!is.na(deseq2_result$log2FoldChange) & deseq2_result$log2FoldChange > 
                                     (1.5) & !is.na(deseq2_result$padj) & deseq2_result$padj<0.01]))

threshold_label[!is.na(deseq2_result$log2FoldChange) & deseq2_result$log2FoldChange < 
                  (-1.5) & !is.na(deseq2_result$padj) & deseq2_result$padj<0.01] <-
  rep("neg",length(threshold_label[!is.na(deseq2_result$log2FoldChange) & deseq2_result$log2FoldChange < (-1.5) &
                                     !is.na(deseq2_result$padj) & deseq2_result$padj<0.01]))
deseq2_result$change<-threshold_label
p<-ggplot(deseq2_result,aes(x=log2FoldChange,y=-log10(pvalue),colour=change,
                            size=threshold_label,label = rownames(deseq2_result)))+
  xlab("log2FoldChange")+ylab("-log10(padj)")+geom_point(size = 1, alpha = 0.3)+xlim(-10,10)+
  scale_color_manual(values =c("red","#0077c8","grey"),breaks=c("pos","neg","no_sig"))+
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.background = element_blank())
p


