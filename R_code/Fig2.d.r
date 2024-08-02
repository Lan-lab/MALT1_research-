### kim data input:---------------------------
kim_data<-as.matrix(read.table("kim_data_raw.csv",header = T,sep=","))
rownames(kim_data)<-kim_data[,1]
kim_data<-kim_data[,2:dim(kim_data)[2]]
sample<-read.table(file="kim_sample_info.csv",header = T,sep=",")
re<-colnames(kim_data)[sample$x==1]
non_re<-colnames(kim_data)[sample$x==0]

plot<-data.frame("data"=c(as.numeric(kim_data["MALT1",intersect(colnames(kim_data),re)]),as.numeric(kim_data["MALT1",intersect(colnames(kim_data),non_re)])),
                 "group"=c(rep("Re",length(intersect(colnames(kim_data),re))),
                           rep("Non-Re",length(intersect(colnames(kim_data),non_re)))))

p_value_1<-t.test(as.numeric(kim_data["MALT1",intersect(colnames(kim_data),re)]),as.numeric(kim_data["MALT1",intersect(colnames(kim_data),non_re)]))$p.value

pdf("2.pdf",width=3,height=4)
ggpubr::ggboxplot(plot, 
                  x="group", y="data",title =paste("p_value",p_value_1,sep=":"),
                  color = "group", add = "jitter")+ 
  ggplot2::scale_color_manual(values = c("#8b8c8d","#cc0000"))+
  ggpubr::rotate_x_text(angle = 45)+ 
  #geom_hline(yintercept = mean(mydata$proportion), linetype=2)+# Add horizontal line at base mean 
  ggpubr::stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = "Tumor")# Pairwise comparison against all                                                                                                                 

dev.off()
write.csv(plot,file="icb_MALT1_expression_Kim.csv")



###input Gide raw data:----------------------------------------------
Gide_data<-as.matrix(read.table("Gide_data_raw.csv",header = T,sep=","))
rownames(Gide_data)<-Gide_data[,1]
Gide_data<-Gide_data[,2:dim(Gide_data)[2]]
sample<-read.table(file="Gide_sample_info.csv",header = T,sep=",")
re<-colnames(Gide_data)[sample$x=="CR"]
non_re<-colnames(Gide_data)[sample$x=="PR"]
plot<-data.frame("data"=c(as.numeric(Gide_data["MALT1",intersect(colnames(Gide_data),re)]),as.numeric(Gide_data["MALT1",intersect(colnames(Gide_data),non_re)])),
                 "group"=c(rep("Re",length(intersect(colnames(Gide_data),re))),
                           rep("Non-Re",length(intersect(colnames(Gide_data),non_re)))))

p_value_1<-t.test(as.numeric(Gide_data["MALT1",intersect(colnames(Gide_data),re)]),as.numeric(Gide_data["MALT1",intersect(colnames(Gide_data),non_re)]))$p.value

pdf("MALT1_expression_NR_vs_R_Gide.pdf",width=3,height=4)
ggpubr::ggboxplot(plot, 
                  x="group", y="data",title =paste("p_value",p_value_1,sep=":"),
                  color = "group", add = "jitter")+ 
  ggplot2::scale_color_manual(values = c("#8b8c8d","#cc0000"))+
  ggpubr::rotate_x_text(angle = 45)+ 
  #geom_hline(yintercept = mean(mydata$proportion), linetype=2)+# Add horizontal line at base mean 
  ggpubr::stat_compare_means(label = "p.format", method = "wilcox.test", ref.group = "Tumor")# Pairwise comparison against all                                                                                                                 

dev.off()
write.csv(plot,file="icb_MALT1_expression_gide.csv")
