###input source data of raw data:-------------------------------
load("./ICB_data/dat.Gide.RData")
dat1<-Gide_data[["TPM"]]
sample_choose1<-Gide_data[["Clinical"]]
sample_choose1$group<-"Excluded"
candidate_gene<-"MALT1"
sample_choose1$group[(dat1[candidate_gene,]) <= (quantile(dat1[candidate_gene,])[3])]<-"low"
sample_choose1$group[dat1[candidate_gene,] > quantile(dat1[candidate_gene,])[3]]<-"high"

sample_choose1$status <- as.numeric(ifelse(sample_choose1$OS> 1825,0,sample_choose1$OS_CNSR))
sample_choose1$time <- ifelse(sample_choose1$OS >1825,1825,sample_choose1$OS)

fit <- survfit(Surv(time, status) ~ group, data = sample_choose1)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = F ,
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("red","blue"))

