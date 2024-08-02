### Mac clustering: -------------------------------------------
control.subset$cell_label3<-control.integrate$cell_label3[match(colnames(control.subset),sub("_1","",colnames(control.integrate)))]
PD.subset$cell_label3<-PD.integrate$cell_label3[match(colnames(PD.subset),sub("_3","",colnames(PD.integrate)))]
WT.subset$cell_label3<-WT.integrate$cell_label3[match(colnames(WT.subset),sub("_2","",colnames(WT.integrate)))]

control.subset_mac<-subset(control.subset,cells=colnames(control.subset)[control.subset$cell_label3=="Mac"])
PD.subset_mac<-subset(PD.subset,cells=colnames(PD.subset)[PD.subset$cell_label3=="Mac"])
WT.subset_mac<-subset(WT.subset,cells=colnames(WT.subset)[WT.subset$cell_label3=="Mac"])

control.subset_mac$sample<-"control"
WT.subset_mac$sample<-"WT"
PD.subset_mac$sample<-"PD"

reference.list_mac <- list()
reference.list_mac[[1]]<-control.subset_mac
reference.list_mac[[2]]<-WT.subset_mac
reference.list_mac[[3]]<-PD.subset_mac
for (i in 1:3){
  reference.list_mac[[i]] <- RunUMAP(reference.list_mac[[i]], dims = 1:30)
}

list.anchors_mac <- FindIntegrationAnchors(object.list = reference.list_mac, dims = 1:30, anchor.features = 3000)
list.integrated_mac <- IntegrateData(anchorset = list.anchors_mac, dims = 1:30)
DefaultAssay(list.integrated_mac) <- "integrated"
# Run the standard workflow for visualization and clustering
list.integrated_mac <- ScaleData(list.integrated_mac, verbose = FALSE)
list.integrated_mac <- RunPCA(list.integrated_mac, npcs = 30, verbose = FALSE)
list.integrated_mac <- FindNeighbors(list.integrated_mac, dims = 1:30)
list.integrated_mac <- FindClusters(list.integrated_mac, resolution = 1,algorithm = 1) #original Louvain algorithm
list.integrated_mac <- RunTSNE(list.integrated_mac, dims.use = 1:30, perplexity = 10,check_duplicates = FALSE)
list.integrated_mac<-RunUMAP(list.integrated_mac, reduction = "pca", dims = 1:30)
DimPlot(list.integrated_mac,group.by = "sample",label = F,reduction = "umap")
DimPlot(list.integrated_mac,label = T)
DimPlot(list.integrated_mac,group.by = "sample",label = F,reduction = "tsne")
save(list.integrated_mac,file="list.integrated_mac.rdata")


####:-------------------------------
mac_merge<-list.integrated_mac
DefaultAssay(mac_merge)<-"RNA"
mac_merge <- NormalizeData(mac_merge, normalization.method = "LogNormalize", scale.factor = 10000)
mac_merge <- FindVariableFeatures(mac_merge, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(mac_merge)
mac_merge <- ScaleData(mac_merge, features = all.genes)
mac_merge <- RunPCA(mac_merge, features = VariableFeatures(object = mac_merge))
mac_merge <- RunUMAP(mac_merge, dims = 1:30)
mac_merge <- FindNeighbors(mac_merge, dims = 1:30)
mac_merge <- FindClusters(mac_merge, resolution = 0.5)
DimPlot(mac_merge,label = T)
DimPlot(mac_merge,label = T,group.by = "sample")
###Fig 3b:----------------------------------------------------------------
DimPlot(subset(mac_merge,cells=colnames(mac_merge)[mac_merge$sample=="control"]),label = T,group.by = "sample")
DimPlot(subset(mac_merge,cells=colnames(mac_merge)[mac_merge$sample=="WT"]),label = T,group.by = "sample")
DimPlot(subset(mac_merge,cells=colnames(mac_merge)[mac_merge$sample=="PD"]),label = T,group.by = "sample")
cluster3_marker<-FindMarkers(mac_merge, ident.1 = 3, min.pct = 0.25)
cluster3_marker[order(cluster3_marker$avg_logFC,decreasing = T)[1:20],]
rownames(cluster3_marker[order(cluster3_marker$avg_logFC,decreasing = T)[1:20],])

#IFN_TAM:----------------------------------------------------------
number<-ceiling(length(IFN_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_IFN_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = IFN_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}

###三个样本合在一起
mac_merge<-list.integrated_mac
DefaultAssay(mac_merge)<-"RNA"


mac_merge <- NormalizeData(mac_merge, normalization.method = "LogNormalize", scale.factor = 10000)
mac_merge <- FindVariableFeatures(mac_merge, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(mac_merge)
mac_merge <- ScaleData(mac_merge, features = all.genes)
mac_merge <- RunPCA(mac_merge, features = VariableFeatures(object = mac_merge))
mac_merge <- RunUMAP(mac_merge, dims = 1:30)
mac_merge <- FindNeighbors(mac_merge, dims = 1:30)
mac_merge <- FindClusters(mac_merge, resolution = 0.5)
DimPlot(mac_merge,label = T)
DimPlot(mac_merge,label = T,group.by = "sample")
DimPlot(subset(mac_merge,cells=colnames(mac_merge)[mac_merge$sample=="control"]),label = T,group.by = "sample")
DimPlot(subset(mac_merge,cells=colnames(mac_merge)[mac_merge$sample=="WT"]),label = T,group.by = "sample")
DimPlot(subset(mac_merge,cells=colnames(mac_merge)[mac_merge$sample=="PD"]),label = T,group.by = "sample")
cluster3_marker<-FindMarkers(mac_merge, ident.1 = 3, min.pct = 0.25)
cluster3_marker[order(cluster3_marker$avg_logFC,decreasing = T)[1:20],]
rownames(cluster3_marker[order(cluster3_marker$avg_logFC,decreasing = T)[1:20],])

### cell markers:-------------------------------------------------
IFN_TAM<-c("Cxcl10","Cd274","Isg15","Cd86","MHCII","Ccl2","Ccl7","Ccl15","Cxcl9","Cxcl11","Tnfsf10","Stat1",
           "Ifit1","Ifit2","Ifit3","Nos2","Rsad2","Ifitm1","Ido1","Cxcl8","Cxcl10")

Reg_TAM<-c("Arg1","Mrc1","Cxc3r1","Apoe","C1qa","Ccl2","Cd63","Clec4d","Gpnmb","Trem2","Hilpda","Hmox1","Il7r",
           "Pf4","Spp1","Vegfa","Itga4","Adgre1","Cd274")

Inflam_TAM<-c("Cxcl1","Cxcl2","Cxcl3","Cxcl5","Cxcl15","Ccl20","Ccl3l1","Il1rn","Il1b","G0s2","Inhba","Spp1","Ccl3","Cxcl1")

Angio_TAM<-c("Arg1","Adam8","Bnip3","Mif","Slc2a1","Vegfa","Areg","Cebpb")

LA_TAM<-c("Acp5","Apoc1","Apoe","C1qa","C1qb","C1qc","Ccl18","Ccl8","Cd163","Cd206","Cd36","Cd63","Ctsb","Ctsd","Ctsl","Cxcl9",
          "Fabp5","Folr2","Gpnmb","Lgals3","Macro",	"Mrc1",	"Lgals3",	"Macro",	"Trem2")


RTM_TAM<-c("Krt79",	"Krt19",	"Car4",	"Bin1",	"Nav3",	"P2ry12",	"Folr2",	"Hes1",	"Lyve1")
Prolif_TAM<-c("Cdk1",	"Mki67",	"Stmn1",	"Top2a",	"Tubb")

M1_1<-c("Il12",	"Il23",	"Tnfa",	"Il6",	"Cd86",	"MHCII",	"Il1b",	"Marco",	"Nos2",	"Il12",	"Cd64",	"Cd80",	"Cxcr10",	"Il23",	"Cxcl9",	"Cxcl10",	"Cxcl11",
        "Cd86","Il1a",	"Il6",	"Ccl5",	"Irf5",	"Irf1",	"Cd40",	"Ido1",	"Kynu",	"Ccr7")

M2_1<-c("Arg1",	"Arg2","IL10","Cd32","Cd163","Cd23","Cd200r1","Pdcd1lg2","Cd274","Marco","Csf1r","Cd206","Il1rn","Il1r2","Il4r","Ccl4","Ccl13","Ccl20","Ccl17","Ccl18","Ccl22","Ccl24","Lyve1","Vegfa",
        "Vegfb",	"Vegfc","Vegfd","Egf","Ctsa",	"Ctsb",	"Ctsc","Ctsd",	"Tgfb1","Tgfb2","Tgfb3","Mmp14","Mmp19","Mmp9","Clec7a",	"Wnt7b",	"Fasl",	"Tnfsf12",	"Tnfsf8",	"Cd276",	"Vtcn1",	"Msr1",	"Fn1",	"Irf4")


#IFN_TAM:----------------------------------------------------------
number<-ceiling(length(IFN_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_IFN_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = IFN_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(IFN_TAM,rownames(mac_merge@assays$RNA))),name = "marker")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker1"))
png(paste("IFN_TAM","_modulescore.png",sep=""),width = 513,height = 415)
ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker1))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

### Reg_TAM:-----------------------------------
number<-ceiling(length(Reg_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_Reg_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = Reg_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(Reg_TAM,rownames(mac_merge@assays$RNA))),name = "marker_Reg_TAM")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_Reg_TAM1"))

png(paste("Reg_TAM","_modulescore.png",sep=""),width = 513,height = 415)
ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_Reg_TAM1))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

###Inflam_TAM:--------------------------------------------------
number<-ceiling(length(Inflam_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_Inflam_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = Inflam_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(Inflam_TAM,rownames(mac_merge@assays$RNA))),name = "marker_Inflam_TAM")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_Inflam_TAM1"))
png(paste("Inflam_TAM","_modulescore.png",sep=""),width = 513,height = 415)

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_Inflam_TAM1))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

###Angio_TAM:--------------------------------------------------
number<-ceiling(length(Angio_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_Angio_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = Angio_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(Angio_TAM,rownames(mac_merge@assays$RNA))),name = "marker_Angio_TAM")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_Angio_TAM1"))

png(paste("Angio_TAM","_modulescore.png",sep=""),width = 513,height = 415)

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_Angio_TAM1))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

###RTM_TAM:--------------------------------------------------
number<-ceiling(length(RTM_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_RTM_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = RTM_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(RTM_TAM,rownames(mac_merge@assays$RNA))),name = "marker_RTM_TAM")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_RTM_TAM1"))
png(paste("RTM_TAM","_modulescore.png",sep=""),width = 513,height = 415)

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_RTM_TAM1))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))

dev.off()
###Prolif_TAM:--------------------------------------------------
number<-ceiling(length(Prolif_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_Prolif_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = Prolif_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(Prolif_TAM,rownames(mac_merge@assays$RNA))),name = "marker_Prolif_TAM")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_Prolif_TAM1"))

png(paste("Prolif_TAM","_modulescore.png",sep=""),width = 513,height = 415)

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_Prolif_TAM1))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

###M1_1:--------------------------------------------------
number<-ceiling(length(M1_1)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_M1_1",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = M1_1[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(M1_1,rownames(mac_merge@assays$RNA))),name = "marker_M1_1")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_M1_11"))

png(paste("M1","_modulescore.png",sep=""),width = 513,height = 415)

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_M1_11))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

###M2_1:--------------------------------------------------
number<-ceiling(length(M2_1)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_M2_1",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = M2_1[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(M2_1,rownames(mac_merge@assays$RNA))),name = "marker_M2_1")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_M2_11"))

png(paste("M2","_modulescore.png",sep=""),width = 513,height = 415)

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_M2_11))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

###LA_TAM:--------------------------------------------------
number<-ceiling(length(LA_TAM)/9)
for (j in 1:number){
  step=9
  file_name1=paste(j,"_LA_TAM",".png",sep="")
  png(file_name1,width = 750,height = 650)
  p_rna<-FeaturePlot(
    object = mac_merge,
    features = LA_TAM[((j-1)*step+1):(j*step)],
    ncol = 3,pt.size = 0.1
  )
  print(p_rna)
  dev.off()
}


mac_merge<-AddModuleScore(mac_merge,features=list(intersect(LA_TAM,rownames(mac_merge@assays$RNA))),name = "marker_LA_TAM")
mydata<-FetchData(mac_merge,vars=c("UMAP_1","UMAP_2","marker_LA_TAM1"))

png(paste("LA_TAM","_modulescore.png",sep=""),width = 513,height = 415)

ggplot(mydata,aes(x=UMAP_1,y=UMAP_2,colour=marker_LA_TAM1))+geom_point(size=1)+
  scale_color_gradientn(values=seq(0,1,0.2),colours=c("blue","cyan","green","yellow","orange","red")) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "TSNE1",y = "TSNE2", title = '') +               
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 15, color = 'black'),
        legend.title = element_text(size = 15, color = 'black'))
dev.off()

###annotation cells:----------------------------------------------
mac_merge$annotation<-""
mac_merge$annotation[mac_merge$seurat_clusters=="0"]<-"C1q_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="1"]<-"Cxcl10_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="2"]<-"Fn1_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="3"]<-"Ifitm1_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="4"]<-"NA"
mac_merge$annotation[mac_merge$seurat_clusters=="5"]<-"Spp1_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="6"]<-"Lpl_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="7"]<-"NA"
mac_merge$annotation[mac_merge$seurat_clusters=="8"]<-"Ace_monocytes"
mac_merge$annotation[mac_merge$seurat_clusters=="9"]<-"Areg_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="10"]<-"Arg1_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="11"]<-"Cd74_macrophage"
mac_merge$annotation[mac_merge$seurat_clusters=="12"]<-"Vcam1_macrophage"
mac_merge2<-subset(mac_merge,cells=colnames(mac_merge)[!(mac_merge$annotation=="NA")])

###cell number:-----------------------------------
data_plot_mac<-data.frame(table(paste(mac_merge2$sample,mac_merge2$annotation,sep="_")))
data_plot_mac[33,1]<-c("control_Cd74_macrophage")
data_plot_mac[33,2]<-0


data_plot_mac$Var1<-factor(c(as.character(data_plot_mac$Var1[-33]),"control_Cd74_macrophage"),levels = c(as.character(factor(data_plot_mac$Var1))[-33],"control_Cd74_macrophage"))
data_plot_mac[33,1]<-c("control_Cd74_macrophage")

data_plot_mac$cell_type<-as.factor(sub("^.*?_(.*$)","\\1",data_plot_mac$Var1))
data_plot_mac$sample<-factor(c(rep("control",10),rep("PD",11),rep("WT",11),"control"),levels = c("control","WT","PD"))

data_plot_mac$cell_por<-(data_plot_mac$Freq/c(rep(512,10),rep(3514,11),
                                              rep(5108,11),512))*100


###Fig S6.b：--------------------------------------------------
ggplot(data=data_plot_mac, aes(x=cell_type,y=Freq,fill=sample)) +geom_bar(stat = 'identity', position = 'dodge')+scale_y_continuous(expand = c(0,0))+
  labs(x = "Cluster",y = "cell number") +theme_bw()+
  scale_fill_manual(values = as.character(c("grey","red","lightblue"))) +scale_y_continuous(expand = c(0,0)) +
  theme(axis.title =element_text(size = 10),
        axis.text.y =element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10,angle = 90,hjust = 0.5,vjust = 0.5,
                                   color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 10),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())
ggsave("11.pdf", units="in", dpi=300, width=6, height=6, device="pdf")


####Fig3.b:---------------------------------------------------------
DimPlot(subset(mac_merge2,cells=colnames(mac_merge2)[mac_merge2$sample=="control"]),label = F,group.by = "annotation")+
  coord_fixed(ratio = (max(mac_merge2@reductions$umap@cell.embeddings[,1])-min(mac_merge2@reductions$umap@cell.embeddings[,1]))/(max(mac_merge2@reductions$umap@cell.embeddings[,2])-min(mac_merge2@reductions$umap@cell.embeddings[,2])))
ggsave("7.pdf", units="in", dpi=300, width=6, height=6, device="pdf")

DimPlot(subset(mac_merge2,cells=colnames(mac_merge2)[mac_merge2$sample=="WT"]),label =F,group.by = "annotation")+
  coord_fixed(ratio = (max(mac_merge2@reductions$umap@cell.embeddings[,1])-min(mac_merge2@reductions$umap@cell.embeddings[,1]))/(max(mac_merge2@reductions$umap@cell.embeddings[,2])-min(mac_merge2@reductions$umap@cell.embeddings[,2])))
ggsave("8.pdf", units="in", dpi=300, width=6, height=6, device="pdf")

DimPlot(subset(mac_merge2,cells=colnames(mac_merge2)[mac_merge2$sample=="PD"]),label = F,group.by = "annotation")+
  coord_fixed(ratio = (max(mac_merge2@reductions$umap@cell.embeddings[,1])-min(mac_merge2@reductions$umap@cell.embeddings[,1]))/(max(mac_merge2@reductions$umap@cell.embeddings[,2])-min(mac_merge2@reductions$umap@cell.embeddings[,2])))
ggsave("9.pdf", units="in", dpi=300, width=6, height=6, device="pdf")










