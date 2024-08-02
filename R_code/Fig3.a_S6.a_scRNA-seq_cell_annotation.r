### 1. control sample:---------------------------------------------------------
# read raw data:---------------------------------------------------------------
control.data <- Read10X(data.dir = "./10x_out/control")
control <- CreateSeuratObject(counts = control.data, project = "control", min.cells = 3, min.features = 200)
# check qc:-----------------
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(control, features = c("nFeature_RNA"))+geom_hline(aes(yintercept=600))+geom_hline(aes(yintercept=5000))
ggsave("control_nfeatures.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
VlnPlot(control, features = c("nCount_RNA"))+geom_hline(aes(yintercept=20000))
ggsave("control_nCount.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
VlnPlot(control, features = c("percent.mt"))+geom_hline(aes(yintercept=10))
ggsave("control_mt.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
control.subset <- subset(control, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA < 20000)
# normalize and cell clustering:-------------------------
control.subset <- NormalizeData(control.subset, normalization.method = "LogNormalize", scale.factor = 10000)
control.subset <- FindVariableFeatures(control.subset, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(control.subset)
control.subset <- ScaleData(control.subset, features = all.genes)
control.subset <- RunPCA(control.subset, features = VariableFeatures(object = control.subset))
control.subset <- RunUMAP(control.subset, dims = 1:30)
control.subset <- FindNeighbors(control.subset, dims = 1:30)
control.subset <- FindClusters(control.subset, resolution = 2)
pdf("control_1.pdf")
DimPlot(control.subset,label = T)
dev.off()
# cell type annotation:-----------------------------
FeaturePlot(control.subset, ncol = 3,features = c("Ptprc",
                                           "Cd68",
                                           "S100a8",
                                           "Flt3",
                                           "Ncr1","Cd3e",
                                           "Cd8a",
                                           "Cd4",
                                           "Cd19"))

control.subset$cell_label<-"NA"
control.subset$cell_label[control.subset$seurat_clusters=="9"|
                       control.subset$seurat_clusters=="11"|
                       control.subset$seurat_clusters=="14"|
                       control.subset$seurat_clusters=="24"|
                       control.subset$seurat_clusters=="2"|
                       control.subset$seurat_clusters=="13"|
                       control.subset$seurat_clusters=="4"|
                       control.subset$seurat_clusters=="8"|
                       control.subset$seurat_clusters=="19"|
                       control.subset$seurat_clusters=="1"|
                       control.subset$seurat_clusters=="27"|
                       control.subset$seurat_clusters=="21"|
                       control.subset$seurat_clusters=="7"|
                       control.subset$seurat_clusters=="6"|
                       control.subset$seurat_clusters=="27"]<-"Mac"

control.subset$cell_label[control.subset$seurat_clusters=="16"|
                        control.subset$seurat_clusters=="23"]<-"Cd4_T"

control.subset$cell_label[control.subset$seurat_clusters=="3"|
                        control.subset$seurat_clusters=="25"]<-"Cd8_T"

control.subset$cell_label[control.subset$seurat_clusters=="29"|
                        control.subset$seurat_clusters=="12"|
                        control.subset$seurat_clusters=="0"]<-"B_cell"

control.subset$cell_label[control.subset$seurat_clusters=="28"|
                        control.subset$seurat_clusters=="10"|
                        control.subset$seurat_clusters=="20"|
                        control.subset$seurat_clusters=="18"]<-"DC"

control.subset$cell_label[control.subset$seurat_clusters=="5"]<-"NK"
control.subset$cell_label[control.subset$seurat_clusters=="26"]<-"Neu"
control.subset$cell_label[control.subset$seurat_clusters=="22"|control.subset$seurat_clusters=="17"
                      |control.subset$seurat_clusters=="15"]<-"other_T"
pdf("control_2.pdf")
DimPlot(control.subset,group.by = "cell_label",label = T)
dev.off()
### 2. PD sample:---------------------------------------------------------
# read raw data:--------------------
PD.data <- Read10X(data.dir = "./10x_out/PD")
PD <- CreateSeuratObject(counts = PD.data, project = "PD", min.cells = 3, min.features = 200)
# check qc:-----------------
PD[["percent.mt"]] <- PercentageFeatureSet(PD, pattern = "^mt-")
VlnPlot(PD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(PD, features = c("nFeature_RNA"))+geom_hline(aes(yintercept=600))+geom_hline(aes(yintercept=5000))
ggsave("PD_nfeatures.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
VlnPlot(PD, features = c("nCount_RNA"))+geom_hline(aes(yintercept=20000))
ggsave("PD_nCount.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
VlnPlot(PD, features = c("percent.mt"))+geom_hline(aes(yintercept=10))
ggsave("PD_mt.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
PD.subset <- subset(PD, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA < 20000)
# normalize and umap:-------------------------
PD.subset <- NormalizeData(PD.subset, normalization.method = "LogNormalize", scale.factor = 10000)
PD.subset <- FindVariableFeatures(PD.subset, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(PD.subset)
PD.subset <- ScaleData(PD.subset, features = all.genes)
PD.subset <- RunPCA(PD.subset, features = VariableFeatures(object = PD.subset))
PD.subset <- RunUMAP(PD.subset, dims = 1:30)
PD.subset <- FindNeighbors(PD.subset, dims = 1:30)
PD.subset <- FindClusters(PD.subset, resolution = 2)
pdf("PD_1.pdf")
DimPlot(PD.subset,label = T)
dev.off()
# cell type annotation:-----------------------------
FeaturePlot(PD.subset, ncol = 3,features = c("Ptprc",
                                             "Cd68",
                                             "S100a8",
                                             "Flt3",
                                             "Ncr1","Cd3e",
                                             "Cd8a",
                                             "Cd4",
                                             "Cd19"))
PD.subset$cell_label<-"NA"
PD.subset$cell_label[PD.subset$seurat_clusters=="25"|
                         PD.subset$seurat_clusters=="9"|
                         PD.subset$seurat_clusters=="1"|
                         PD.subset$seurat_clusters=="17"|
                         PD.subset$seurat_clusters=="4"|
                         PD.subset$seurat_clusters=="5"|
                         PD.subset$seurat_clusters=="15"|
                         PD.subset$seurat_clusters=="8"|
                         PD.subset$seurat_clusters=="12"|
                         PD.subset$seurat_clusters=="19"]<-"Mac"
PD.subset$cell_label[PD.subset$seurat_clusters=="3"|
                         PD.subset$seurat_clusters=="13"]<-"Cd4_T"
PD.subset$cell_label[PD.subset$seurat_clusters=="6"|
                         PD.subset$seurat_clusters=="14"|
                         PD.subset$seurat_clusters=="20"|
                         PD.subset$seurat_clusters=="16"]<-"Cd8_T"
PD.subset$cell_label[PD.subset$seurat_clusters=="0"|
                         PD.subset$seurat_clusters=="2"|
                         PD.subset$seurat_clusters=="18"|
                         PD.subset$seurat_clusters=="27"|
                         PD.subset$seurat_clusters=="26"]<-"B_cell"
PD.subset$cell_label[PD.subset$seurat_clusters=="10"|
                         PD.subset$seurat_clusters=="11"|
                         PD.subset$seurat_clusters=="24"|
                         PD.subset$seurat_clusters=="22"]<-"DC"
PD.subset$cell_label[PD.subset$seurat_clusters=="7"]<-"NK"
PD.subset$cell_label[PD.subset$seurat_clusters=="23"]<-"Neu"
PD.subset$cell_label[PD.subset$seurat_clusters=="21"]<-"other_T"
pdf("PD_2.pdf")
DimPlot(PD.subset,group.by = "cell_label",label = T)
dev.off()

### 3. WT sample:---------------------------------------------------------
# read raw data:--------------------
WT.data <- Read10X(data.dir = "./10x_out/WT")
WT <- CreateSeuratObject(counts = WT.data, project = "WT", min.cells = 3, min.features = 200)
# check qc:-----------------
WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^mt-")
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(WT, features = c("nFeature_RNA"))+geom_hline(aes(yintercept=600))+geom_hline(aes(yintercept=5000))
ggsave("WT_nfeatures.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
VlnPlot(WT, features = c("nCount_RNA"))+geom_hline(aes(yintercept=20000))
ggsave("WT_nCount.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
VlnPlot(WT, features = c("percent.mt"))+geom_hline(aes(yintercept=10))
ggsave("WT_mt.pdf", units="in", dpi=300, width=4, height=5, device="pdf")
WT.subset <- subset(WT, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA < 20000)
# normalize and umap:-------------------------
WT.subset <- NormalizeData(WT.subset, normalization.method = "LogNormalize", scale.factor = 10000)
WT.subset <- FindVariableFeatures(WT.subset, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(WT.subset)
WT.subset <- ScaleData(WT.subset, features = all.genes)
WT.subset <- RunPCA(WT.subset, features = VariableFeatures(object = WT.subset))
WT.subset <- RunUMAP(WT.subset, dims = 1:30)
WT.subset <- FindNeighbors(WT.subset, dims = 1:30)
WT.subset <- FindClusters(WT.subset, resolution = 2)
pdf("WT_1.pdf")
DimPlot(WT.subset,label = T)
dev.off()
FeaturePlot(WT.subset, ncol = 3,features = c("Ptprc",
                                              "Cd68",
                                              "S100a8",
                                              "Flt3",
                                              "Ncr1","Cd3e",
                                              "Cd8a",
                                              "Cd4",
                                              "Cd19"))

WT.subset$cell_label<-"NA"
WT.subset$cell_label[WT.subset$seurat_clusters=="9"|
                         WT.subset$seurat_clusters=="11"|
                         WT.subset$seurat_clusters=="14"|
                         WT.subset$seurat_clusters=="24"|
                         WT.subset$seurat_clusters=="2"|
                         WT.subset$seurat_clusters=="13"|
                         WT.subset$seurat_clusters=="4"|
                         WT.subset$seurat_clusters=="8"|
                         WT.subset$seurat_clusters=="19"|
                         WT.subset$seurat_clusters=="1"|
                         WT.subset$seurat_clusters=="27"|
                         WT.subset$seurat_clusters=="21"|
                         WT.subset$seurat_clusters=="7"|
                         WT.subset$seurat_clusters=="6"|
                         WT.subset$seurat_clusters=="27"]<-"Mac"

WT.subset$cell_label[WT.subset$seurat_clusters=="16"|
                         WT.subset$seurat_clusters=="23"]<-"Cd4_T"

WT.subset$cell_label[WT.subset$seurat_clusters=="3"|
                         WT.subset$seurat_clusters=="25"]<-"Cd8_T"

WT.subset$cell_label[WT.subset$seurat_clusters=="29"|
                         WT.subset$seurat_clusters=="12"|
                         WT.subset$seurat_clusters=="0"]<-"B_cell"

WT.subset$cell_label[WT.subset$seurat_clusters=="28"|
                         WT.subset$seurat_clusters=="10"|
                         WT.subset$seurat_clusters=="20"|
                         WT.subset$seurat_clusters=="18"]<-"DC"

WT.subset$cell_label[WT.subset$seurat_clusters=="5"]<-"NK"
WT.subset$cell_label[WT.subset$seurat_clusters=="26"]<-"Neu"
WT.subset$cell_label[WT.subset$seurat_clusters=="22"|WT.subset$seurat_clusters=="17"
                       |WT.subset$seurat_clusters=="15"]<-"other_T"
pdf("WT_2.pdf")
DimPlot(WT.subset,group.by = "cell_label",label = T)
dev.off()

### 4.Integrate the three samples of WT, PD and control:---------------------------------------
control.subset$sample<-"control"
WT.subset$sample<-"WT"
PD.subset$sample<-"PD"

reference.list <- list()
reference.list[[1]]<-control.subset
reference.list[[2]]<-WT.subset
reference.list[[3]]<-PD.subset
for (i in 1:3){
  reference.list[[i]] <- RunUMAP(reference.list[[i]], dims = 1:30)
}

list.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 3000)
list.integrated <- IntegrateData(anchorset = list.anchors, dims = 1:30)
DefaultAssay(list.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
list.integrated <- ScaleData(list.integrated, verbose = FALSE)
list.integrated <- RunPCA(list.integrated, npcs = 30, verbose = FALSE)
list.integrated <- FindNeighbors(list.integrated, dims = 1:30)
list.integrated <- FindClusters(list.integrated, resolution = 1,algorithm = 1) #original Louvain algorithm
list.integrated<-RunUMAP(list.integrated, reduction = "pca", dims = 1:30)
pdf("integrate.pdf")
DimPlot(list.integrated,group.by = "cell_label",label = T)
dev.off()
png("integrate_featureplot.png",width = 800,height = 600)
FeaturePlot(list.integrated, ncol = 3,features = c("Ptprc",
                                             "Cd68",
                                             "S100a8",
                                             "Flt3",
                                             "Ncr1","Cd3e",
                                             "Cd8a",
                                             "Cd4",
                                             "Cd19"))
dev.off()
### qc check:----------------------------------
data_qc_plot<-data.frame("cell"=colnames(list.integrated),"nFeature"=list.integrated$nFeature_RNA,
                         "nCount"=list.integrated$nCount_RNA,
                         "mt"=list.integrated$percent.mt,
                         "umap1"=list.integrated@reductions$umap@cell.embeddings[,1],
                         "umap2"=list.integrated@reductions$umap@cell.embeddings[,2])

pdf("total_cell_qc_nfeatures.pdf")
ggplot(data_qc_plot,aes(x = umap1,y = umap2)) +
  geom_point(mapping=aes(color=nFeature),size=0.01)+
  scale_color_gradient(low = "grey", high = "red", na.value = NA) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "umap1",y = "umap2", title = '') +               
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
pdf("total_cell_qc_ncounts.pdf")
ggplot(data_qc_plot,aes(x = umap1,y = umap2)) +
  geom_point(mapping=aes(color=nCount),size=0.01)+
  scale_color_gradient(low = "grey", high = "red", na.value = NA) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "umap1",y = "umap2", title = '') +               
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
pdf("total_cell_qc_mtDNA.pdf")
ggplot(data_qc_plot,aes(x = umap1,y = umap2)) +
  geom_point(mapping=aes(color=mt),size=0.01)+
  scale_color_gradient(low = "grey", high = "red", na.value = NA) +
  #scale_color_manual(values = as.character(color_new3$color))+
  ###scale_color_manual设置point的颜色和顺序
  # geom_text(aes(label = integrate_cluster), data = class_avg_seurat,label.size = 1)+
  labs(x = "umap1",y = "umap2", title = '') +               
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
pdf("total_cell_Bcell_markergene.pdf",width = 10,height = 3)
FeaturePlot(list.integrated, ncol = 3,features = c("Cd19","Cd79a" ,"Ms4a1"))
dev.off()

###In the new version of manuscript, we exluded the proliferating cells to make the results more accurate：---------------------
FeaturePlot(list.integrated, ncol = 2,features = c("Cdk1","Mki67"))
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
list.integrated <- CellCycleScoring(list.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(list.integrated,group.by ="Phase")

list.integrated <- FindClusters(list.integrated, resolution = 3,algorithm = 1) #original Louvain algorithm
DimPlot(list.integrated,group.by = "integrated_snn_res.3",label = T)
###delete proliferating cell cluster:--------------------
list.integrated<-subset(list.integrated,
                         cells=colnames(list.integrated)[!(list.integrated$integrated_snn_res.3==27|
                                                             list.integrated$integrated_snn_res.3==40|
                                                             list.integrated$integrated_snn_res.3==18)])

list.integrated<-subset(list.integrated,
                         cells=colnames(list.integrated)[!(list.integrated$integrated_snn_res.1==13)])

DefaultAssay(list.integrated) <- "integrated"
list.integrated <- RunPCA(list.integrated, npcs = 30, verbose = FALSE)
list.integrated <- FindNeighbors(list.integrated, dims = 1:30)
list.integrated <- FindClusters(list.integrated, resolution = 2,algorithm = 1) #original Louvain algorithm
list.integrated<-RunUMAP(list.integrated, reduction = "pca", dims = 1:30)
pdf("integrate.pdf")
DimPlot(list.integrated,group.by = "integrated_snn_res.3",label = T)
dev.off()
control.integrate<-subset(list.integrated,cells = colnames(list.integrated)[list.integrated$sample=="control"])
DimPlot(control.integrate,group.by = "cell_label",label =T)
WT.integrate<-subset(list.integrated,cells = colnames(list.integrated)[list.integrated$sample=="WT"])
DimPlot(WT.integrate,group.by = "cell_label",label =T)
PD.integrate<-subset(list.integrated,cells = colnames(list.integrated)[list.integrated$sample=="PD"])
DimPlot(PD.integrate,group.by = "cell_label",label =T)
###Cell ratio:-----------------------------------------
data_plot<-data.frame(table(paste(list.integrated$sample,list.integrated$cell_label,sep="_")))
data_plot$sample<-factor(c(rep("control",8),rep("PD",8),rep("WT",8)),levels = c("control","WT","PD"))
data_plot$cell_type<-sub("^.*?_(.*$)","\\1",data_plot$Var1)
data_plot$cell_por<-(data_plot$Freq/c(rep(length(colnames(control.integrate)),8),rep(length(colnames(control.integrate)),8),
                                      rep(length(colnames(control.integrate)),8)))*100

###1：
ggplot(data_plot, aes(x=sample,y=Freq,fill=cell_type)) +
  geom_bar(stat='identity',position="fill")+
  # scale_fill_manual(values = as.character(c("grey",color[5],color[8]))) +
  theme_bw()+labs(x = "Sample",y = "value", 
                  title = 'Cell proportion') +
  theme(axis.title =element_text(size = 10),
        axis.text.y =element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10,angle = 90,hjust = 0.5,vjust = 0.5,
                                   color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 10),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank())
ggsave("1.pdf", units="in", dpi=300, width=5, height=6, device="pdf")

