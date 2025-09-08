
library(data.table)

  samples <- c()
  
  library(Seurat)
  library(tidyverse)
  

  CT <- Read10X(data.dir = './')
  PR <- Read10X(data.dir = './')
  

  scRNA1 <- CreateSeuratObject(counts = , project = "")
  scRNA2 <- CreateSeuratObject(counts = , project = "")
  sce.all = merge(scRNA1, y = c(scRNA2), add.cell.ids = c("", ""),
                  project = '', merge.data = TRUE)
  

  sce.all@meta.data$patient=sce.all@meta.data$orig.ident
  sce.all@meta.data$orig.ident=stringr::str_remove(sce.all@meta.data$orig.ident,'[0-9]')
  
  scRNA=sce.all
  dir.create('Pipeline')
  setwd('./Pipeline/')
  dir.create('QC')
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
  HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
  HB.genes <- HB.genes[!is.na(HB.genes)] 
  scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
  col.num <- length(levels(scRNA@active.ident))
  violin <- VlnPlot(scRNA,
                    features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                    cols =rainbow(col.num), 
                    pt.size = 0.01,
                    ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  violin
  ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
  plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
  pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
  ggsave("QC/pearplot_before_qc.pdf", plot = pearplot, width = 12, height = 5) 

  minGene=200
  maxGene=2500
  pctMT=10
  pctHB=3
  
  scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)
  
  col.num <- length(levels(scRNA@active.ident))
  violin <-VlnPlot(scRNA,
                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                   cols =rainbow(col.num), 
                   pt.size = 0.1, 
                   ncol = 4) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
  ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 
  scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  dir.create("cluster")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
  top10 <- head(VariableFeatures(scRNA), 10) 
  plot1 <- VariableFeaturePlot(scRNA) 
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
  plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
  plot
  ggsave("cluster/VariableFeatures.pdf", plot = plot, width = 8, height = 6) 

  scale.genes <-  VariableFeatures(scRNA)
  scRNA <- ScaleData(scRNA, features = scale.genes)
  
  scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
  plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
  plot2 <- ElbowPlot(scRNA, ndims=20, reduction="pca") 
  plotc <- plot1+plot2
  plotc
  ggsave("cluster/pca.pdf", plot = plotc, width = 8, height = 4) 

  pc.num=1:10
  
  scRNA <- FindNeighbors(scRNA, dims = pc.num) 
  scRNA <- FindClusters(scRNA)
  table(scRNA@meta.data$seurat_clusters)
  metadata <- scRNA@meta.data
  cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
  write.csv(cell_cluster,'cluster/cell_cluster.csv',row.names = F)
  scRNA <- RunUMAP(scRNA, dims = pc.num)
  embed_umap <- Embeddings(scRNA, 'umap')
  write.csv(embed_umap,'cluster/embed_umap.csv') 
  plot2 = DimPlot(scRNA, reduction = "umap",label = T) 
  plot2
  ggsave("cluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
  library(patchwork)
  plotc <- plot1+plot2+ plot_layout(guides = 'collect')
  plotc
  ggsave("cluster/tSNE_UMAP.pdf", plot = plotc, width = 8, height = 4)
  dir.create('cell_identify')

  library(SingleR)
  refdata <- HumanPrimaryCellAtlasData()
  
  testdata <- GetAssayData(scRNA, slot="data")
  clusters <- scRNA@meta.data$seurat_clusters
  cellpred <- SingleR(test = testdata, ref = refdata,
                      labels =refdata$label.main,
                      method = "cluster", clusters = clusters, 
                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
  
  celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
  write.csv(celltype,"cell_identify/celltype_singleR.csv",row.names = F)
  scRNA@meta.data$celltype = "NA"
  for(i in 1:nrow(celltype)){
    scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
  
  p2 = DimPlot(scRNA, group.by="celltype", label=T, label.size=4.5, reduction='umap',split.by = 'orig.ident')
  p2
  ggsave("cell_identify/UMAP_celltype.pdf", p2, width=7 ,height=6)
 Idents(scRNA)=scRNA$celltype 
  load('rf_lasso_genes.Rdata')
VlnPlot(scRNA,features = genes,group.by = 'celltype',pt.size = 0)
VlnPlot(scRNA,features = genes,group.by = 'orig.ident',pt.size = 0)


save(scRNA,file ='scRNA_anno.RDS')

table(scRNA@meta.data$celltype)

tissue_cell=paste0(scRNA@meta.data$orig.ident,"_",scRNA@meta.data$celltype)
scRNA <- AddMetaData(scRNA,tissue_cell,'tissue_cell')
table(scRNA@meta.data$tissue_cell)
ap <- AverageExpression(scRNA,features = genes,group.by ='tissue_cell' ,slot = 'data')[[1]]

rowgroup=genes
colgroup=colnames(ap)
for( i in rowgroup){
  for(j in colgroup){
    sobj <- subset(scRNA,tissue_cell==j)
    ap[i,j] <- as.numeric(table(sobj@assays$RNA@data[i,]>0)[2])/length(sobj@assays$RNA@data[i,])*100
  }}

ap=ap[,c('CT_DC','PR_DC','CT_Keratinocytes','PR_Keratinocytes','CT_Monocyte','PR_Monocyte','CT_NK_cell','CT_NK_cell','CT_T_cells','PR_T_cells')]
p=pheatmap::pheatmap(ap,display_numbers = T,  
                   color = colorRampPalette(c(rep("white",1), rep("firebrick3",1)))(100),
                   cluster_rows = F,
                   cluster_cols = F,angle_col = 45,main = 'Proportion')
p
ggsave(filename = 'pheatmap_proportion.pdf',plot = p,width = 6,height = 5)

scRNA@meta.data$tissue_cell=paste0(scRNA@meta.data$orig.ident,'_',scRNA@meta.data$celltype)
library(Seurat)

ae=AverageExpression(scRNA,assays = 'RNA',group.by = 'tissue_cell',features = genes)
ae=as.data.frame(ae$RNA)
ae=log2(ae+1)
ae=na.omit(ae)
ae=ae[,c('CT_DC','PR_DC','CT_Keratinocytes','PR_Keratinocytes','CT_Monocyte','PR_Monocyte','CT_NK_cell','CT_NK_cell','CT_T_cells','PR_T_cells')]
p=pheatmap::pheatmap(ae,display_numbers = T,  
                   color = colorRampPalette(c(rep("white",1), rep("firebrick3",1)))(100),
                   cluster_rows = F,
                   cluster_cols = F,angle_col = 45 ,main='Expression')
p
ggsave(filename = 'pheatmap_expression.pdf',plot = p,width = 6,height = 5)

p2=FeaturePlot(scRNA, features = c("SDC4", "PRDM1"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T,
            split.by = 'orig.ident',label = T) 
p2
ggsave(filename = 'SDC4_PRDM1.pdf',plot = p2,width = 10,height = 5)


p3=FeaturePlot(scRNA, features = c("NPTN"), cols =c("lightgrey", "orange"),pt.size = .9, order=T,
               split.by = 'orig.ident',label = T) 
p3
ggsave(filename = 'NPTN.pdf',plot = p3,width = 8,height = 4)

load('GS.Rdata')
scRNA <-AddModuleScore(scRNA, features= list,name = names(list))
names(x = scRNA[[]])
Cells.sub <- subset(scRNA@meta.data, celltype=='DC')
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
Idents(scRNAsub)=scRNAsub$orig.ident
DC_CT <- subset(scRNAsub, idents = c("CT"))
DC_PR <- subset(scRNAsub, idents = c("PR"))
VlnPlot(scRNAsub, features="M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34", group.by = "orig.ident", pt.size = 0,) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Oxidative score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
wilcox.test(DC_CT$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, DC_PR$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, alternative = "two.sided") #p-value < 2.2e-16
table(scRNA$celltype)
Cells.sub <- subset(scRNA@meta.data, celltype=='Keratinocytes')
scRNAsub <- subset(scRNA, cells=row.names(Cells.sub))
Idents(scRNAsub)=scRNAsub$orig.ident
Keratinocytes_CT <- subset(scRNAsub, idents = c("CT"))
Keratinocytes_PR <- subset(scRNAsub, idents = c("PR"))

VlnPlot(scRNAsub, features="M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34", group.by = "orig.ident", pt.size = 0,) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Oxidative score", x="") + theme(legend.position="right") +  
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")
wilcox.test(Keratinocytes_CT$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, Keratinocytes_PR$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34, alternative = "two.sided") #p-value < 2.2e-16

scRNA$Oxidative_score=scRNA$M5936_HALLMARK_OXIDATIVE_PHOSPHORYLATION34
p=FeaturePlot(scRNA, features = c("SDC4", "Oxidative_score"), cols =c("lightgrey", "orange", "cyan4"),pt.size = .9,  blend = TRUE, order=T,
               split.by = 'orig.ident',label = T) 
p
ggsave(filename = 'Oxidative_score.pdf',plot = p,width = 10,height = 5)
