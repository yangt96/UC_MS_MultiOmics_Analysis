
gc()
dev.off()

scRNA=readRDS('./')
load('')
gene=unique(table1$Gene)
gene
library(Seurat)
library(viridis)
DotPlot(scRNA,features = gene,cols = c('#dadada','#bc3c29'))
scRNA_T=readRDS('./')
library(viridis)
FeaturePlot(scRNA_T,features = '',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'),split.by = 'tissue_type')
FeaturePlot(scRNA,features = '',label = T,pt.size = 0.5,order = T,cols = c('#dadada','#bc3c29'),split.by = 'tissue_type')

library(mixtools)

library(GeneSwitches)
library(SingleCellExperiment)


scRNA_cm=subset(scRNA_T,T_celltype=='')

cellinfo <- scRNA_cm@meta.data

sce <- as.SingleCellExperiment(scRNA_cm)

library(slingshot)
library(RColorBrewer)
library(SingleCellExperiment)
library(Seurat)

sce_slingshot <- slingshot(sce , clusterLabels = 'T_celltype', reducedDim = 'UMAP', 
                           start.clus = c(3,5), shrink = 0.2)





dev.off()

cl1 <- cellinfo$T_celltype
plot(reducedDims(sce_slingshot)$UMAP,col = brewer.pal(12,"Paired")[cl1],pch=16,asp=1)


igraph::igraph.options(sparsematrices = FALSE)

lines(SlingshotDataSet(sce_slingshot), lwd=2, type = 'lineages', col = 'black')


legend("right",legend = unique(sce$T_celltype),
       col = unique(brewer.pal(12,"Paired")[cl1]),inset=c(3,2,4), pch = 16)

allexpdata <- as.matrix(scRNA_cm@assays$RNA@data);dim(allexpdata)
allcells<-colData(sce_slingshot);dim(allcells)

allcells$slingPseudotime_1

cells <- allcells[!is.na(allcells$slingPseudotime_1),];dim(cells)
expdata <- allexpdata[,rownames(cells)];dim(expdata)

expdata <- expdata[apply(expdata > 0,1,sum) >= 5,];dim(expdata)

rd_UMAP <- Embeddings(object = scRNA_cm, reduction = "umap");dim(rd_UMAP)#原 object = seu3obj.integrated
rd_UMAP <- rd_UMAP[rownames(cells), ];dim(rd_UMAP)
all(rownames(rd_UMAP) == colnames(expdata))

library(mixtools)
library(GeneSwitches)
library(SingleCellExperiment)



sce <- SingleCellExperiment(assays = List(expdata = expdata))

colData(sce)$Pseudotime <- cells$slingPseudotime_1

reducedDims(sce) <- SimpleList(UMAP = rd_UMAP)
sce_p1 <- sce

h <- hist(assays(sce_p1)$expdata, breaks = 200, plot = FALSE)
{plot(h, freq = FALSE, xlim = c(0,2), ylim = c(0,1), main = "Histogram of gene expression",
      xlab = "Gene expression", col = "darkgoldenrod2", border = "grey")
  abline(v=0.2, col="blue")}


sce_p1 <- binarize_exp(sce_p1, fix_cutoff = TRUE, binarize_cutoff = 0.2)

sce_p1 <- find_switch_logistic_fastglm(sce_p1, downsample = TRUE, show_warning = FALSE, zero_ratio = 0.65, ds_cutoff = 0.65)

table(rowData(sce_p1)$prd_quality)

sg_allgenes <- filter_switchgenes(sce_p1, allgenes = TRUE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92);dim(sg_allgenes)


sg_gtypes <- filter_switchgenes(sce_p1, allgenes = FALSE, r2cutoff = 0.01, topnum = 25, zero_pct = 0.92,
                                genelists = gs_genelists);dim(sg_gtypes)#, genetype = c("Surface proteins", "TFs"))

sg_vis <- rbind(sg_gtypes, sg_allgenes[setdiff(rownames(sg_allgenes), rownames(sg_gtypes)),]);dim(sg_vis)
gl=gene

intersect(sg_vis$geneID, gl)
sg_my <- rowData(sce_p1)[gl,];head(sg_my)
sg_my$feature_type <- "Mendelian genes"
sg_vis <- rbind(sg_vis, sg_my)
plot_timeline_ggplot(sg_vis, timedata = sce_p1$Pseudotime, txtsize = 3.5)


a=sce_p1@assays@data$expdata['',] 
b=sce_p1$Pseudotime

df=data.frame(gene=a,time=b)

ggstatsplot::ggscatterstats(data=df,x='time',y='gene')
dev.off()

gc()

scRNA=readRDS('./')
scRNA_other=subset(scRNA,celltype != 'T_cells')
rm(scRNA)
gc()
scRNA_T=readRDS('./')


gc()


scRNA_CM=subset(scRNA_T,T_celltype=='CD4_REG')
scRNA_CM$gene_group=ifelse(scRNA_CM@assays$RNA@counts['PDIA6',]>0,'PDIA6+REG','PDIA6-REG')

scRNA_otherT=subset(scRNA_T,T_celltype != 'CD4_REG')

scRNA_other$gene_group =scRNA_other$celltype
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_chat=merge(scRNA_CM,c(scRNA_other,scRNA_otherT))

rm(scRNA_CM,scRNA_otherT)
rm(scRNA_T)
rm(scRNA_other)
gc()

scRNA_chat_SLE=subset(scRNA_chat,tissue_type=='UC')

gc()
set.seed(123)
a=sample(1:ncol(scRNA_chat_SLE),6000)
scRNA_chat_SLE=scRNA_chat_SLE[,a]

meta =scRNA_chat_SLE@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat_SLE@assays$RNA@data)
identical(colnames(data_input),rownames(meta))

library(CellChat)
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "gene_group")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)

cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= T, sources.use = c('PDIA6+REG','PDIA6-REG'),
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           targets.use = c('PDIA6+REG','PDIA6-REG'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble


scRNA_T=readRDS('./')

gc()
scRNA_CM=subset(scRNA_T,T_celltype=='')
scRNA_CM$gene_group=ifelse(scRNA_CM@assays$RNA@counts['',]>0,'','')

scRNA_otherT=subset(scRNA_T,T_celltype != '')
scRNA_otherT$gene_group=scRNA_otherT$T_celltype

scRNA_metab=merge(scRNA_CM,c(scRNA_otherT))

gc()

rm(scRNA_chat)
gc()
scRNA_metab_SLE=subset(scRNA_metab,tissue_type=='UC')

set.seed(123)
a=sample(1:ncol(scRNA_metab_SLE),2000)
scRNA_metab_SLE=scRNA_metab_SLE[,a]



library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA_metab_SLE<-sc.metabolism.Seurat(obj = scRNA_metab_SLE, method = 'AUCell', imputation = F, ncores = 2, metabolism.type = "KEGG")

input.pathway <- rownames(scRNA_metab_SLE@assays[["METABOLISM"]][["score"]])[61:90]
DotPlot.metabolism(obj =scRNA_metab_SLE,
                   pathway = input.pathway, phenotype = "gene_group", norm = "y")


gc()

library(Seurat)

Idents(scRNA_CM)=scRNA_CM$gene_group
df=FindAllMarkers(scRNA_CM,only.pos = T,logfc.threshold =0.5)
write.csv(df,'IFIT3_marker.csv',quote = F)

gc()

library(data.table)
rt=fread('./GSE112087_counts-matrix-EnsembIDs-GRCh37.p10.txt',data.table = F)
rownames(rt)=rt$V1
rt$V1=NULL


