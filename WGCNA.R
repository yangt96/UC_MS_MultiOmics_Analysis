
library(limma)
setwd('d:/')
library(sva)
mergeFile="merge.preNorm.txt"            
normalizeFile="merge.normalzie.txt"      
files=c('','')      
geneList=list()
for(i in 1:length(files)){
  fileName=files[i]
  rt=read.table(fileName, header=T, sep="\t", check.names=F)
  header=unlist(strsplit(fileName, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  fileName=files[i]
  header=unlist(strsplit(fileName, "\\.|\\-"))
  rt=read.table(fileName, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))

 
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(header[1],ncol(rt)))
}

allTabOut=rbind(geneNames=colnames(allTab), allTab)
write.table(allTabOut, file=mergeFile, sep="\t", quote=F, col.names=F)


normalizeTab=ComBat(allTab, batchType, par.prior=TRUE)
normalizeTab=rbind(geneNames=colnames(normalizeTab), normalizeTab)
write.table(normalizeTab, file=normalizeFile, sep="\t", quote=F, col.names=F)


library(ggplot2)       

rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))


data.pca=prcomp(data)
pcaPredict=predict(data.pca)


bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]



p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
ggsave(p,file = 'PCA_pre.pdf',width = 5.1,height = 3.8)


rt=read.table('merge.normalzie.txt', header=T, sep="\t", check.names=F, row.names=1)
data=t(rt)
Project=gsub("(.*?)\\_.*", "\\1", rownames(data))
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))


data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

bioCol=c("#bc2c29",'#0072b5',"#20854e","#ef7c1c","#EE4C97","#FF9900","#20854E","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121")
bioCol=bioCol[1:length(levels(factor(Project)))]



p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = Type)) +
  scale_colour_manual(name="",  values=bioCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(p)
ggsave(p,file = 'PCA_norm.pdf',width = 5.1,height = 3.8)



library(WGCNA)
options(stringsAsFactors = FALSE)
setwd('d:/')
exprSet = read.table("merge.normalzie.txt",header = T,check.names = F,row.names = 1)

A=grep('',colnames(exprSet))
B=grep('',colnames(exprSet))
exprSet=exprSet[,c(A,B)]

load('')
load('')


group_list=c(group_list1,group_list3)

anno=data.frame(row.names=colnames(exprSet),group=group_list)

WGCNA_matrix = t(exprSet[order(apply(exprSet,1,mad), decreasing = T)[1:nrow(exprSet)],])
datExpr0 <- WGCNA_matrix  
datExpr0 <- as.data.frame(datExpr0)

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK 

if (!gsg$allOK)
{
  
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

  sampleTree = hclust(dist(datExpr0), method = "average")
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)


clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust) 
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]


anno=anno[rownames(datExpr),,drop=F]


sampleTree = hclust(dist(datExpr), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers (after)", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

anno$=ifelse(anno$group=='',1,0)
anno$=ifelse(anno$group=='',1,0)

datTraits =anno[,-1]
datExpr=datExpr[rownames(datTraits),]
sampleNames = rownames(datExpr)

traitRows = match(sampleNames, rownames(datTraits))  


enableWGCNAThreads()   
powers = c(1:30)       
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,blockSize = 100000)

dev.off()
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red") 

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

sft 
softPower =sft$powerEstimate 

softPower
adjacency = adjacency(datExpr, power = softPower)

net = blockwiseModules(datExpr, power = softPower,
                       TOMType = "unsigned", minModuleSize = 60,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)

table(net$colors)


mergedColors = labels2colors(net$colors)


plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]


nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)



textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


identical(rownames(datTraits),rownames(datExpr))


PR = as.data.frame(datTraits$Psoriasis)
names(PR) = "Psoriasis"

geneTraitSignificance = as.data.frame(cor(datExpr, PR, use = "p"))

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(PR), sep="")
names(GSPvalue) = paste("p.GS.", names(PR), sep="")  
module = "blue"
column = match(module, modNames) 
moduleGenes = moduleColors==module
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for MS",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



MS = as.data.frame(datTraits$MS)
names(MS) = "MS"

geneTraitSignificance = as.data.frame(cor(datExpr, MS, use = "p"))

GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(MS), sep="")
names(GSPvalue) = paste("p.GS.", names(MS), sep="")  

module = "blue"
column = match(module, modNames) 
moduleGenes = moduleColors==module
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for MS",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


module = "blue"
probes = names(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]
modProbes
write.table(modProbes,file ='blue_gene.txt',row.names = F,col.names = F,quote=F)
