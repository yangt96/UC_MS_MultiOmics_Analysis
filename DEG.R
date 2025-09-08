setwd('d:/')
library(GEOquery)
gset=getGEO('',destdir = '.',getGPL = F)
exprSet=exprs(gset[[1]])
exprSet=as.data.frame(exprSet)
library(limma)
boxplot(exprSet)
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet)
exprSet=as.data.frame(exprSet)
library(data.table)
anno = fread("./",data.table = F,sep = '\t')
colnames(anno)
anno = anno[,c(1,11)]
anno = anno[!anno$`Gene Symbol`== "",]
tmp = rownames(exprSet)%in% anno[,1]
exprSet = exprSet[tmp,]
dim(exprSet)
match(rownames(exprSet),anno$ID)
anno = anno[match(rownames(exprSet),anno$ID),]
match(rownames(exprSet),anno$ID)
dim(exprSet)
dim(anno)
tail(sort(table(anno[,2])), n = 12L)
{
  MAX = by(exprSet, anno[,2],
           function(x) rownames(x)[ which.max(rowMeans(x))])
  MAX = as.character(MAX)
  exprSet = exprSet[rownames(exprSet) %in% MAX,]
  rownames(exprSet) = anno[ match(rownames(exprSet), anno[,1] ),2]
}
dim(exprSet)


exprSet[1:5,1:5]

dim(exprSet)

exprSet_anno=exprSet
exprSet4=exprSet
pdata=pData(gset[[1]])
pdata$source_name_ch1 <- gsub("", "", pdata$source_name_ch1)


pdata$characteristics_ch1.2 <- gsub("", "", pdata$characteristics_ch1.2, fixed = TRUE)

group_list=pdata$`characteristics_ch1.2`
group_list1=factor(group_list,levels = c('',''))
library(limma)
design=model.matrix(~ group_list1)
n_samples <- ncol(exprSet)
n_levels <- length(levels(group_list1))
group_list <- factor(c(rep('Normal', 13), rep('Ulcerative Colitis', 30)))
design <- model.matrix(~ group_list)
fit <- lmFit(exprSet, design)
fit <- eBayes(fit) 
fit=lmFit(exprSet,design)
fit=eBayes(fit) 
allDiff_DEG=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=0.05) 
allDiff_pr1_up=allDiff_pr1[allDiff_pr1$logFC>0.5,]
allDiff_pr1_down=allDiff_pr1[allDiff_pr1$logFC< -0.5,]
save(exprSet,group_list1,allDiff_pr1_up,allDiff_pr1_down,file='')
write.table(exprSet,file='',quote = F,sep = '\t',col.names = NA)
load('')
