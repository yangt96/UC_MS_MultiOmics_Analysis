
library("limma")         

source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "merge.normalzie.txt", perm=100, QN=F)

library(limma)
library(vioplot)
immuneFile=""     

pFilter=0.05          

immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

gene=read.csv('')
gene=gene$GENE
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

rt=data[gene,]
load('')
a=grep('',colnames(rt))
b=grep('',colnames(rt))
rt_pr=rt[,c(a,b)]
group1=ifelse(group_list1=='','','')
group2=ifelse(group_list4=='','' )
group=c(group1,group2)

anno=data.frame(row.names = colnames(rt_pr),group=group)

ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_pr=rt_pr[,ss]

group=anno$group

lowName=row.names(anno)[anno[,1]=='']
highName=row.names(anno)[anno[,1]=='']

lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)
outTab=data.frame()
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
bioCol=c("#0072b5","#bc3c29","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
for(i in 1:ncol(rt)){
  if(sd(rt[1:lowNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
    rt[(lowNum+1),i]=0.00001
  }
  lowData=rt[1:lowNum,i]
  highData=rt[(lowNum+1):(lowNum+highNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col=bioCol[1])
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col=bioCol[2])
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("", ""),
       lwd=4.5,bty="n",cex=1.5,
       col=bioCol[1:2])
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)

write.table(outTab,file="diff_result_pr.txt",sep="\t",row.names=F,quote=F)
library(Hmisc)
library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"   

pFilter=0.05         

immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])
gene=read.csv('')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

rt=data[gene,]

a=grep('',colnames(rt))
b=grep('',colnames(rt))
rt_pr=rt[,c(a,b)]
group1=ifelse(group_list1=='')
group2=ifelse(group_list2=='')
group=c(group1,group2)

anno=data.frame(row.names = colnames(rt_pr),group=group)

ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_pr=rt_pr[,ss]

nc =cbind(immune,t(rt_pr))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
m=as.data.frame(m)
m =dplyr::filter(m,m$CYP2W1 != 'NaN')
p = rcorr(nc)$P[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
p =p[rownames(m),]
library(dplyr)
tmp <- matrix(ifelse(p < 0.001, '***',
                     ifelse(p < 0.01, '**',
                            ifelse(p < 0.05, '*', ''))), nrow = nrow(p))


library(pheatmap)
pheatmap(t(m),
         display_numbers =t(tmp),
         angle_col =45,
         color = colorRampPalette(c("#0072b5", "white", "#bc3c29"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)

library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"      
pFilter=0.05          

immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

gene=read.csv('')
gene=gene$GENE
data=read.table('',header = T,check.names = F,row.names = 1)

rt=data[gene,]


a=grep('',colnames(rt))
rt_MS=rt[,a]
anno=data.frame(row.names = colnames(rt_MS),group=group_list3)

ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_MS=rt_MS[,ss]

group=anno$group
lowName=row.names(anno)[anno[,1]=='']
highName=row.names(anno)[anno[,1]=='']
lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)

outTab=data.frame()
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
bioCol=c("#20854e","#e18727","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
for(i in 1:ncol(rt)){
  if(sd(rt[1:lowNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
    rt[(lowNum+1),i]=0.00001
  }
  lowData=rt[1:lowNum,i]
  highData=rt[(lowNum+1):(lowNum+highNum),i]
  vioplot(lowData,at=3*(i-1),lty=1,add = T,col=bioCol[1])
  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col=bioCol[2])
  wilcoxTest=wilcox.test(lowData,highData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(lowData,highData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("CT", "MS"),
       lwd=4.5,bty="n",cex=1.5,
       col=bioCol[1:2])
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)


library(Hmisc)

library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"   

pFilter=0.05          

immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)

immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

gene=read.csv('')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

rt=data[gene,]


a=grep('',colnames(rt))
rt_MS=rt[,a]
anno=data.frame(row.names = colnames(rt_MS),group=group_list3)

ss=intersect(rownames(anno),rownames(immune))
anno=anno[ss,,drop=F]
immune=immune[ss,]
rt_MS=rt_MS[,ss]

nc =cbind(immune,t(rt_MS))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
m = as.data.frame(m)
m =dplyr::filter(m,m$CYP2W1 != 'NaN')
p = rcorr(nc)$P[1:ncol(immune),(ncol(nc)-length(gene)+1):ncol(nc)]
library(dplyr)
tmp = matrix(case_when(p<0.001~'***',
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))


pheatmap(t(m),
         display_numbers =t(tmp),
         angle_col =45,
         color = colorRampPalette(c("#20854e", "white", "#e18727"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)


load('GS.Rdata')

library(GSVA)
library(Biobase)

data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)

uni_matrix=data

list= list


gsva_matrix<- gsva(as.matrix(uni_matrix), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

rownames(gsva_matrix)=sub('^......','',rownames(gsva_matrix))

save(gsva_matrix,file = 'gsva_array.Rdata')
load('gsva_array.Rdata')

metabolism=gsva_matrix[c('HALLMARK_HYPOXIA',
                         'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
                         'HALLMARK_ADIPOGENESIS',
                         'HALLMARK_XENOBIOTIC_METABOLISM',
                         'HALLMARK_FATTY_ACID_METABOLISM',
                         'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                         'HALLMARK_GLYCOLYSIS',
                         'HALLMARK_HEME_METABOLISM',
                         'HALLMARK_BILE_ACID_METABOLISM'),]

gene=read.csv('')
gene=gene$GENE
rt=data[gene,]

a=grep('',colnames(rt))
b=grep('',colnames(rt))
rt_pr=rt[,c(a,b)]
group1=ifelse(group_list1=='')
group2=ifelse(group_list2=='' )
group=c(group1,group2)

anno=data.frame(row.names = colnames(rt_pr),group=group)

metabolism_pr=metabolism[,rownames(anno)]
metabolism_pr=t(metabolism_pr)

library(Hmisc)
nc =cbind(metabolism_pr,t(rt_pr))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(metabolism_pr),(ncol(nc)-length(gene)+1):ncol(nc)]
rownames(m)=stringr::str_remove(rownames(m),pattern = 'HALLMARK_')
p = rcorr(nc)$P[1:ncol(metabolism_pr),(ncol(nc)-length(gene)+1):ncol(nc)]
library(dplyr)
tmp = matrix(case_when(p<0.001~'***',
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))


library(pheatmap)
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#0072b5", "white", "#bc3c29"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)


load('gsva_array.Rdata')

metabolism=gsva_matrix[c('HALLMARK_HYPOXIA',
                         'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
                         'HALLMARK_ADIPOGENESIS',
                         'HALLMARK_XENOBIOTIC_METABOLISM',
                         'HALLMARK_FATTY_ACID_METABOLISM',
                         'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                         'HALLMARK_GLYCOLYSIS',
                         'HALLMARK_HEME_METABOLISM',
                         'HALLMARK_BILE_ACID_METABOLISM'),]

gene=read.csv('genes.csv')
gene=gene$x
data=read.table('merge.normalzie.txt',header = T,check.names = F,row.names = 1)
rt=data[gene,]

a=grep('',colnames(rt))
rt_MS=rt[,a]
anno=data.frame(row.names = colnames(rt_MS),group=group_list3)

metabolism_MS=metabolism[,rownames(anno)]
metabolism_MS=t(metabolism_MS)

nc =cbind(metabolism_MS,t(rt_MS))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(metabolism_MS),(ncol(nc)-length(gene)+1):ncol(nc)]
rownames(m)=stringr::str_remove(rownames(m),pattern = 'HALLMARK_')
p = rcorr(nc)$P[1:ncol(metabolism_MS),(ncol(nc)-length(gene)+1):ncol(nc)]
library(dplyr)
tmp = matrix(case_when(p<0.001~'***',
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))


library(pheatmap)
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#20854E", "white", "#E18727"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)

