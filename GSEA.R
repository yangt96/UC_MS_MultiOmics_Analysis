library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

gene=""     
expFile=""        
gmtFile=""    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

aL=data[,data[gene,]<median(data[gene,]),drop=F]     #?ͱaH=data[,data[gene,]>=median(data[gene,]),drop=F]    #?߱nL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=meanH-meanL
#???FC=sort(logFC, decreasing=T)
genes=names(logFC)

#??�=read.gmt(gmtFile)

#GSEGSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
#???ab=kkTab[kkTab$pvalue<0.05,]
#kkTte.table(kkTab,file="GSEAp="\t",quote=F,row.names = F)

#???Num=5     #??p=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in high ", gene, " group"))
	pdf(file="GSEA.highExp.pdf", width=6.5, height=5.5)
	print(gseaplot)
	dev.off()
}

#???Num=5     #??own=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title=paste0("Enriched in low ", gene, " group"))
	pdf(file="GSEA.lowExp.pdf", width=6.5, height=5.5)
	print(gseaplot)
	dev.off()
}



####