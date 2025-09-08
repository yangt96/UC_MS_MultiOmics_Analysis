
library(rms)
library(rmda)

inputFile=""    
geneFile=""       

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
row.names(data)=gsub("-", "_", row.names(data))

geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]n row.names(data)){
	data[i,]=ifelse(as.numeric(data[i,])>median(as.numeric(data[i,])), "High", "Low")
}

#??ȡ??�ata)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)

#???ݴ??atadist(rt)
options(datadist="ddist")

#????ģ?mogram(lrmModel, fun=plogis,
	fun.at=c(0.001,0.1,0.5,0.9,0.99),
	lp=F, funlabel="Risk of Disease")
#???????mo.pdf", width=8, height=6)
plot(nomo)
dev.off()

#????У?ibration.pdf", width=5.5, height=5.5)
plot(cali,
	xlab="Predicted probability",
	ylab="Actual probability", sub=F)
dev.off()


######rt_curve(rtStype ~ EPB41L3 + TPST2 + NR3C2 + PDIA6 + KIAA1715 + AP3D1, 
                     data = rt,
                     family = binomial(link = 'logit'),
                     thresholds = seq(0, 1, by = 0.01),
                     confidence.intervals = 0.95)
pdf(file="DCA.pdf", width=6, height=6)
plot_decision_curve(dc,
                    curve.names="nomogram",
                    xlab="Threshold probability",
                    cost.benefit.axis=TRUE,
                    col="#561215",
                    confidence.intervals=FALSE,
                    standa ize=FALSE)

dev.off()Video so?