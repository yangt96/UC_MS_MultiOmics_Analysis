write.table(data.frame(ID=immune), file="CIBERSORT-Resultss", sep="\t", quote=FALSE, row.names=FALSE)

data = data[c('EPB41L3','TPST2','NR3C2','PDIA6','GMCL1','LNPK','AP3D1'),]
data = t(data)

Control = read.table("S2.txt", header=FALSE, sep="\t", check.names=FALSE)

Treat = read.table("S1.txt", header=FALSE, sep="\t", check.names=FALSE)

data = data[c(Control[,1], Treat[,1]),]

conNum = length(rownames(Control))

treatNum = length(rownames(Treat))
group = c(rep("Control", conNum), rep("Treat", treatNum))

group = c(rep("Control", conNum), rep("Treat", treatNum))


rt = cbind(as.data.frame(data), Type=group)


ddist = datadist(rt)
options(datadist="ddist")

paste(colnames(data), collapse="+")

IrmModel = lrm(Type ~ EPB41L3+TPST2+NR3C2+PDIA6+GMCL1+KIAA1715+AP3D1, data=rt, x=T, y=T,maxit=1000)

nomo = nomogram(IrmModel, fun=plogis,
                fun.at=c(0.1,0.5,0.99),
                lp=FALSE, funlabel="Risk of Disease")


pdf("Nom.pdf", width=11, height=6)
plot(nomo)
dev.off()


cali <- calibrate(IrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=6, height=6)
plot(cali, 
     xlab="Predicted probability",
     ylab="Actual probability", sub=FALSE)
dev.off()
rt$rtStype <- ifelse(rt$Type == "HC", 0, 1)

dc <- decision_curve(rtStype ~ EPB41L3 + TPST2 + NR3C2 + PDIA6 + KIAA1715 + AP3D1, 
                     data = rt,
                     family = binomial(link = 'logit'),
                     thresholds = seq(0, 1, by = 0.01),
                     confidence.intervals = 0.95)

rtStype <- ifelse(rt$Type == "HC", 0, 1)
dc <- decision_curve(rtStype ~ EPB41L3+TPST2+NR3C2+PDIA6+KIAA1715+AP3D1, data=rt,
                     family = binomial(link = 'logit'),
                     thresholds= seq(0, 1, by = 0.01),
                     confidence.intervals = 0.95)

pdf(file="DCA.pdf", width=6, height=6)
plot_decision_curve(dc,
                    curve.names="nomogram",
                    xlab="Threshold probability",
                    cost.benefit.axis=TRUE,
                    col="#561215",
                    confidence.intervals=FALSE,
                    standardize=FALSE)

dev.off()
outcome <- as.numeric(outcome)

rt$Type = factor(rt$Type, levels = c("", ""))


rt$rtStype <- ifelse(rt$Type == "", 0, 1)


dc <- decision_curve(rtStype ~ EPB41L3+TPST2+NR3C2+PDIA6+GMCL1+KIAA1715+AP3D1, data=rt,
                     family = binomial(link = 'logit'),
                     thresholds= seq(0, 1, by = 0.01),
                     confidence.intervals = 0.95)


