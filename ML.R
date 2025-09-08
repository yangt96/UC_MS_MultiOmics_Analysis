load()
load()

rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)


up=read.table('up.txt',header = F)
down=read.table('down.txt',header = F)
diff_gene=c(up$V1,down$V1)


diff_gene=diff_gene[diff_gene %in% rownames(rt)]

setwd('d:/')

data=exprSet[diff_gene,]
data=as.data.frame(t(data))
group_list5

levels(group_list5)[levels(group_list5) == ""] <- ""

library(randomForest)
set.seed(245083)
rf=randomForest(x = data,y = group,ntree = 1000)

plot(rf, main="Random forest", lwd=2)

optionTrees=which.min(rf$err.rate[,1])
optionTrees
set.seed(245083)
rf2=randomForest(x = data,y = group, ntree=optionTrees)
importance=importance(x=rf2)
dev.off()
varImpPlot(rf2, main="")
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0])    
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)
dev.off()
library(dplyr)
library(glmnet)
x=as.matrix(data[,rfGenes])
y=unlist(group)
set.seed(12345)
fit <- glmnet(x,y,alpha=1,family="binomial")
plot(fit,xvar="lambda",label=F)
set.seed(123456)
cvfit <- cv.glmnet(x,y,alpha=1,family="binomial",nfolds = 10)
plot(cvfit)
coef =coef(fit,s = cvfit$lambda.min)
index = which(coef !=0)
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
geneCoef  
write.csv(geneCoef,file ='genecoef.csv',quote = F,row.names = F)
genes=geneCoef[,1][2:length(lassoGene)]
genes
save(genes,file ='rf_lasso_genes.Rdata')
load('rf_lasso_genes.Rdata')
library(pROC)
for (i in genes) {
  pdf(paste0(i,'.pdf'),width = 4,height = 4)
  plot.roc(y,as.numeric(data[,i]),print.auc=T)
  dev.off()
}

rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)

load('rf_lasso_genes.Rdata')
rt=rt[genes,]

a=grep('',colnames(rt))
X_train=t(rt[,a])
load('')
Y_train=ifelse(group_list5 =='HC',0,1)

b=grep('',colnames(rt))
X_test=t(rt[,b])
load('')
Y_test=ifelse(group_list1 == 'control',0,1)

library(xgboost)
library(rBayesianOptimization)

dtrain <- xgb.DMatrix(data = X_train, label = Y_train)
cv_folds <- KFold(Y_train, nfolds = 4, stratified = TRUE, seed = 0)
xgb_cv_bayes <- function(eta, max.depth, min_child_weight, subsample) {
  cv <- xgb.cv(
    params = list(
      booster = "gbtree",
      eta = eta,
      max_depth = max.depth,
      min_child_weight = min_child_weight,
      subsample = subsample,
      colsample_bytree = 0.6,
      lambda = 1,
      alpha = 0,
      objective = "binary:logistic",
      eval_metric = "auc"
    ),
    data = dtrain,
    nround = 30,
    folds = cv_folds,
    prediction = TRUE,
    showsd = TRUE,
    early.stop.round = 5,
    maximize = TRUE,
    verbose = 0
  )
  list(
    Score = cv$evaluation_log[, max(test_auc_mean)], Pred = cv$pred
  )
}
OPT_Res <- BayesianOptimization(
  xgb_cv_bayes,
  bounds = list(
    eta = c(0.01L, 0.05L, 0.1L, 0.3L),
    max.depth = c(6L, 8L, 12L),
    min_child_weight = c(1L, 10L),
    subsample = c(0.5, 0.8, 1)),
  init_grid_dt = NULL,
  init_points = 10,
  n_iter = 10,
  acq = "ucb",
  kappa = 2.576,
  eps = 0.0,
  verbose = TRUE
)
params <- list(
  "eta" = unname(OPT_Res$Best_Par["eta"]),
  "max_depth" = unname(OPT_Res$Best_Par["max.depth"]),
  "colsample_bytree" = 1,
  "min_child_weight" = unname(OPT_Res$Best_Par["min_child_weight"]),
  "subsample"= unname(OPT_Res$Best_Par["subsample"]),
  "objective"="binary:logistic",
  "gamma" = 1,
  "lambda" = 1,
  "alpha" = 0,
  "max_delta_step" = 0,
  "colsample_bylevel" = 1,
  "eval_metric"= "auc",
  "set.seed" = 176
)


watchlist <- list("train" = dtrain)
nround = 30
xgb.model <- xgb.train(params, dtrain, nround, watchlist)
dtest <- xgb.DMatrix(data = X_test)

XGB_train_Predictions <- predict(object = xgb.model, newdata = dtrain, type = 'prob')

XGB_Predictions <- predict(object = xgb.model, newdata = dtest, type = 'prob')


pROC::plot.roc(as.factor(Y_train),XGB_train_Predictions,print.auc=T)
pROC::plot.roc(as.factor(Y_test),XGB_Predictions,print.auc=T)


XGB_train_Predictions_bi <- ifelse(XGB_train_Predictions > 0.5, '1', '0')
XGB_Predictions_bi <- ifelse(XGB_Predictions > 0.5, '1', '0')
library(caret)
con_XGB_train <- confusionMatrix(as.factor(XGB_train_Predictions_bi), as.factor(Y_train))
con_XGB_test <- confusionMatrix(as.factor(XGB_Predictions_bi),as.factor(Y_test))
library(tidyverse)
XGB_pred_outcome <- cbind(as.numeric(XGB_train_Predictions), as.numeric(Y_train)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))

XGB_fg <- XGB_pred_outcome %>% 
  filter(outcome == 1) %>%
  pull(predictions)

XGB_bg <- XGB_pred_outcome %>% 
  filter(outcome == 0) %>%
  pull(predictions)

XGB_pr <- PRROC::pr.curve(
  scores.class0 = XGB_fg, scores.class1 = XGB_bg, curve = TRUE)


theme_bluewhite <- function (base_size = 11, base_family = "serif") {
  theme_bw() %+replace% 
    theme(
      text = element_text(family = "serif"),
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = "grey97"),
      panel.border = element_rect(color = "#bc3c29", fill = NA, size = 1), ##05014a
      axis.line = element_line(color = "grey97"),
      axis.ticks = element_line(color = "grey25"),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "grey25", size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 15, hjust = 0.5),
      strip.background = element_rect(fill = '#0072b5'),
      strip.text = element_text(size = 10, colour = 'white'), # changes the facet wrap text size and color
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}
XGB_pr_curve <- XGB_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "XGBoost")
library(ggplot2)
XGB_prPlot <-
  ggplot(XGB_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("XGBoost")) +
  labs(title= paste0("Precision-Recall Curves, AUC=",round(XGB_pr$auc.integral,3)), 
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite() 

XGB_prPlot
library(pROC)

XGB_roc <- roc(as.numeric(Y_train),as.numeric(XGB_train_Predictions))

XGB_rocPlot <- 
  ggroc(XGB_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = "Light GBM") +
  labs(title= paste0("ROC Curves, AUC=",round(XGB_roc$auc,3)), 
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() 

XGB_rocPlot
library(ggplot2)
library(patchwork)
(XGB_rocPlot + XGB_prPlot) + plot_annotation(title = "Performance in training set", tag_levels = "A")

XGB_pred_outcome <- cbind(as.numeric(XGB_Predictions), as.numeric(Y_test)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))
XGB_fg <- XGB_pred_outcome %>% 
  filter(outcome == 1) %>%
  pull(predictions)

XGB_bg <- XGB_pred_outcome %>% 
  filter(outcome == 0) %>%
  pull(predictions)

XGB_pr <- PRROC::pr.curve(
  scores.class0 = XGB_fg, scores.class1 = XGB_bg, curve = TRUE)
theme_bluewhite <- function (base_size = 11, base_family = "serif") {
  theme_bw() %+replace% 
    theme(
      text = element_text(family = "serif"),
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = "grey97"),
      panel.border = element_rect(color = "#bc3c29", fill = NA, size = 1), ##05014a
      axis.line = element_line(color = "grey97"),
      axis.ticks = element_line(color = "grey25"),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "grey25", size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 15, hjust = 0.5),
      strip.background = element_rect(fill = '#0072b5'),
      strip.text = element_text(size = 10, colour = 'white'), # changes the facet wrap text size and color
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

XGB_pr_curve <- XGB_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "XGBoost")

library(ggplot2)
XGB_prPlot <-
  ggplot(XGB_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("XGBoost")) +
  labs(title= paste0("Precision-Recall Curves, AUC=",round(XGB_pr$auc.integral,3)), 
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite()

XGB_prPlot
library(pROC)

XGB_roc <- roc(as.numeric(Y_test),as.numeric(XGB_Predictions))

XGB_rocPlot <- 
  ggroc(XGB_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = "XGBoost") +
  labs(title= paste0("ROC Curves, AUC=",round(XGB_roc$auc,3)), 
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite()
XGB_rocPlot

(XGB_rocPlot + XGB_prPlot) + plot_annotation(title = "Performance in validation set", tag_levels = "A")

genes
write.csv(genes,file ='genes.csv',quote = F)


load('MS.Rdata')

rt=read.table('merge.preNorm.txt', header=T, sep="\t", check.names=F, row.names=1)

load('rf_lasso_genes.Rdata')

rt=rt[genes,]

a=grep('',colnames(rt))
X_train=t(rt[,a])
Y_train=ifelse(group_list3 =='',0,1)

library(xgboost)
library(rBayesianOptimization)

dtrain <- xgb.DMatrix(data = X_train, label = Y_train)
cv_folds <- KFold(Y_train, nfolds = 10, stratified = TRUE, seed = 0)
xgb_cv_bayes <- function(eta, max.depth, min_child_weight, subsample) {
  cv <- xgb.cv(
    params = list(
      booster = "gbtree",
      eta = eta,
      max_depth = max.depth,
      min_child_weight = min_child_weight,
      subsample = subsample,
      colsample_bytree = 0.6,
      lambda = 1,
      alpha = 0,
      objective = "binary:logistic",
      eval_metric = "auc"
    ),
    data = dtrain,
    nround = 10,
    folds = cv_folds,
    prediction = TRUE,
    showsd = TRUE,
    early.stop.round = 5,
    maximize = TRUE,
    verbose = 0
  )
  list(
    Score = cv$evaluation_log[, max(test_auc_mean)], Pred = cv$pred
  )
}

OPT_Res <- BayesianOptimization(
  xgb_cv_bayes,
  bounds = list(
    eta = c(0.01L, 0.05L, 0.1L, 0.3L),
    max.depth = c(6L, 8L, 12L),
    min_child_weight = c(1L, 10L),
    subsample = c(0.5, 0.8, 1)),
  init_grid_dt = NULL,
  init_points = 10,
  n_iter = 10,
  acq = "ucb",
  kappa = 2.576,
  eps = 0.0,
  verbose = TRUE
)
params <- list(
  "eta" = unname(OPT_Res$Best_Par["eta"]),
  "max_depth" = unname(OPT_Res$Best_Par["max.depth"]),
  "colsample_bytree" = 1,
  "min_child_weight" = unname(OPT_Res$Best_Par["min_child_weight"]),
  "subsample"= unname(OPT_Res$Best_Par["subsample"]),
  "objective"="binary:logistic",
  "gamma" = 1,
  "lambda" = 1,
  "alpha" = 0,
  "max_delta_step" = 0,
  "colsample_bylevel" = 1,
  "eval_metric"= "auc",
  "set.seed" = 176
)

watchlist <- list("train" = dtrain)
nround = 20
xgb.model <- xgb.train(params, dtrain, nround, watchlist)


XGB_train_Predictions <- predict(object = xgb.model, newdata = dtrain, type = 'prob')

pROC::plot.roc(as.factor(Y_train),XGB_train_Predictions,print.auc=T)

XGB_train_Predictions_bi <- ifelse(XGB_train_Predictions > 0.5, '1', '0')


library(caret)

con_XGB_train <- confusionMatrix(as.factor(XGB_train_Predictions_bi), as.factor(Y_train))
library(tidyverse)
XGB_pred_outcome <- cbind(as.numeric(XGB_train_Predictions), as.numeric(Y_train)) %>% 
  data.frame() %>% 
  setNames(c("predictions", "outcome"))

XGB_fg <- XGB_pred_outcome %>% 
  filter(outcome == 1) %>%
  pull(predictions)

XGB_bg <- XGB_pred_outcome %>% 
  filter(outcome == 0) %>%
  pull(predictions)

XGB_pr <- PRROC::pr.curve(
  scores.class0 = XGB_fg, scores.class1 = XGB_bg, curve = TRUE)

theme_bluewhite <- function (base_size = 11, base_family = "serif") {
  theme_bw() %+replace% 
    theme(
      text = element_text(family = "serif"),
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = "grey97"),
      panel.border = element_rect(color = "#bc3c29", fill = NA, size = 1), ##05014a
      axis.line = element_line(color = "grey97"),
      axis.ticks = element_line(color = "grey25"),
      axis.title = element_text(size = 10),
      axis.text = element_text(color = "grey25", size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 15, hjust = 0.5),
      strip.background = element_rect(fill = '#0072b5'),
      strip.text = element_text(size = 10, colour = 'white'), # changes the facet wrap text size and color
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

XGB_pr_curve <- XGB_pr$curve %>% 
  data.frame() %>% 
  mutate(X4 = "XGBoost")

library(ggplot2)
XGB_prPlot <-
  ggplot(XGB_pr_curve,aes(x = X1, y = X2, color = X4)) +
  geom_line() +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = c("XGBoost")) +
  labs(title= paste0("Precision-Recall Curves, AUC=",round(XGB_pr$auc.integral,3)), 
       y = "Precision",
       x = "Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) +
  theme_bluewhite() 

XGB_prPlot
library(pROC)

XGB_roc <- roc(as.numeric(Y_train),as.numeric(XGB_train_Predictions))

XGB_rocPlot <- 
  ggroc(XGB_roc,legacy.axes = TRUE, linetype = 1) +
  geom_abline(show.legend = TRUE, alpha = 0.3) +
  scale_color_viridis_d(option = "D", name = "Model",
                        labels = "Light GBM") +
  labs(title= paste0("ROC Curves, AUC=",round(XGB_roc$auc,3)), 
       y = "Sensitivity",
       x = "1-Specificity") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme_bluewhite() 

XGB_rocPlot

library(ggplot2)
library(patchwork)
(XGB_rocPlot + XGB_prPlot) + plot_annotation(title = "Performance in MS cohort", tag_levels = "A")

