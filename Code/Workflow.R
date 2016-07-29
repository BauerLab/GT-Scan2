library(MASS)
library(ROCR)
library(caTools)
library(ggplot2)

set.seed(11)

#------------------------------------
#Functions
#------------------------------------

### Calculate combined score ###
calcCombined = function(df, pen) {
  df$Combined = 0
  for (i in c(1:nrow(df))) {
    if (is.na(df[i,'WU.CRISPR'])) { df[i,'Combined'] = df[i,'sgRNAscorer.Hg']-pen }
    else { df[i,'Combined'] = (df[i,'sgRNAscorer.Hg']+df[i,'WU.CRISPR'])/2 }
  }
  return(df)
}

### CV testing ###
cvTest = function(df, fold) {
  dat = df
  nTest = round(nrow(dat)/fold,0)
  pred = c(0,0)
  class = c(0)
  
  for (k in (c(1:fold))) {
    dat.test = dat[sample(nrow(dat), nTest),]
    class = c(class, dat.test$Class)
    dat.train = dat[!(row.names(dat) %in% row.names(dat.test)),]
    cv.glm = glm(Class~., data=dat.train)
    for (i in c(1:nrow(dat.test))) {
      predi = predict(cv.glm, dat.test[i,])
      pred = rbind(pred, c(i,predi))
    }
  }
  pred = pred[-1,]
  class = class[-1]
  pred.to.roc = pred[,2]
  pred.rocr = prediction(pred.to.roc, as.factor(class))
  perf.rocr = performance(pred.rocr, measure='auc', x.measure='cutoff')
  AUC = deparse(round(as.numeric(perf.rocr@y.values),3))
  perf.tpr.rocr = performance(pred.rocr, 'tpr','fpr')
  perf.precrec.rocr = performance(pred.rocr, 'prec','rec')
  perf.sens.rocr = performance(pred.rocr, 'tpr','tnr')
  res = list()
  res[[1]] = AUC
  res[[2]] = perf.tpr.rocr
  res[[3]] = perf.precrec.rocr
  res[[4]] = perf.sens.rocr
  names(res) = c('AUC','ROC','PrecRec','SensSpe')
  return(res)
}

### Calculate area under Prec/Rec curve ###
calcAUPRC = function(list) {
  prerec.df = data.frame(Recall=list[[3]]@x.values[[1]], Precision=list[[3]]@y.values[[1]])
  prerec.df[is.na(prerec.df)] <- 0
  auprc = trapz(prerec.df$Recall, prerec.df$Precision)
  auprc = round(auprc,3)
  return(auprc)
}

### Find optimal cut-off threshold ###
findThr = function(list) {
  thrRes = list
  thr = round(thrRes@alpha.values[[1]][which.max(thrRes@x.values[[1]]+thrRes@y.values[[1]])],3)
  tpr = round(thrRes@y.values[[1]][which.max(thrRes@x.values[[1]]+thrRes@y.values[[1]])],3)
  tnr = round(thrRes@x.values[[1]][which.max(thrRes@x.values[[1]]+thrRes@y.values[[1]])],3)
  return(c(thr,tpr,tnr))
}

### Plot sensitivity vs specificity ###
plotSensSpec = function(cvRes, modelName, thr) {
  par(mar=c(5,4,4,4)+0.3)
  plot(cvRes[[4]]@alpha.values[[1]], cvRes[[4]]@y.values[[1]], col='red', lwd=3, type='l', xlab='', ylab='',
       main=paste('Sensitivity/Specificity of ', modelName), ylim=c(0,1))
  par(new=TRUE)
  plot(cvRes[[4]]@alpha.values[[1]], cvRes[[4]]@x.values[[1]], col='blue', lwd=3,
       type='l', bty='n', xlab='', ylab='', ylim=c(0,1), axes=FALSE)
  par(new=TRUE)
  plot(cvRes[[4]]@alpha.values[[1]], cvRes[[4]]@x.values[[1]]+cvRes[[4]]@y.values[[1]], col='black', lwd=3,
       type='l', bty='n', xlab='', ylab='', ylim=c(0,2), axes=FALSE)
  axis(4, ylim=c(0,2), col='black')
  mtext('Sensitivity/Specificity', side=2, line=2)
  mtext('Sensitivity + Specificity', side=4, line=2)
  mtext('Cut-off Threshold', side=1, line=2)
  abline(v=thr[1], lty=2)
  legend('bottomright', pt.cex=1, cex=0.5, y.intersp=0.75,
         legend=c(paste('Sensitivity: Opt=',thr[2]),
                  paste('Specificity: Opt=',thr[3]),
                  paste('Sensitivity + Specificity: Thr=',thr[1])),
         col=c('red','blue','black'), lwd=3)
}

### Test model on validation set ###
valTest = function(df, model) {
  pred = c(0,0)
  for (i in c(1:nrow(df))) {
    predi = predict(model, df[i,])
    pred = rbind(pred, c(i,predi))
  }
  pred = pred[-1,]
  pred.to.roc = pred[,2]
  pred.rocr = prediction(pred.to.roc, as.factor(df$Class))
  perf.rocr = performance(pred.rocr, measure='auc', x.measure='cutoff')
  AUC = deparse(round(as.numeric(perf.rocr@y.values),3))
  perf.tpr.rocr = performance(pred.rocr,'tpr','fpr')
  perf.precrec.rocr = performance(pred.rocr,'prec','rec')
  res=list()
  res[[1]]=AUC
  res[[2]]=perf.tpr.rocr
  res[[3]]=perf.precrec.rocr
  names(res)=c('AUC','ROC','PrecRec')
  return(res)
}

### Calc TPR (precision), TNR, Recall, F-score ###
calcStats = function(model, thr, dat) {
  pred = predict(model, dat)
  pred.df = data.frame(row.names=row.names(dat), Actual=dat$Class, Predicted=pred)
  for (i in c(1:nrow(pred.df))) {
    if (pred.df[i,'Predicted'] >= thr[1]) { pred.df[i,'Predicted'] = 1 }
    else { pred.df[i,'Predicted'] = 0 }
  }
  
  p00 = 0 #TN
  p01 = 0 #FP
  p10 = 0 #FN
  p11 = 0 #TP
  
  for (i in c(1:nrow(pred.df))) {
    if (pred.df[i,'Actual'] == 0 && pred.df[i,'Predicted'] == 0) { p00 = p00 + 1 }
    else if (pred.df[i,'Actual'] == 0 && pred.df[i,'Predicted'] == 1) { p01 = p01 + 1 }
    else if (pred.df[i,'Actual'] == 1 && pred.df[i,'Predicted'] == 0) { p10 = p10 + 1 }
    else if (pred.df[i,'Actual'] == 1 && pred.df[i,'Predicted'] == 1) { p11 = p11 + 1 }
  }
  
  prec = p11/(p11+p01)
  rec  = p11/(p11+p10)
  npv  = p00/(p00+p10)
  fscore = 2*((prec*rec)/(prec+rec))
  tnr = p00/(p00+p01)
  
  res = list()
  res[[1]] = prec
  res[[2]] = rec
  res[[3]] = npv
  res[[4]] = fscore
  res[[5]] = tnr
  names(res) = c('Precision','Recall','NPV','F-score','TNR')
  return(res)
}

#------------------------------------
#Workflow
#------------------------------------

### This script is set to run from the git main folder ###

c.dat = read.table(file='Data/Chari_dataset.txt',    sep='\t', header=T, row.names=1) #Training data set
d.dat = read.table(file='Data/Doench_dataset.txt',   sep='\t', header=T, row.names=1)
h.dat = read.table(file='Data/Horlbeck_dataset.txt', sep='\t', header=T, row.names=1)
t.dat = read.table(file='Data/Heigwer_dataset.txt',  sep='\t', header=T, row.names=1)

#------------------------------------
#Remove variables that are 
#only one level from training set
#------------------------------------
c.dat = as.data.frame(Filter(function(x) length(unique(x))>1, c.dat))

#------------------------------------
#Build models
#------------------------------------
#Full model           - fm.glm
#Chromatin only       - ch.glm
#sgRNAscorer.Hg only  - sg.glm
#WU.CRISPR only       - wu.glm
#------------------------------------
fm.glm = glm(Class~. - Actual  - sgRNAscorer - WU.CRISPR, data=c.dat)
ch.glm = glm(Class~. - Actual  - sgRNAscorer - WU.CRISPR - Combined, data=c.dat)
sg.glm = glm(Class~sgRNAscorer, data=c.dat)
wu.glm = glm(Class~WU.CRISPR, data=c.dat)

#------------------------------------
#Reduce features of chromatin models
#------------------------------------
fm.glm = stepAIC(fm.glm, trace=TRUE)
ch.glm = stepAIC(ch.glm, trace=TRUE)

#------------------------------------
#Collect model variables
#------------------------------------
fm.vars = all.vars(as.formula(fm.glm))
ch.vars = all.vars(as.formula(ch.glm))
fm.vars = fm.vars[-1]
ch.vars = ch.vars[-1]

#------------------------------------
#Perform 10-fold CV on models
#------------------------------------
fm.cvRes = cvTest(subset(c.dat, select=c(fm.vars,'Class')),10)
ch.cvRes = cvTest(subset(c.dat, select=c(ch.vars,'Class')),10)
sg.cvRes = cvTest(subset(c.dat, select=c('sgRNAscorer.Hg','Class')),10)
wu.cvRes = cvTest(subset(c.dat, select=c('WU.CRISPR','Class')),10)

#------------------------------------
#Calculate area under Pre/Rec curve
#------------------------------------
fm.auprc = calcAUPRC(fm.cvRes)
ch.auprc = calcAUPRC(ch.cvRes)
sg.auprc = calcAUPRC(sg.cvRes)
wu.auprc = calcAUPRC(wu.cvRes)

#------------------------------------
#Plot ROC and Pre/Rec curves
#------------------------------------
pdf('Results/Model.10fold-CV.ROC.pdf')
plot(fm.cvRes[[2]], col='blue',      lwd=3, main='10-fold CV ROC curves')
plot(ch.cvRes[[2]], col='red',       lwd=3, add=TRUE)
plot(sg.cvRes[[2]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.cvRes[[2]], col='orange',    lwd=3, add=TRUE)
plot(fm.cvRes[[2]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUC=',     fm.cvRes[[1]]),
                paste('Chromatin Only: AUC=', ch.cvRes[[1]]),
                paste('sgRNAscorer: AUC=',    sg.cvRes[[1]]),
                paste('WU.CRISPR: AUC=',      wu.cvRes[[1]])),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75) 
dev.off()

pdf('Results/Model.10fold-CV.PreRec.pdf')
plot(fm.cvRes[[3]], col='blue',      lwd=3, ylim=c(0,1), main='10-fold CV Precision/Recall curves')
plot(ch.cvRes[[3]], col='red',       lwd=3, add=TRUE)
plot(sg.cvRes[[3]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.cvRes[[3]], col='orange',    lwd=3, add=TRUE)
plot(fm.cvRes[[3]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUPRC=',     fm.auprc),
                paste('Chromatin only: AUPRC=', ch.auprc),
                paste('sgRNAscorer: AUPRC=',    sg.auprc),
                paste('WU-CRISPR: AUPRC=',      wu.auprc)),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75)
dev.off()

#------------------------------------
#Validating on Doench dataset
#------------------------------------
fm.d.valRes = valTest(d.dat, fm.glm)
ch.d.valRes = valTest(d.dat, ch.glm)
sg.d.valRes = valTest(d.dat, sg.glm)
wu.d.valRes = valTest(d.dat, wu.glm)

fm.d.auprc = calcAUPRC(fm.d.valRes)
ch.d.auprc = calcAUPRC(ch.d.valRes)
sg.d.auprc = calcAUPRC(sg.d.valRes)
wu.d.auprc = calcAUPRC(wu.d.valRes)

pdf('Results/Model.Doench_validation.ROC.pdf')
plot(fm.d.valRes[[2]], col='blue',      lwd=3, main='Doench validation ROC ')
plot(ch.d.valRes[[2]], col='red',       lwd=3, add=TRUE)
plot(sg.d.valRes[[2]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.d.valRes[[2]], col='orange',    lwd=3, add=TRUE)
plot(fm.d.valRes[[2]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUC=',     fm.d.valRes[[1]]),
                paste('Chromatin Only: AUC=', ch.d.valRes[[1]]),
                paste('sgRNAscorer: AUC=',    sg.d.valRes[[1]]),
                paste('WU.CRISPR: AUC=',      wu.d.valRes[[1]])),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75)
dev.off()

pdf('Results/Model.Doench_validation.PreRec.pdf')
plot(fm.d.valRes[[3]], col='blue',      lwd=3, ylim=c(0,1), main='Doench validation Precision/Recall curves')
plot(ch.d.valRes[[3]], col='red',       lwd=3, add=TRUE)
plot(sg.d.valRes[[3]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.d.valRes[[3]], col='orange',    lwd=3, add=TRUE)
plot(fm.d.valRes[[3]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUPRC=',     fm.d.auprc),
                paste('Chromatin only: AUPRC=', ch.d.auprc),
                paste('sgRNAscorer: AUPRC=',    sg.d.auprc),
                paste('WU-CRISPR: AUPRC=',      wu.d.auprc)),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75)
dev.off()

#------------------------------------
#Validating on Horlbeck dataset
#------------------------------------
fm.h.valRes = valTest(h.dat, fm.glm)
ch.h.valRes = valTest(h.dat, ch.glm)
sg.h.valRes = valTest(h.dat, sg.glm)
wu.h.valRes = valTest(h.dat, wu.glm)

fm.h.auprc = calcAUPRC(fm.h.valRes)
ch.h.auprc = calcAUPRC(ch.h.valRes)
sg.h.auprc = calcAUPRC(sg.h.valRes)
wu.h.auprc = calcAUPRC(wu.h.valRes)

pdf('Results/Model.Horlbeck_validation.ROC.pdf')
plot(fm.h.valRes[[2]], col='blue',      lwd=3, main='Horlbeck validation ROC ')
plot(ch.h.valRes[[2]], col='red',       lwd=3, add=TRUE)
plot(sg.h.valRes[[2]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.h.valRes[[2]], col='orange',    lwd=3, add=TRUE)
plot(fm.h.valRes[[2]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUC=',     fm.h.valRes[[1]]),
                paste('Chromatin Only: AUC=', ch.h.valRes[[1]]),
                paste('sgRNAscorer: AUC=',    sg.h.valRes[[1]]),
                paste('WU.CRISPR: AUC=',      wu.h.valRes[[1]])),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75)
dev.off()

pdf('Results/Model.Horlbeck_validation.PreRec.pdf')
plot(fm.h.valRes[[3]], col='blue',      lwd=3, ylim=c(0,1), main='Horlbeck validation Precision/Recall curves')
plot(ch.h.valRes[[3]], col='red',       lwd=3, add=TRUE)
plot(sg.h.valRes[[3]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.h.valRes[[3]], col='orange',    lwd=3, add=TRUE)
plot(fm.h.valRes[[3]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUPRC=',     fm.h.auprc),
                paste('Chromatin only: AUPRC=', ch.h.auprc),
                paste('sgRNAscorer: AUPRC=',    sg.h.auprc),
                paste('WU-CRISPR: AUPRC=',      wu.h.auprc)),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75)
dev.off()

#------------------------------------
#Validating on Heigwer dataset
#------------------------------------
fm.t.valRes = valTest(t.dat, fm.glm)
ch.t.valRes = valTest(t.dat, ch.glm)
sg.t.valRes = valTest(t.dat, sg.glm)
wu.t.valRes = valTest(t.dat, wu.glm)

fm.t.auprc = calcAUPRC(fm.t.valRes)
ch.t.auprc = calcAUPRC(ch.t.valRes)
sg.t.auprc = calcAUPRC(sg.t.valRes)
wu.t.auprc = calcAUPRC(wu.t.valRes)

pdf('Results/Model.Heigwer_validation.ROC.pdf')
plot(fm.t.valRes[[2]], col='blue',      lwd=3, main='Heigwer validation ROC ')
plot(ch.t.valRes[[2]], col='red',       lwd=3, add=TRUE)
plot(sg.t.valRes[[2]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.t.valRes[[2]], col='orange',    lwd=3, add=TRUE)
plot(fm.t.valRes[[2]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUC=',     fm.t.valRes[[1]]),
                paste('Chromatin Only: AUC=', ch.t.valRes[[1]]),
                paste('sgRNAscorer: AUC=',    sg.t.valRes[[1]]),
                paste('WU.CRISPR: AUC=',      wu.t.valRes[[1]])),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75)
dev.off()

pdf('Results/Model.Heigwer_validation.PreRec.pdf')
plot(fm.t.valRes[[3]], col='blue',      lwd=3, ylim=c(0,1), main='Heigwer validation Precision/Recall curves')
plot(ch.t.valRes[[3]], col='red',       lwd=3, add=TRUE)
plot(sg.t.valRes[[3]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.t.valRes[[3]], col='orange',    lwd=3, add=TRUE)
plot(fm.t.valRes[[3]], col='blue',      lwd=3, add=TRUE) #Plot again to bring to front
legend('bottomright', title=expression(bold('Model')), pt.cex=1, cex=0.5,
       legend=c(paste('Full Model: AUPRC=',     fm.t.auprc),
                paste('Chromatin only: AUPRC=', ch.t.auprc),
                paste('sgRNAscorer: AUPRC=',    sg.t.auprc),
                paste('WU-CRISPR: AUPRC=',      wu.t.auprc)),
       col=c('blue','red','darkgreen','orange'),
       lwd=3, bty='n', y.intersp=0.75)
dev.off()

#------------------------------------
#Find optimal cut-off threshold
#------------------------------------
fm.thr = findThr(fm.cvRes[[4]])
ch.thr = findThr(ch.cvRes[[4]])
sg.thr = findThr(sg.cvRes[[4]])
wu.thr = findThr(wu.cvRes[[4]])

#------------------------------------
#Plot Sensitivity vs Specificity
#------------------------------------
pdf('./Results/SensSpec.FullModel.pdf')
plotSensSpec(fm.cvRes, 'Full Model',     fm.thr)
dev.off()

pdf('./Results/SensSpec.ChromatinOnly.pdf')
plotSensSpec(ch.cvRes, 'Chromatin Only', ch.thr)
dev.off()

pdf('./Results/SensSpec.sgRNAscorer.pdf')
plotSensSpec(sg.cvRes, 'sgRNAscorer',    sg.thr)
dev.off()

pdf('./Results/SensSpec.WU-CRISPR.pdf')
plotSensSpec(wu.cvRes, 'WU-CRISPR',      wu.thr)
dev.off()

#------------------------------------
#Calc stats on validation datasets
#------------------------------------
### Doench dataset ###
fm.d.stats = calcStats(fm.glm, fm.thr, d.dat)
ch.d.stats = calcStats(ch.glm, ch.thr, d.dat)
sg.d.stats = calcStats(sg.glm, sg.thr, d.dat)
wu.d.stats = calcStats(wu.glm, wu.thr, d.dat)

### Horlbeck dataset ###
fm.h.stats = calcStats(fm.glm, fm.thr, h.dat)
ch.h.stats = calcStats(ch.glm, ch.thr, h.dat)
sg.h.stats = calcStats(sg.glm, sg.thr, h.dat)
wu.h.stats = calcStats(wu.glm, wu.thr, h.dat)

### Heigwer dataset ###
fm.t.stats = calcStats(fm.glm, fm.thr, t.dat)
ch.t.stats = calcStats(ch.glm, ch.thr, t.dat)
sg.t.stats = calcStats(sg.glm, sg.thr, t.dat)
wu.t.stats = calcStats(wu.glm, wu.thr, t.dat)

#------------------------------------
#Plot average of stats across validations sets
#------------------------------------
prec.df = data.frame(Precision=c(fm.d.stats$Precision,fm.h.stats$Precision,fm.t.stats$Precision,
                                 ch.d.stats$Precision,ch.h.stats$Precision,sg.t.stats$Precision,
                                 sg.d.stats$Precision,sg.h.stats$Precision,sg.t.stats$Precision,
                                 wu.d.stats$Precision,wu.h.stats$Precision,wu.t.stats$Precision),
                     Model=c('Full Model','Full Model','Full Model',
                             'Chromatin only','Chromatin only','Chromatin only',
                             'sgRNAscorer','sgRNAscorer','sgRNAscorer',
                             'WU-CRISPR','WU-CRISPR','WU-CRISPR'))

recc.df = data.frame(Recall=c(fm.d.stats$Recall,fm.h.stats$Recall,fm.t.stats$Recall,
                              ch.d.stats$Recall,ch.h.stats$Recall,sg.t.stats$Recall,
                              sg.d.stats$Recall,sg.h.stats$Recall,sg.t.stats$Recall,
                              wu.d.stats$Recall,wu.h.stats$Recall,wu.t.stats$Recall),
                     Model=c('Full Model','Full Model','Full Model',
                             'Chromatin only','Chromatin only','Chromatin only',
                             'sgRNAscorer','sgRNAscorer','sgRNAscorer',
                             'WU-CRISPR','WU-CRISPR','WU-CRISPR'))

npv.df = data.frame(NPV=c(fm.d.stats$NPV,fm.h.stats$NPV,fm.t.stats$NPV,
                          ch.d.stats$NPV,ch.h.stats$NPV,sg.t.stats$NPV,
                          sg.d.stats$NPV,sg.h.stats$NPV,sg.t.stats$NPV,
                          wu.d.stats$NPV,wu.h.stats$NPV,wu.t.stats$NPV),
                    Model=c('Full Model','Full Model','Full Model',
                            'Chromatin only','Chromatin only','Chromatin only',
                            'sgRNAscorer','sgRNAscorer','sgRNAscorer',
                            'WU-CRISPR','WU-CRISPR','WU-CRISPR'))

fscore.df = data.frame(F.score=c(fm.d.stats$`F-score`,fm.h.stats$`F-score`,fm.t.stats$`F-score`,
                                 ch.d.stats$`F-score`,ch.h.stats$`F-score`,sg.t.stats$`F-score`,
                                 sg.d.stats$`F-score`,sg.h.stats$`F-score`,sg.t.stats$`F-score`,
                                 wu.d.stats$`F-score`,wu.h.stats$`F-score`,wu.t.stats$`F-score`),
                       Model=c('Full Model','Full Model','Full Model',
                               'Chromatin only','Chromatin only','Chromatin only',
                               'sgRNAscorer','sgRNAscorer','sgRNAscorer',
                               'WU-CRISPR','WU-CRISPR','WU-CRISPR'))

tnr.df = data.frame(TNR=c(fm.d.stats$TNR,fm.h.stats$TNR,fm.t.stats$TNR,
                          ch.d.stats$TNR,ch.h.stats$TNR,sg.t.stats$TNR,
                          sg.d.stats$TNR,sg.h.stats$TNR,sg.t.stats$TNR,
                          wu.d.stats$TNR,wu.h.stats$TNR,wu.t.stats$TNR),
                    Model=c('Full Model','Full Model','Full Model',
                            'Chromatin only','Chromatin only','Chromatin only',
                            'sgRNAscorer','sgRNAscorer','sgRNAscorer',
                            'WU-CRISPR','WU-CRISPR','WU-CRISPR'))

prec.df$Model   <- factor(prec.df$Model,   levels=c('Full Model','sgRNAscorer','WU-CRISPR','Chromatin only'))
recc.df$Model   <- factor(recc.df$Model,   levels=c('Full Model','sgRNAscorer','WU-CRISPR','Chromatin only'))
npv.df$Model    <- factor(tnr.df$Model,    levels=c('Full Model','sgRNAscorer','WU-CRISPR','Chromatin only'))
fscore.df$Model <- factor(fscore.df$Model, levels=c('Full Model','sgRNAscorer','WU-CRISPR','Chromatin only'))
tnr.df$Model    <- factor(tnr.df$Model,    levels=c('Full Model','sgRNAscorer','WU-CRISPR','Chromatin only'))

prec.plot = ggplot(prec.df, aes(x=Model, y=Precision, fill='red')) + geom_boxplot() +
  labs(title='Average Precision across validation sets') +
  scale_y_continuous(limits=c(0,1)) +
  guides(fill=FALSE) +
  geom_jitter(width=0.3)

recc.plot = ggplot(recc.df, aes(x=Model, y=Recall, fill='red')) + geom_boxplot() +
  labs(title='Average Recall across validation sets') +
  scale_y_continuous(limits=c(0,1)) +
  guides(fill=FALSE) +
  geom_jitter(width=0.3)

npv.plot = ggplot(npv.df, aes(x=Model, y=NPV, fill='red')) + geom_boxplot() +
  labs(title='Average NPV across validation sets') +
  scale_y_continuous(limits=c(0,1)) +
  guides(fill=FALSE) +
  geom_jitter(width=0.3)

fscore.plot = ggplot(fscore.df, aes(x=Model, y=F.score, fill='red')) + geom_boxplot() +
  labs(title='Average F-score across validation sets') +
  scale_y_continuous(limits=c(0,1)) +
  guides(fill=FALSE) +
  geom_jitter(width=0.3)

tnr.plot = ggplot(tnr.df, aes(x=Model, y=TNR, fill='red')) + geom_boxplot() +
  labs(title='Average TNR across validation sets') +
  scale_y_continuous(limits=c(0,1)) +
  guides(fill=FALSE) +
  geom_jitter(width=0.3)

pdf('./Results/Stats.precision.pdf')
prec.plot
dev.off()

pdf('./Results/Stats.recall.pdf')
recc.plot
dev.off()

pdf('./Results/Stats.NPV.pdf')
npv.plot
dev.off()

pdf('./Results/Stats.Fscore.pdf')
fscore.plot
dev.off()

pdf('./Results/Stats.Spec.pdf')
tnr.plot
dev.off()

