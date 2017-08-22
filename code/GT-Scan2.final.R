library(randomForest)
library(ggplot2)
library(ROCR)
library(gplots)
library(reshape)
library(MASS)
library(caTools)

set.seed(18) # Seed is set so that the results are always reproducible

setwd('~/Desktop/GTscan_repository/')

#------------------------------------
# Functions
#------------------------------------

### Construct color vectors for heatmaps
### Colored based on class
color_map = function(classification) {
  
  if (classification=="1")
    "#FF0000" # Red
  else
    "#0000FF" # Blue
  
} #end color.map = function(classification)

### Construct color vectors for heatmaps
### Colored based on classifcaiont (TP, FP etc)
color_map2 = function(classification) {
  if (classification=="TP")
    "#FF0000" #Red
  else if (classification=="FN")
    "#FFB90F" #Orange
  else if (classification=="TN")
    "#0000FF" #Blue
  else
    "#6495ED" #Ligth Blue
}

### Calc TPR (precision), TNR, Recall, F-score ###
glm_calcStats = function(model, thr, dat) {
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
  acc = (p11+p00)/(p11+p00+p01+p10)
  
  res = list()
  res[[1]] = prec
  res[[2]] = rec
  res[[3]] = npv
  res[[4]] = fscore
  res[[5]] = tnr
  res[[6]] = acc
  names(res) = c('Precision','Recall','NPV','F-score','TNR','Accuracy')
  return(res)
}

### Calc TPR (precision), TNR, Recall, F-score ###
rf_calcStats = function(model, thr, dat) {
  pred = predict(model, dat, type='prob')
  pred.df = data.frame(row.names=row.names(dat), Actual=dat$Class, Predicted=pred[,2])
  
  for (i in c(1:nrow(pred.df))) {
    if (pred.df[i,'Predicted'] >= thr[1]) { pred.df[i,'Predicted'] = 1 }
    else { pred.df[i,'Predicted'] = 0}
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
  acc = (p11+p00)/(p11+p00+p01+p10)
  
  res = list()
  res[[1]] = prec
  res[[2]] = rec
  res[[3]] = npv
  res[[4]] = fscore
  res[[5]] = tnr
  res[[6]] = acc
  names(res) = c('Precision','Recall','NPV','F-score','TNR','Accuracy')
  return(res)
}

### Find optimal cut-off threshold ###
findThr = function(list) {
  thrRes = list
  thr = round(thrRes@alpha.values[[1]][which.max(thrRes@x.values[[1]]+thrRes@y.values[[1]])],3)
  tpr = round(thrRes@y.values[[1]][which.max(thrRes@x.values[[1]]+thrRes@y.values[[1]])],3)
  tnr = round(thrRes@x.values[[1]][which.max(thrRes@x.values[[1]]+thrRes@y.values[[1]])],3)
  return(c(thr,tpr,tnr))
}

### Calculate combined score
calcCombined = function(df, pen) {
  df$Combined = 0
  for (i in c(1:nrow(df))) {
    if (is.na(df[i,'WU.CRISPR'])) { df[i,'Combined'] = df[i,'sgRNAscorer.Hg']-pen }
    else { df[i,'Combined'] = (df[i,'sgRNAscorer.Hg']+df[i,'WU.CRISPR'])/2 }
  }
  return(df)
}

### Build sequential RFs to find best signature
seq_rf = function(train, imp.ord, reps) {
  dat.t = train
  imp = imp.ord
  genes = row.names(subset(imp, Rank<=1))
  dat.sub = as.data.frame(dat.t[,c(genes)])
  colnames(dat.sub) = genes
  row.names(dat.sub) = row.names(dat.t)
  dat.sub$Class = dat.t$Class
  oob = 0
  print(1)
  for (j in c(1:reps)) {
    rf = randomForest(data=dat.sub, Class~., ntree=1001, type='classification', mtry=1)
    oob = oob + as.numeric(rf$err.rate[1001,1])
  }
  oob = oob/reps
  error = data.frame(Variables=1, OOB=oob)
  for (i in c(2:nrow(imp))) {
    print(i)
    genes = row.names(subset(imp, Rank<=i))
    dat.sub = dat.t[,c(genes)]
    dat.sub$Class = dat.t$Class
    oob = 0
    for (j in c(1:reps)) {
      rf = randomForest(data=dat.sub, Class~., ntree=1001, type='classification', mtry=i)
      oob = oob + as.numeric(rf$err.rate[1001,1])
    }
    oob = oob/reps
    error = rbind(error, c(i,oob))
  }
  return(error)
}

#Construct ROC for public models
constructROC = function(pred, class) {
  pred.to.roc = pred
  pred.rocr = prediction(pred.to.roc, as.factor(class))
  perf.rocr = performance(pred.rocr, measure='auc', x.measure='cutoff')
  AUC = deparse(round(as.numeric(perf.rocr@y.values),3))
  perf.tpr.rocr = performance(pred.rocr, 'tpr', 'fpr')
  perf.precrec.rocr = performance(pred.rocr, 'prec', 'rec')
  res = list()
  res[[1]] = AUC
  res[[2]] = perf.tpr.rocr
  res[[3]] = perf.precrec.rocr
  names(res) = c('AUC','ROC','PrecRec')
  return(res)
}

### CV testing on RF models
rf_cvTest = function(df, fold, m) {
  dat = df
  nTest = round(nrow(dat)/fold,0)
  pred = c(0,0)
  class = c(0)
  
  for (k in (c(1:fold))) {
    dat.test = dat[sample(nrow(dat), nTest),]
    class = c(class, as.numeric(as.character(dat.test$Class)))
    dat.train = dat[!(row.names(dat) %in% row.names(dat.test)),]
    cv.rf = randomForest(Class~., data=dat.train, ntree=10001, type='prob')
    predi = predict(cv.rf, dat.test, type='prob')
    pred = rbind(pred, predi)
  }
  pred = pred[-1,]
  class = class[-1]
  pred.to.roc = pred[,2]
  pred.rocr = prediction(pred.to.roc, as.factor(class))
  perf.rocr = performance(pred.rocr, measure='auc', x.measure='cutoff')
  AUC = deparse(round(as.numeric(perf.rocr@y.values),3))
  perf.tpr.rocr = performance(pred.rocr,'tpr','fpr')
  perf.precrec.rocr = performance(pred.rocr,'prec','rec')
  perf.sens.rocr = performance(pred.rocr,'tpr','tnr')
  res = list()
  res[[1]] = AUC
  res[[2]] = perf.tpr.rocr
  res[[3]] = perf.precrec.rocr
  res[[4]] = perf.sens.rocr
  names(res) = c('AUC','ROC','PrecRec','SensSpe')
  return(res)
}

### Calculate area under Prec/Rec curve
calcAUPRC = function(list) {
  prerec.df = data.frame(Recall=list[[3]]@x.values[[1]], Precision=list[[3]]@y.values[[1]])
  prerec.df[is.na(prerec.df)] <- 0
  auprc = trapz(prerec.df$Recall, prerec.df$Precision)
  auprc = round(auprc,3)
  return(auprc)
}

### Test RF model on validation set
rf_valTest = function(df, model) {
  pred = c(0,0)
  for (i in c(1:nrow(df))) {
    predi = predict(model, df[i,], type='prob')
    pred = rbind(pred, predi)
  }
  pred = pred[-1,]
  pred.to.roc = pred[,2]
  pred.rocr = prediction(pred.to.roc, as.factor(df$Class))
  perf.rocr = performance(pred.rocr, measure='auc', x.measure='cutoff')
  AUC = deparse(round(as.numeric(perf.rocr@y.values),3))
  perf.tpr.rocr = performance(pred.rocr,'tpr','fpr')
  perf.precrec.rocr = performance(pred.rocr,'prec','rec')
  res = list()
  res[[1]] = AUC
  res[[2]] = perf.tpr.rocr
  res[[3]] = perf.precrec.rocr
  names(res)=c('AUC','ROC','PrecRec')
  return(res)
}

### Function to convert FPKM to binary
formatExpression = function(df, cutoff) {
  dat = df
  for (i in c(1:nrow(dat))) {
    if (dat[i,'FPKM']<cutoff) { dat[i,'FPKM']=0 }
    else { dat[i,'FPKM'] = 1 }
  }
  return(dat)
}

#------------------------------------
# Workflow
#------------------------------------

### Read in datasets
c.dat = read.delim(sep='\t', header=T, row.names=1, file='data/Chari.featureTable.txt')     # Training dataset
h.dat = read.delim(sep='\t', header=T, row.names=1, file='data/Horlbeck.featureTable.txt')  # Validation set 1

### Read in scores from other models
c.sco = read.delim(sep='\t', header=T, row.names=1, file='data/Chari.scores.txt')
h.sco = read.delim(sep='\t', header=T, row.names=1, file='data/Horlbeck.scores.txt')

### Format datasets
c.dat$Class = as.factor(c.dat$Class)
h.dat$Class = as.factor(h.dat$Class)

### Perform feature selection
rf.rf = randomForest(Class ~ .,                   data=c.dat, importance=T, ntree=10001) # Full model
ch.rf = randomForest(Class ~ . - Combined - FPKM, data=c.dat, importance=T, ntree=10001) # Chromatin only model

rf.imp = as.data.frame(varImpPlot(rf.rf)) # Collect feature importance
ch.imp = as.data.frame(varImpPlot(ch.rf))

rf.imp = rf.imp[order(-rf.imp$MeanDecreaseGini),]
ch.imp = ch.imp[order(-ch.imp$MeanDecreaseGini),]

rf.imp$Rank = c(1:nrow(rf.imp))
ch.imp$Rank = c(1:nrow(ch.imp))

rf.error = seq_rf(c.dat, rf.imp, 3) # Build sequential models, in order of feature importance
ch.error = seq_rf(c.dat, ch.imp, 3)

### Plot results of sequential error
rf.errorPlot = ggplot(rf.error, aes(x=Variables, y=OOB)) +
  geom_point() +
  labs(x='Number of features', y='Average OOB error') 
rf.errorPlot

rf.min  = as.numeric(row.names(subset(rf.error, OOB<=min(rf.error$OOB))))[1] # Find minimum number of features for best model
rf.vars = row.names(subset(rf.imp, Rank<=rf.min))

ch.min  = as.numeric(row.names(subset(ch.error, OOB<=min(ch.error$OOB))))[1]
ch.vars = row.names(subset(ch.imp, Rank<=ch.min))

### Perform 10-fold CV
rf.cvRes = rf_cvTest(subset(c.dat, select=c(rf.vars,'Class')),10,rf.min)

### Construct ROC curves for other models for comparison
ch.cvRes = rf_cvTest(subset(c.dat, select=c(ch.vars,'Class')),10,ch.min)
wu.cvRes = constructROC(c.sco$WU.CRISPR, c.sco$Class) # sgRNA only model
tr.cvRes = constructROC(c.dat$FPKM,      c.dat$Class) # Transcription only model

### Plot ROC results
plot(rf.cvRes[[2]], col='blue',   lwd=3)
plot(ch.cvRes[[2]], col='red',    lwd=3, add=TRUE)
plot(wu.cvRes[[2]], col='orange', lwd=3, add=TRUE)
plot(tr.cvRes[[2]], col='black',  lwd=3, add=TRUE)

### Construct final model
rf = randomForest(Class ~ ., data=c.dat[,c(rf.vars,'Class')], ntree=10001)

### Validate on Horlbeck dataset
rf.h.valRes = rf_valTest(h.dat, rf)
sg.h.valRes = constructROC(h.sco$sgRNAscorer.Hg, h.sco$Class)
wu.h.valRes = constructROC(h.sco$WU.CRISPR,      h.sco$Class)
az.h.valRes = constructROC(h.sco$Azimuth,        h.sco$Class)

### Plot Precision/Recall results
plot(rf.h.valRes[[3]], col='blue',      lwd=3, ylim=c(0.5,1), xlim=c(0,0.5))
plot(sg.h.valRes[[3]], col='darkgreen', lwd=3, add=TRUE)
plot(wu.h.valRes[[3]], col='orange',    lwd=3, add=TRUE)
plot(az.h.valRes[[3]], col='purple',    lwd=3, add=TRUE)
