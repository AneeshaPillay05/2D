library(geomorph)
library(MFPCA)
library(factoextra)
library(devtools)
getwd()
setwd("C:/Users/HP/Desktop/DR ARPAH/tpsdata/dorsal")

filelist<- list.files(pattern = ".TPS")
dorsal<-readmulti.tps(filelist)

setwd("C:/Users/HP/Desktop/DR ARPAH/tpsdata/jaws")
filelist1<- list.files(pattern = ".TPS")
jaws<-readmulti.tps(filelist1)

setwd("C:/Users/HP/Desktop/DR ARPAH/tpsdata/lateral")
filelist2<- list.files(pattern = ".TPS")
lateral<-readmulti.tps(filelist2)


#Generalised Procrustes Analysis for each view
proc1<-gpagen(dorsal)
proc2<-gpagen(jaws)
proc3<-gpagen(lateral)

#Combine all three angles
comb.lm <- combine.subsets(dorsal = proc1, jaws=proc2, lateral=proc3, gpa = TRUE)

#PCA Using GM Method###
shrew.PCA<- gm.prcomp(comb.lm$coords)
summary(shrew.PCA)
groupnum=c(rep("red",30),rep("blue",30),rep("orange",30))
plot(shrew.PCA$x[,1:3], col=groupnum, pch=19,cex=1.5,xlab="PC 1",ylab="PC 2")
legend("bottomright", inset = c(0.04, 0.1),
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares       

#LDA using GM Method######

library(MASS)
classlabel <- groupnum

##CLASSICAL

MVscore.train <- shrew.PCA$x[,1:3]
colnames(MVscores )=c("C1","C2","C3")
model <- lda(classlabel~., data=data.frame(MVscores))
model
ldamodel <- MVscores %*% model$scaling
rownames(ldamodel) <- groupnum
par(mfrow = c(1,1))
plot(ldamodel,col=groupnum, pch=19,xlab="LD 1",ylab="LD 2",cex=1.5)
legend("topright",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares 

#Convert Standardised Landmark Data to Multivariate Functional Data
library(fda)
library(funData)
library(refund)
library(MFPCA)

#all coordinates (all three angles) are combined based on separate dimensions

shrew.FUN<- list()

for (j in 1:2) {
  objdorsal= seq(1:dim(proc1$coords[,j,])[1])
  
  objXdorsal =t(as.matrix(proc1$coords[,j,]))#each dim
  
  
  shrew.FUN[[j]] = funData(argvals =objdorsal, X=objXdorsal)
  
}

for (j in 3:4) {
  objjaws= seq(1:dim(proc2$coords[,(j-2),])[1])
  
  
  objXjaws =t(as.matrix(proc2$coords[,(j-2),]))#each dim
  
  shrew.FUN[[j]] = funData(argvals =objjaws, X=objXjaws)
  
}

for (j in 5:6) {
  objlateral= seq(1:dim(proc3$coords[,(j-4),])[1])
  
  
  objXlateral =t(as.matrix(proc3$coords[,(j-4),]))#each dim
  
  shrew.FUN[[j]] = funData(argvals =objlateral, X=objXlateral)
  
}

XS=multiFunData(shrew.FUN)

#MFPCA
dFPCA <- MFPCA(mFData=XS, M = 30, uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"),list(type = "uFPCA"),list(type = "uFPCA"),list(type = "uFPCA"),list(type = "uFPCA")), fit=TRUE)

cumsum(dFPCA$values)/sum(dFPCA$values)
plot(dFPCA$scores[,1],dFPCA$scores[,2], col=groupnum,pch=19, cex=2,xlab="MFPC 1", ylab="MFPC 2" )
legend("topleft",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares       


#######LDA###############
library(MASS)
classlabel <- groupnum

##CLASSICAL

MVscores_comb <- shrew.PCA$x[,1:3]
colnames(MVscores_comb)=c("C1","C2","C3")
model_comb <- lda(classlabel~., data=data.frame(MVscores_comb))
model_comb
ldamodel_comb <- MVscores.train_comb  %*% model_comb$scaling
rownames(ldamodel_comb) <- groupnum

#LD1-LD2
par(mfrow = c(1,1))
plot(ldamodel_comb,col=groupnum, pch=19, xlab="LD 1",ylab="LD 2",cex=1.5)
legend("topright",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares 

##FUNCTIONAL (Using MFPC Scores)
p=3
scores=data.frame(dFPCA$scores[,1:p])
scores=dFPCA$scores
for (i in 1:122) {
  scores[i,]=dFPCA$scores[i,]*dFPCA[["normFactors"]]
}
fdpcascores <- scores[,1:3]
colnames(fdpcascores)=c("C1","C2","C3")
model <- lda(classlabel ~ ., data.frame(fdpcascores))
model

ldamodel <- fdpcascores %*% model$scaling
rownames(ldamodel) <- groupnum
plot(ldamodel,col=groupnum, pch=19,xlab="FLD 1",ylab="FLD 2",cex=1.5)
legend("bottomleft",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares 

##FUNCTIONAL (Using Predicted Data)

fdpcascores <- dFPCA$scores[,1:3]
colnames(fdpcascores)=c("C1","C2","C3")
model <- lda(classlabel ~ ., data.frame(fdpcascores))
model
#type model to check lda output

ldamodel <- fdpcascores %*% model$scaling
rownames(ldamodel) <- groupnum
plot(ldamodel,col=groupnum, pch=19,xlab="FLD 1",ylab="FLD 2",cex=1.5)
legend("bottomleft",
        legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
        fill=c("red","blue","orange"),       # Color of the squares
        border = "black") # Color of the border of the squares 




#Machine Learning using predicted data
p=3
predPCA<-predict(dFPCA)
pcascores=cbind(predPCA[[1]]@X[,1:p],predPCA[[2]]@X[,1:p])
pcaFDGM<-data.frame(pcascores, y = as.factor(groupnum))
pcaGM<-data.frame(shrew.PCA$x[,1:p], y = as.factor(groupnum))

#Naive Bayes-The categorical Naive Bayes classifier is 
#suitable for classification with discrete features that are categorically distributed. 
library(e1071)
library(caret)

vector1 <- c();
vector2 <- c();
M=20
for (i in 1:M) {
  ## Generate training and test data
  trainID<-sample(1:nrow(pcaFDGM),63)
  testID<-setdiff(1:nrow(pcaFDGM),trainID)
  sample_trainFDGM<-pcaFDGM[trainID,]
  sample_testFDGM<-pcaFDGM[testID,]
  sample_trainGM<-pcaGM[trainID,]
  sample_testGM<-pcaGM[testID,]
  
  ## Fit a model on training data
  #Classification Techniques
  #Naive Bayes (FDGM Vs GM)
  tcontrol <- trainControl(method = "cv", number = 10)
  modelnbFDGM <- train(y~., data = sample_trainFDGM, method = "nb", trControl = tcontrol)
  modelnbGM <- train(y~., data = sample_trainGM, method = "nb", trControl = tcontrol)
  
  ## Predict on test  data
  #NB
  spNBFDGM <- predict(modelnbFDGM, sample_testFDGM)
  spNBGM <- predict(modelnbGM, sample_testGM)
  test_tableFDGM<- table(spNBFDGM,sample_testFDGM$y)
  correct_classFDGM<-sum(diag(test_tableFDGM))/27
  test_tableGM<- table(spNBGM,sample_testGM$y)#81.5% of correct classification
  correct_classGM<-sum(diag(test_tableGM))/27
  
  ## Save the result in a vector
  vector1 <- c(vector1,correct_classFDGM)
  vector2<- c(vector2,correct_classGM)
  
}
mean(vector1)
sd(vector1)
mean(vector2)
sd(vector2)

#SVM
vector3 <- c();
vector4 <- c();
#M=2
for (i in 1:M) {
  ## Generate training and test data
  trainID<-sample(1:nrow(pcaFDGM),63)
  testID<-setdiff(1:nrow(pcaFDGM),trainID)
  sample_trainFDGM<-pcaFDA[trainID,]
  sample_testFDGM<-pcaFDGM[testID,]
  sample_trainGM<-pcaGM[trainID,]
  sample_testGM<-pcaGM[testID,]
  
  ## Fit a model on training data
  #Classification Techniques
  tcontrol <- trainControl(method = "cv", number = 10)
  modelsvmFDGM <- train(y~., data = sample_trainFDGM, method = "svmRadial", trControl =  tcontrol,  preProcess = c("center","scale"))
  modelsvmGM <- train(y~., data = sample_trainGM, method = "svmRadial", trControl =  tcontrol,  preProcess = c("center","scale"))
  
  ## Predict on test  data
  spsvmFDGM <- predict(modelsvmFDGM, sample_testFDGM)
  spsvmGM <- predict(modelsvmGM, sample_testGM)
  test_tableFDGM<- table(spsvmFDGM,sample_testFDGM$y)#81.5% of correct classification
  correct_classFDGM<-sum(diag(test_tableFDGM))/27
  test_tableGM<- table(spsvmFDGM,sample_testGM$y)#81.5% of correct classification
  correct_classGM<-sum(diag(test_tableGM))/27
  
  ## Save the result in a vector
  vector3 <- c(vector3,correct_classFDGM)
  vector4<- c(vector4,correct_classGM)
  
}
mean(vector3)
sd(vector3)
mean(vector4)
sd(vector4)

#RF

vector5 <- c();
vector6 <- c();
#M=2
for (i in 1:M) {
  ## Generate training and test data
  trainID<-sample(1:nrow(pcaFDGM),63)
  testID<-setdiff(1:nrow(pcaFDGM),trainID)
  sample_trainFDGM<-pcaFDGM[trainID,]
  sample_testFDGM<-pcaFDGM[testID,]
  sample_trainGM<-pcaGM[trainID,]
  sample_testGM<-pcaGM[testID,]
  
  ## Fit a model on training data
  #Classification Techniques
  tcontrol <- trainControl(method = "cv", number = 10)
  modelRFFDGM <- train(y~., data = sample_trainFDGM, method = "rf", ntree = 100, importance = T, trControl = tcontrol)
  modelRFGM <- train(y~., data = sample_trainGM, method = "rf", ntree = 100, importance = T, trControl = tcontrol)
  
  ## Predict on test  data
  
  spRFFDGM <- predict(modelRFFDGM, sample_testFDGM)
  spRFMV <- predict(modelRFGM, sample_testGM)
  test_tableFDGM<- table(spRFFDGM,sample_testFDGM$y)#81.5% of correct classification
  correct_classFDGM<-sum(diag(test_tableFDGM))/27
  test_tableGM<- table(spRFMV,sample_testGM$y)#81.5% of correct classification
  correct_classGM<-sum(diag(test_tableGM))/27
  
  ## Save the result in a vector
  vector5 <- c(vector5,correct_classFDGM)
  vector6<- c(vector6,correct_classGM)
  
}


#GLM


vector7 <- c();
vector8 <- c();

for (i in 1:M) {
  ## Generate training and test data
  trainID<-sample(1:nrow(pcaFDGM),63)
  testID<-setdiff(1:nrow(pcaFDGM),trainID)
  sample_trainFDGM<-pcaFDA[trainID,]
  sample_testFDGM<-pcaFDGM[testID,]
  sample_trainGM<-pcaGM[trainID,]
  sample_testMV<-pcaMV[testID,]
  
  ## Fit a model on training data
  #Classification Techniques
  tcontrol <- trainControl(method = "cv", number = 10)
  modelglmFDGM <- train(y~., data = sample_trainFDGM, method = "glmnet",  trControl = tcontrol)
  modelglmGM <- train(y~., data = sample_trainFDGM, method = "glmnet",  trControl = tcontrol)
  
  ## Predict on test  data
  
  spGLMFDGM <- predict(modelglmFDA, sample_testFDGM)
  spGLMGM <- predict(modelglmGM, sample_testGM)
  test_tableFDGM<- table(spGLMFDGM,sample_testFDGM$y)#81.5% of correct classification
  correct_classFDGM<-sum(diag(test_tableFDGM))/27
  test_tableGM<- table(spGLMGM,sample_testMV$y)#81.5% of correct classification
  correct_classMV<-sum(diag(test_tableMV))/27
  
  ## Save the result in a vector
  vector7 <- c(vector1,correct_classFDGM)
  vector8<- c(vector2,correct_classGM)
}
mean(vector7)
sd(vector7)
mean(vector8)
sd(vector8)