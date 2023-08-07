library(geomorph)
library(MFPCA)
library(factoextra)
library(devtools)
getwd()
setwd("C:/Users/HP/Desktop/DR ARPAH/tpsdata/dorsal")

#Open TPS File
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
plot(proc1$consensus,col="black",pch=19,cex=2)
plot(proc2$consensus,col="black",pch=19,cex=2)
plot(proc3$consensus,col="black",pch=19,cex=2)

#PCA Using GM Method###
#Dorsal 
shrew.PCA<- gm.prcomp(proc1$coords)
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
# Dorsal View

library(fda)
library(funData)
library(refund)
library(MFPCA)

shrew.FUN<- list()

for (j in 1:2) {
  objdorsal= seq(1:dim(proc1$coords[,j,])[1])
  
  objXdorsal =t(as.matrix(proc1$coords[,j,]))#each dim
  
  
  shrew.FUN[[j]] = funData(argvals =objdorsal, X=objXdorsal)
  
}
par(mfrow = c(1, 2))
plot(shrew.FUN[[1]],ylab='converted X-coordinates',main="Dimension 1")
plot(shrew.FUN[[2]],ylab='converted Y-coordinates',main="Dimension 2")
XS <- multiFunData(do.call("list",shrew.FUN)) ## multivariate smoothing

#MFPCA Using FDGM Method###
dFPCA <- MFPCA(mFData=XS, M = 9, uniExpansions = list(list(type = "uFPCA"),
                                                      list(type = "uFPCA")), fit=TRUE)
plot(dFPCA)
summary(dFPCA)
scores=dFPCA$scores
plot(scores[,1:3],col=groupnum, pch=19,xlab="MFPC 1",ylab="MFPC 2",cex=1.5)

##FUNCTIONAL
scores=data.frame(pcascores)

##Functional Linear Discriminant Analysis (Dorsal)
scores=dFPCA$scores
for (i in 1:25) {
  scores[i,]=dFPCA$scores[i,]*dFPCA[["normFactors"]]
}
fdpcascores <- scores[,1:3]
colnames(fdpcascores)=c("C1","C2","C3")
model <- lda(classlabel ~ ., data.frame(fdpcascores))
model
ldamodel <- fdpcascores %*% model$scaling
rownames(ldamodel) <- groupnum
plot(ldamodel,col=groupnum, pch=19,xlab="FLD 1",ylab="FLD 2",cex=1.5)
legend("topleft",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares 

######################
#JAW
######################

#PCA Using GM Method###
shrew.PCA_jaw<- gm.prcomp(proc2$coords)
summary(shrew.PCA_jaw)
groupnum=c(rep("red",30),rep("blue",30),rep("orange",30))
plot(shrew.PCA_jaw$x[,1:3], col=groupnum, pch=19,cex=1.5,xlab="PC 1",ylab="PC 2")
legend("bottomright", inset = c(0.04, 0.1),
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares       

#LDA using GM Method######

library(MASS)
classlabel <- groupnum

##CLASSICAL

MVscores_jaw <- shrew.PCA_jaw$x[,1:3]
colnames(MVscores_jaw)=c("C1","C2","C3")
model_jaw <- lda(classlabel~., data=data.frame(MVscores_jaw))
model_jaw
ldamodel_jaw <- MVscores_jaw  %*% model_jaw$scaling
rownames(ldamodel_jaw) <- groupnum
plot(ldamodel_jaw,col=groupnum, pch=19,xlab="LD 1",ylab="LD 2",cex=1.5)
legend("topright",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares 


#Convert Standardised Landmark Data to Multivariate Functional Data
# Jaw View

library(fda)
library(funData)
library(refund)
library(MFPCA)

shrew.FUN_jaw<- list()

for (j in 1:2) {
  obj_jaw= seq(1:dim(proc2$coords[,j,])[1])
  
  objX_jaw =t(as.matrix(proc2$coords[,j,]))#each dim
  
  
  shrew.FUN_jaw[[j]] = funData(argvals =obj_jaw, X=objX_jaw)
  
}
par(mfrow = c(1, 2))
plot(shrew.FUN_jaw[[1]],ylab='converted X-coordinates',main="Dimension 1")
plot(shrew.FUN_jaw[[2]],ylab='converted Y-coordinates',main="Dimension 2")
XS_jaw <- multiFunData(do.call("list",shrew.FUN_jaw)) ## multivariate smoothing

#MFPCA Using FDGM Method###
dFPCA_jaw <- MFPCA(mFData=XS_jaw, M = 9, uniExpansions = list(list(type = "uFPCA"),
                                                      list(type = "uFPCA")), fit=TRUE)
plot(dFPCA_jaw)
summary(dFPCA_jaw)
scores_jaw=dFPCA_jaw$scores
plot(scores_jaw[,1:3],col=groupnum, pch=19,xlab="MFPC 1",ylab="MFPC 2",cex=1.5)

##FUNCTIONAL
scores_jaw=data.frame(pcascores_jaw)

##Functional Linear Discriminant Analysis (Dorsal)
scores_jaw=dFPCA_jaw$scores
for (i in 1:50) {
  scores_jaw[i,]=dFPCA_jaw$scores[i,]*dFPCA_jaw[["normFactors"]]
}
fdpcascores_jaw <- scores_jaw[,1:3]
colnames(fdpcascores_jaw)=c("C1","C2","C3")
model_jaw <- lda(classlabel ~ ., data.frame(fdpcascores_jaw))
model_jaw
ldamodel_jaw <- fdpcascores_jaw %*% model_jaw$scaling
rownames(ldamodel_jaw) <- groupnum
plot(ldamodel_jaw,col=groupnum, pch=19,xlab="FLD 1",ylab="FLD 2",cex=1.5)
legend("topleft",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares

######################
#LATERAL
######################

#PCA Using GM Method###
shrew.PCA_lat<- gm.prcomp(proc3$coords)
summary(shrew.PCA_lat)
groupnum=c(rep("red",30),rep("blue",30),rep("orange",30))
plot(shrew.PCA_lat$x[,1:3], col=groupnum, pch=19,cex=1.5,xlab="PC 1",ylab="PC 2")
legend("bottomright", inset = c(0.04, 0.1),
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares       

#LDA using GM Method######

library(MASS)
classlabel <- groupnum

##CLASSICAL

MVscores_lat <- shrew.PCA_lat$x[,1:3]
colnames(MVscores_jaw)=c("C1","C2","C3")
model_lat <- lda(classlabel~., data=data.frame(MVscores_lat))
model_lat
ldamodel_lat <- MVscores_lat  %*% model_lat$scaling
rownames(ldamodel_lat) <- groupnum
plot(ldamodel_lat,col=groupnum, pch=19,xlab="LD 1",ylab="LD 2",cex=1.5)
legend("topright",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares 


#Convert Standardised Landmark Data to Multivariate Functional Data
# Lateral View

library(fda)
library(funData)
library(refund)
library(MFPCA)

shrew.FUN_lat<- list()

for (j in 1:2) {
  obj_lat= seq(1:dim(proc3$coords[,j,])[1])
  
  objX_lat =t(as.matrix(proc3$coords[,j,]))#each dim
  
  
  shrew.FUN_lat[[j]] = funData(argvals =obj_lat, X=objX_lat)
  
}
par(mfrow = c(1, 2))
plot(shrew.FUN_lat[[1]],ylab='converted X-coordinates',main="Dimension 1")
plot(shrew.FUN_lat[[2]],ylab='converted Y-coordinates',main="Dimension 2")
XS_lat <- multiFunData(do.call("list",shrew.FUN_lat)) ## multivariate smoothing

#MFPCA Using FDGM Method###
dFPCA_lat <- MFPCA(mFData=XS_lat, M = 9, uniExpansions = list(list(type = "uFPCA"),
                                                              list(type = "uFPCA")), fit=TRUE)
plot(dFPCA_lat)
summary(dFPCA_lat)
scores_lat=dFPCA_lat$scores
plot(scores_lat[,1:3],col=groupnum, pch=19,xlab="MFPC 1",ylab="MFPC 2",cex=1.5)

##FUNCTIONAL
scores_lat=data.frame(pcascores_lat)

##Functional Linear Discriminant Analysis (Dorsal)
scores_lat=dFPCA_lat$scores
for (i in 1:47) {
  scores_lat[i,]=dFPCA_lat$scores[i,]*dFPCA_lat[["normFactors"]]
}
fdpcascores_lat <- scores_lat[,1:3]
colnames(fdpcascores_lat)=c("C1","C2","C3")
model_lat <- lda(classlabel ~ ., data.frame(fdpcascores_lat))
model_lat
ldamodel_lat <- fdpcascores_lat %*% model_lat$scaling
rownames(ldamodel_lat) <- groupnum
plot(ldamodel_lat,col=groupnum, pch=19,xlab="FLD 1",ylab="FLD 2",cex=1.5)
legend("topleft",
       legend = c("Suncus murinus", "Crocidura malayana","Crocidura monticola"),
       fill=c("red","blue","orange"),       # Color of the squares
       border = "black") # Color of the border of the squares

#####################################################################
####MACHINE LEARNING ALGORITHMS######################################
#Repeat steps using PC components obtained from jaw and lateral views
#####################################################################
library(e1071)
library(caret)

#DORSAL
p=3
predPCA <- predict(dFPCAdorsal)
pcascores=cbind(predPCA[[1]]@X[,1:p],predPCA[[2]]@X[,1:p])
pcaFDGM<-data.frame(pcascores, y = as.factor(groupnum))
pcaGM<-data.frame(shrew.PCA$x[,1:p], y = as.factor(groupnum))

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

#####################################################################
#Repeat steps using PC components obtained from jaw and lateral views
#####################################################################