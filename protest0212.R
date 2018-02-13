train=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_pair_prediction_u_18f_train",header = TRUE)
test=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_pair_prediction_u_18f_test",header = TRUE)
val=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_pair_prediction_u_18f_val",header = TRUE)

library(glinternet)
library(e1071)
library(randomForest)
library(glmnet)

ligand=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\Ligm1_data.txt",header = TRUE)
receptor=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\1f3cAB_data.txt",header = TRUE)
ligand
receptor
train[1:100,]
sum(receptor$resabsEA==0)

unique(train[,1])
unique(test[,1])
unique(val[,1])

colnames(ligand)=c("ProtName_L","resSeq_L","resNo_L","resabsEA_L","resrelEA_L","resEC_L","resIC_L","resEV_L","resH1_L","resH2_L","respKa1_L","respKa2_L")
colnames(receptor)=c("ProtName_R","resSeq_R","resNo_R","resabsEA_R","resrelEA_R","resEC_R","resIC_R","resEV_R","resH1_R","resH2_R","respKa1_R","respKa2_R")
colnames(train)
protest=merge(ligand,receptor)
protest1=merge(ligand,receptor)

#########################
trainf=which(train$Flag==1)
trainn=which(train$Flag==0)
sampletimes=1
set.seed(212)
sampletrain=sample(trainn,replace=F,size=length(trainf)*sampletimes)
max(sampletrain)
trainset=train[c(trainf,sampletrain),]

testf=which(test$Flag==1)
testn=which(test$Flag==0)
sampletimes=1
set.seed(212)
sampletest=sample(testn,replace=F,size=length(testf)*sampletimes)
max(sampletest)
testset=test[c(testf,sampletest),]

valf=which(val$Flag==1)
valn=which(val$Flag==0)
sampletimes=1
set.seed(212)
sampleval=sample(valn,replace=F,size=length(valf)*sampletimes)
max(sampleval)
valset=val[c(valf,sampleval),]

rm(train,test,val)
rm(sampletest,sampletrain,sampleval,testf,testn,trainf,trainn,valf,valn)
#############################

train=rbind(trainset,valset,testset)
rm(trainset,testset,valset)

###################################

protest[1,]
train[1,]
protest=protest[,c(1,14,2,16:24,4:12)]
colnames(protest)[1:3]=c("ProtName","resName_R","resName_L")
#############################
levels(protest$resName_L)=levels(protest$resName_R)
levels(train$resName_L)
levels(train$resName_R)

train[,2]=as.numeric(train[,2])-1
train[,3]=as.numeric(train[,3])-1

protest[,2]=as.numeric(protest[,2])-1
protest[,3]=as.numeric(protest[,3])-1

################################glinternet

numLevels = c(rep(1, 18))
fit1 = glinternet(train[,c(2:3,7:9,11:18,20:24)],train[,25],lambda=1,numLevels, family="binomial")
fit1 = glinternet(train[,c(2:3,7:9,11:18,20:24)],train[,25],lambda=fit1$lambda[2]/20,numLevels,screenLimit = 20, family="binomial")

prevalgli=predict(fit1,as.matrix(protest[,c(2:7,9:16,18:21)]),type="response")[,2]
FP=10
glinum=which(prevalgli>sort(prevalgli,decreasing = T)[FP])
protest1[glinum,]

#################################gl
cvfit1=cv.glmnet(as.matrix(train[,c(2:3,7:9,11:18,20:24)]),train[,25], family="binomial")
fit1 = glmnet(as.matrix(train[,c(2:3,7:9,11:18,20:24)]),train[,25], family="binomial",lambda=cvfit1$lambda.min)

prevalgl=predict(fit1,as.matrix(protest[,c(2:7,9:16,18:21)]),type="response")
FP=10
glnum=which(prevalgl>sort(prevalgl,decreasing = T)[FP])
protest1[glnum,]

##################################rf
rf1=randomForest(train[,c(2:3,7:9,11:18,20:24)], train[,25])

prerf=predict(rf1,as.matrix(protest[,c(2:7,9:16,18:21)]),type="response")
FP=10
rfnum=which(prerf>sort(prerf,decreasing = T)[FP])
protest1[rfnum,]
#####################################svm
sv2=svm(train[,c(2:3,7:9,11:18,20:24)],train[,25])

prevalsvm1=predict(sv2,as.matrix(protest[,c(2:7,9:16,18:21)]),type="response")
FP=10
svmnum=which(prevalsvm1>sort(prevalsvm1,decreasing = T)[FP])
protest1[svmnum,]

###################################尝试分类计算
rm(list=ls())
train=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_pair_prediction_u_18f_train",header = TRUE)
test=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_pair_prediction_u_18f_test",header = TRUE)
val=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_pair_prediction_u_18f_val",header = TRUE)
library(glinternet)
library(e1071)
library(randomForest)
library(glmnet)
ligand=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\Ligm1_data.txt",header = TRUE)
receptor=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\1f3cAB_data.txt",header = TRUE)
colnames(ligand)=c("ProtName_L","resSeq_L","resNo_L","resabsEA_L","resrelEA_L","resEC_L","resIC_L","resEV_L","resH1_L","resH2_L","respKa1_L","respKa2_L")
colnames(receptor)=c("ProtName_R","resSeq_R","resNo_R","resabsEA_R","resrelEA_R","resEC_R","resIC_R","resEV_R","resH1_R","resH2_R","respKa1_R","respKa2_R")
protest=merge(ligand,receptor)
protest1=merge(ligand,receptor)
protest=protest[,c(1,14,2,16:24,4:12)]
colnames(protest)[1:3]=c("ProtName","resName_R","resName_L")

###########################找出氨基酸对数少于20000的
a=c()
for(i in 1:length(unique(train[,1])))
a[i]=sum(train[,1]==unique(train[,1])[i])
aa=which(a<=20000)
aaa=which(train[,1]%in%unique(train[,1])[aa])
train=train[aaa,]
b=c()
for(i in 1:length(unique(test[,1])))
  b[i]=sum(test[,1]==unique(test[,1])[i])
bb=which(b<=20000)
bbb=which(test[,1]%in%unique(test[,1])[bb])
test=test[bbb,]
d=c()
for(i in 1:length(unique(val[,1])))
  d[i]=sum(val[,1]==unique(val[,1])[i])
dd=which(d<=20000)
ddd=which(val[,1]%in%unique(val[,1])[dd])
val=val[ddd,]
rm(a,aa,aaa,b,bb,bbb,d,dd,ddd)
################
#########################
trainf=which(train$Flag==1)
trainn=which(train$Flag==0)
sampletimes=1
set.seed(212)
sampletrain=sample(trainn,replace=F,size=length(trainf)*sampletimes)
max(sampletrain)
trainset=train[c(trainf,sampletrain),]

testf=which(test$Flag==1)
testn=which(test$Flag==0)
set.seed(212)
sampletest=sample(testn,replace=F,size=length(testf)*sampletimes)
max(sampletest)
testset=test[c(testf,sampletest),]

valf=which(val$Flag==1)
valn=which(val$Flag==0)
set.seed(212)
sampleval=sample(valn,replace=F,size=length(valf)*sampletimes)
max(sampleval)
valset=val[c(valf,sampleval),]

rm(train,test,val)
rm(sampletest,sampletrain,sampleval,testf,testn,trainf,trainn,valf,valn)
train=rbind(trainset,valset,testset)
rm(trainset,testset,valset)
train[,2]=as.numeric(train[,2])-1
train[,3]=as.numeric(train[,3])-1
protest[,2]=as.numeric(protest[,2])-1
protest[,3]=as.numeric(protest[,3])-1
####################预测过程同第一步，不重复。

####################预测界面点
train=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_prediction_u_9f_train",header = TRUE)
test=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_prediction_u_9f_test",header = TRUE)
val=read.table("D:/旧电脑资料/E盘/尹建鑫/蛋白质数据/20161126数据/Interface_residue_prediction_u_9f_val",header = TRUE)

library(glinternet)
library(e1071)
library(randomForest)
library(glmnet)

ligand=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\Ligm1_data.txt",header = TRUE)
receptor=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\1f3cAB_data.txt",header = TRUE)
colnames(ligand)=c("ProtName","resName","resSeq","absEA","relEA","EC","IC","EV","H1","H2","pKa1","pKa2")
colnames(receptor)=c("ProtName","resName","resSeq","absEA","relEA","EC","IC","EV","H1","H2","pKa1","pKa2")

train=rbind(train,test,val)
rm(test,val)
train[,2]=as.numeric(train[,2])-1
ligand[,2]=as.numeric(ligand[,2])-1
receptor[,2]=as.numeric(receptor[,2])-1
ligand1=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\Ligm1_data.txt",header = TRUE)
receptor1=read.table("C:\\Users\\10630\\Desktop\\蛋白质竞赛\\1f3cAB_data.txt",header = TRUE)

train[1,]
ligand[1,]
receptor[1,]

######################
numLevels = c(rep(1, 9))
fit1 = glinternet(train[,c(2,6:8,10:14)],train[,15],lambda=1,numLevels, family="binomial")
fit1 = glinternet(train[,c(2,6:8,10:14)],train[,15],lambda=fit1$lambda[2]/20,numLevels,screenLimit = 20, family="binomial")

preglil=predict(fit1,as.matrix(ligand[,c(2,4:7,9:12)]),type="response")[,2]
preglir=predict(fit1,as.matrix(receptor[,c(2,4:7,9:12)]),type="response")[,2]
FP=5
glilnum=which(preglil>=sort(preglil,decreasing = T)[FP])
glirnum=which(preglir>=sort(preglir,decreasing = T)[FP])
ligand1[glilnum,]
receptor1[glirnum,]
#################################gl
cvfit1=cv.glmnet(as.matrix(train[,c(2,6:8,10:14)]),train[,15], family="binomial")
fit1 = glmnet(as.matrix(train[,c(2,6:8,10:14)]),train[,15], family="binomial",lambda=cvfit1$lambda.min)

pregll=predict(fit1,as.matrix(ligand[,c(2,4:7,9:12)]),type="response")
preglr=predict(fit1,as.matrix(receptor[,c(2,4:7,9:12)]),type="response")
FP=5
gllnum=which(pregll>=sort(pregll,decreasing = T)[FP])
glrnum=which(preglr>=sort(preglr,decreasing = T)[FP])
ligand1[gllnum,]
receptor1[glrnum,]

##################################rf
rf1=randomForest(train[,c(2,6:8,10:14)], train[,15])

prerfl=predict(rf1,as.matrix(ligand[,c(2,4:7,9:12)]),type="response")
prerfr=predict(rf1,as.matrix(receptor[,c(2,4:7,9:12)]),type="response")
FP=5
rflnum=which(prerfl>=sort(prerfl,decreasing = T)[FP])
rfrnum=which(prerfr>=sort(prerfr,decreasing = T)[FP])
ligand1[rflnum,]
receptor1[rfrnum,]

#####################################svm
sv2=svm(train[,c(2,6:8,10:14)],train[,15])

presvml=predict(sv2,as.matrix(ligand[,c(2,4:7,9:12)]),type="response")
presvmr=predict(sv2,as.matrix(receptor[,c(2,4:7,9:12)]),type="response")
FP=5
svmlnum=which(presvml>=sort(presvml,decreasing = T)[FP])
svmrnum=which(presvmr>=sort(presvmr,decreasing = T)[FP])
ligand1[svmlnum,]
receptor1[svmrnum,]

