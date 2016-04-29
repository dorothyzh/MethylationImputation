
library(devtools)
install_github("dorothyzh/MethyImpute2",subdir="MethyImpute")

library(plotrix)
library(coxme)
library(kinship2)
library(nlme)
# library(MethyImpute2)
library(MetImpute)
data(dat1)
data(dat2)
dim(dat1)
dim(dat2)
qc_frac = 1-5/dim(dat2)[2]
Seq <- MetIm(sequence = t(dat1), microarray = t(dat2), lambda=0.3, cut=10, cvfold=0, use.mvalue = T, qc_frac = qc_frac)

par(mfrow=c(2,1))
image(dat1)
image(t(Seq))

# SIMULATION
n = 1000
n.seq = 100
m = 200
m.array = 20 

# generate m-by-n "true" data by MVN
M <- MVN(mu, Sigma) # raw M-value
B <- inv-logit(C) # raw beta-value

seq<-M[order(sample(1:n, n.seq)), ]
array<-M[, order(sample(1:m, m.array))] 200中20 90%
qc_frac = 1-5/dim(array)[2]
Seq <- MetIm(sequence = t(seq), microarray = t(array), lambda=0.3, cut=10, cvfold=0, use.mvalue = T, qc_frac = qc_frac)

Accuracy = compare(Seq, B)

# plot results
setwd("~/Desktop/uab/R_package/simulation/")
seq_o=read.csv("../Fwd__Methylation_Imputation_Project/Build R Package MetImpute/Sequencing.csv")
array_o=read.csv("../Fwd__Methylation_Imputation_Project/Build R Package MetImpute/MicroArray.csv")
cov2cor(cov(seq))## cov2cor() scales a covariance matrix by its diagonal
##           to become the correlation matrix
#install.packages("MVN")
#install.packages("arm")
#install.packages("compare")
library(arm)
library(MVN)
library(MASS)
library(compare)
Sigma <- matrix(c(10,4,4,2),2,2)
Sigma
mvrnorm(n=2,mu=c(0,0),Sigma)

n <- 100
x1 <- rnorm (n)
x2 <- rbinom (n, 1, .5)
b0 <- 1
b1 <- 1.5
b2 <- 2
Inv.logit <- invlogit(b0+b1*x1+b2*x2)
plot(b0+b1*x1+b2*x2, Inv.logit)

Sigma2 <- matrix(c(10,3,3,2),2,2)
Sigma2
mu=rep(0,2)
mvrnorm(n=1000, mu=mu, Sigma)

install.packages("Metrics")
library(Metrics)
install.packages("DAAG")
library(DAAG)
install.packages("matrixcalc")
library(matrixcalc)
#install.packages("genPositiveDefMat")
#library(genPositiveDefMat)
install.packages("clusterGeneration")
library(clusterGeneration)
source("https://bioconductor.org/biocLite.R")
biocLite("wateRmelon") 
library(wateRmelon)
library(MetImpute)
#genPositiveDefMat
#beta2m
#m2beta

se(actual, predicted)
se(Sigma,Sigma2)
##########################################################################################################################
############################------------------0427----------------------##################################################
##########################################################################################################################
M=log2(Beta/(1-Beta))

nmz = 1000 #number of CpGs, here e.g. 1000
var_a = 0.5 #variance 0.5.
M <- mvrnorm(200, mu = rep(0, nmz), Sigma = diag(nmz)*var_a)  ####200samples


n = 1000  ###1000samples
n.seq = 100
m = 200
m.array = 20 

seq<-M[order(sample(1:n, n.seq)), ] ###sequencing
dim(seq) # 100 1000

array<-M[, order(sample(1:m, m.array))]  ###array    200中20 90% 
dim(array) #200  20


####sequencing
noise1=rnorm(dim(seq)[1]*dim(seq)[2], mean = 0, sd = 1) ##n observations
Epsilon1=1
seq1=seq+Epsilon1*noise1   #####add random noise

seq2=m2beta(seq1)     #####convert M-value to beta value

####microarray
noise2=rnorm(dim(array)[1]*dim(array)[2], mean = 0, sd = 1) ##n observations
Epsilon2=1
array1=array+Epsilon2*noise2

array2=m2beta(array1) ####m-value to beta value

s=1.05
t=0.05
array3=s*array2+t  ###make the scale of array the same with sequencing, 0.9<s<1.1, -0.1<b<0.1

#####see results
dat1=as.data.frame(seq2)
dim(dat1)
dat2=as.data.frame(array3)
dim(dat2)
qc_frac = 1-5/dim(dat2)[2]
seq_impute <- MetIm(sequence = t(dat1), microarray = t(dat2), lambda=0.3, cut=10, cvfold=0, use.mvalue = F, qc_frac = qc_frac)

actual=seq2
predicted=seq_impute
se(actual, predicted)  #####compare(M,M2) ####using element-wise squared error.




################draft
Sigma=matrix(c(1:4),nrow=2,ncol=2)  ###row---samples, col---CpGs


sigma=genPositiveDefMat(dim=5)
,covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
alphad=1, eta=1, rangeVar=c(1,10), lambdaLow=1, ratioLambda=10)

s<-matrix(1:25,5)
s <- matrix(c(1:16),nrow=4, ncol = 4)
s[lower.tri(s)] = t(s)[lower.tri(s)]
sigma=s
is.positive.definite(sigma)
mu=c(10, 5, 3,1)
M=mvrnorm(n=1000, mu=mu, sigma)  ####no missing data, raw data, M-value~MVN distribution

#I want to simualte data from -N- (e.g. = 3) multivariate distributions with -fam- (e.g. 10) persons with sigma tau.a.
nmz = 1000 #number of families, here e.g. 10
var_a = 0.5 #tau.g in the script. a variance of tau.a (e.g. 0.5).
a2_mz <- mvrnorm(200, mu = rep(0, nmz), Sigma = diag(nmz)*var_a)
#a1_mz<-array(dim=c(dim(a2_mz),3))
#for(i in 1:3)  a1_mz[,,i]<-mvrnorm(3,t(a2_mz)[,i],diag(nmz)*var_a)



######################################## examination ---35 line begin in MetIm.R
#################old data
data(dat1)
data(dat2)
dim(dat1)
dim(dat2)
qc_frac = 1-5/dim(dat2)[2]

sequence = t(dat1)
microarray = t(dat2)
lambda=0.3
cut=10
cvfold=0
use.mvalue = T
qc_frac = qc_frac


#################new data
dat1=as.data.frame(seq2)
dim(dat1)
dat2=as.data.frame(array3)
dim(dat2)
qc_frac = 1-5/dim(dat2)[2]

sequence = t(dat1)
microarray = t(dat2)
lambda=0.3
cut=10
cvfold=0
use.mvalue = F
qc_frac = qc_frac


