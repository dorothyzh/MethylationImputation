# M=log2(Beta/(1-Beta))
setwd("~/Desktop/uab/R_package/simulation/element-wise_difference/")
# library(MVN)
library(MASS)
library(MetImpute)
# library(wateRmelon)
#source("https://bioconductor.org/biocLite.R")
#biocLite("impute")
library(impute)

#install.packages("mice")
library(mice)

###############################################################################################
######################################------------functions
M2Beta <-
  function (M) 
  {
    return((2^M)/(2^M + 1))
  }

Beta2M <-
  function (B) 
  {
    return(log2(B/(1 - B)))
  }

Beta_M <- function(data,...){
  data <- as.matrix(data)
  id <- which(is.na(as.vector(data)))
  Range <- range(as.vector(data)[-id])
  if(Range[1]<0 || Range[2]>1)stop("The range of values of input dataset musr be between 0 and 1 !")
  Res <- log(data/(1-data))
  # Res[Res==-Inf] = 
  return(Res)
}
####


###############################################################################################
######################################-----------------simulation
data(dat1)
data(dat2)

na_count <-sapply(as.data.frame(dat1), function(y) sum(length(which(is.na(y)))))
dat1=dat1[,-which(na_count>100)]
length(na_count)
dim(dat1)  #124 153
#dat1_1=impute.knn(dat1,colmax=1,rowmax=1,maxp=nrow(dat1)*length(dat1)) #An expression matrix with genes in the rows, samples in the columns

####imputation
dat1_1=mice(as.matrix(dat1),m=1,maxit=1,meth='sample',seed=500)
dat1=complete(dat1_1)


dat1=dat1*0.99+0.005
dat1=Beta2M(dat1)

#dat1_1=mice(as.matrix(dat1),m=1,maxit=1,meth='sample',seed=500)
#dat1=complete(dat1_1)
#covariance=cov(dat2,use="pairwise.complete.obs") ###--work
#covariance=cov(dat2,use="everything") ###work
means <- rowMeans(dat1, na.rm=T)
covariance=cov(t(dat1),use="pairwise.complete.obs")  #######work, maybe use this one 
#covariance=cov(dat2,use="na.or.complete") ###work 
dim(covariance)

covariance[is.na(covariance)] = 0
#ridge = 4
#diag(covariance) <- diag(covariance)+ridge

n = dim(dat2)[2]
n.seq = 0.1*n
m = dim(covariance)[1]
m.array = 0.1*m
Epsilon1=.1
Epsilon2=.1
s=1.05
t=0.05
#M <- mvrnorm(m, mu = rep(0, n), Sigma = diag(n)*var_a) ####can run
M <- mvrnorm(n, mu = means, Sigma = covariance) 

# M=t(M)
dim(M)  # 992 200

M <- data.frame(M)
rownames(M)=as.character(1:n)
colnames(M)=as.character(68606653:(68606652+ncol(M)))
# M <- data.frame(M)

dim(M)

seq<-M[(sample(1:n, n.seq)), ] ###sequencing
dim(seq) # 99 200

array<-M[, (sample(1:m, m.array))]  ###array    200中20 90% 
dim(array) #992  20

################################################
########----------sequencing-------------#######
################################################
noise1=rnorm(dim(seq)[1]*dim(seq)[2], mean = 0, sd = 1) ##n observations
seq1=seq+Epsilon1*noise1   #####add random noise

#par(mfrow=c(2,1))
#View(seq)
#View(seq1)

seq2=M2Beta(seq1)     #####convert M-value to beta value

##################################################
##########---------microarray-----------##########
##################################################
noise2=rnorm(dim(array)[1]*dim(array)[2], mean = 0, sd = 1) ##n observations
array1=array+Epsilon2*noise2

array2=M2Beta(array1) ####m-value to beta value

array3=s*array2+t  ###make the scale of array the same with sequencing, 0.9<s<1.1, -0.1<b<0.1

###############################################################################################
######################################-----------------see results
dat2=as.data.frame(seq2)
dim(dat2)
dat1=as.data.frame(array3)
dim(dat1)

qc_frac = 1-5/dim(dat1)[1]

seq_impute <- MetIm(sequence = dat2, microarray = dat1, lambda=0.3, cut=10, cvfold=0, use.mvalue = F, qc_frac = qc_frac)
# seq_imputet <- MetIm(sequence = t(dat1), microarray = t(dat2), lambda=0.3, cut=10, cvfold=0, use.mvalue = F, qc_frac = qc_frac)

dim(seq_impute)  ##992 200

#actual=M
#predicted=seq_impute
#se(actual, predicted)  #####compare(M,M2) ####using element-wise squared error.


# element-wise difference
seq_impute1=Beta2M(seq_impute)
diff <- seq_impute1-(as.matrix(M))
total.se <- sum(sum(diff^2,na.rm=T),na.rm=T)
total.se  ####247381.3----49721.68
total.se/(n*m)  ####0.25----0.25

###############################################################################################
######################################-----------------plot
C=as.matrix(seq2)
D=as.matrix(array3)

ROW1 = as.numeric(rownames(C))
ROW2 = as.numeric(rownames(D))
COL1 = as.numeric(colnames(C))
COL2 = as.numeric(colnames(D))

ROW = sort(union(ROW1,ROW2))
COL = sort(union(COL1,COL2))
nrow = length(ROW)
ncol = length(COL)

C_star = matrix(NA,nrow,ncol)
colnames(C_star) = COL
rownames(C_star) = ROW
D_star = C_star

id1.row = match(ROW1,ROW)
id1.col = match(COL1,COL)
id2.row = match(ROW2,ROW)
id2.col = match(COL2,COL)  

C_star[id1.row,id1.col] = as.matrix(C)
D_star[id2.row,id2.col] = as.matrix(D)

png("simulation2.png",width=10,height=10,units="in",res=600)
par(mfrow=c(2,2))
image(C); image(C_star); image(D); image(D_star)
dev.off()

png("simulation3.png",width=10,height=10,units="in",res=600)
par(mfrow=c(2,2))
image(as.matrix(t(M)),main="M"); 
image(t(C_star),main="sequence"); 
image(t(D_star),main="array"); 
image(as.matrix(t(seq_impute1)),main=paste("Reconstruction by mimpute with MEE:", sprintf("%.4f", total.se/(n*m))))
dev.off()

# vCpG
CpG_var = sapply(M, var)
vCpG <- CpG_var > 0.005;
vCpG <- CpG_var > 12;
# diff <- seq_impute1-(as.matrix(M))
total.vCpG.se <- sum(sum(diff[,vCpG]^2,na.rm=T),na.rm=T)
total.vCpG.se  ####247381.3----49721.68
total.vCpG.se/(n*sum(vCpG, na.rm=T))  ####0.25----0.25

# R2
r2 = sapply(1:length(vCpG), function(i){ 
  if (!is.na(vCpG[i]) && vCpG[i]) {
    cor(M[,i],seq_impute1[,i])^2
  }
}
)

unlist(r2)
a=M[,!is.na(vCpG) & vCpG]

plot(M[,"68606703"], seq_impute1[,"68606703"])

plot(M[,77], seq_impute1[,77])

image(as.matrix(M))
hist(as.matrix(M),breaks=100)
hist(seq_impute1,breaks=100)

