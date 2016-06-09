
###########################################################################################
#######function
simul=function(ep1,ep2,outputName){

library(MASS)
source("MethylationImputation/MethyImpute/R/MetIm.R")

######################################
####-------sub functions
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


###########################################
##########-----------------simulation
load("MethylationImputation/MethyImpute/data/dat2.rda")
load("MethylationImputation/MethyImpute/data/dat1.rda")
means <- rowMeans(dat1, na.rm=T)
covariance=cov(t(dat1),use="pairwise.complete.obs")  #######work, maybe use this one 
dim(covariance)  #124 124

covariance[is.na(covariance)] = 0
ridge = 1e-3
diag(covariance) <- diag(covariance)+ridge

n = dim(dat2)[2]
n.seq = 0.1*n
m = dim(covariance)[1]
m.array = 0.1*m
Epsilon1=ep1
Epsilon2=ep2
s=1.05
t=0.05
M <- mvrnorm(n, mu = means, Sigma = covariance) 

dim(M)  # 992 124

M <- data.frame(M)
rownames(M)=as.character(1:n)
colnames(M)=as.character(68606653:(68606652+ncol(M)))

seq<-M[(sample(1:n, n.seq)), ] ###sequencing
dim(seq) # 99 124

array<-M[, (sample(1:m, m.array))]  ###array    200中20 90% 
dim(array) #992  12

################################################
########----------sequencing-------------#######
################################################
noise1=rnorm(dim(seq)[1]*dim(seq)[2], mean = 0, sd = 1) ##n observations
seq1=seq+Epsilon1*noise1   #####add random noise

seq2=M2Beta(seq1)     #####convert M-value to beta value

##################################################
##########---------microarray-----------##########
##################################################
noise2=rnorm(dim(array)[1]*dim(array)[2], mean = 0, sd = 1) ##n observations
array1=array+Epsilon2*noise2

array2=M2Beta(array1) ####m-value to beta value

array3=s*array2+t  ###make the scale of array the same with sequencing, 0.9<s<1.1, -0.1<b<0.1

###############################################
############-----------------see results
dat1=as.data.frame(seq2)  #### 99 124
dim(dat1)
dat2=as.data.frame(array3)  #####992 12
dim(dat2)

qc_frac = 1-5/dim(dat2)[1]



ALLlambda=seq(from=0.01,to=0.5,by=0.01)
ALLlambda

results=data.frame(lambda=ALLlambda,totalSePerElement=rep(0,length(ALLlambda)))
for (i in 1:length(ALLlambda)){
  lambda=ALLlambda[i]
  seq_impute <- MetIm(sequence = dat1, microarray = dat2, lambda=lambda, cut=10, cvfold=0, use.mvalue = F, qc_frac = qc_frac)
  dim(seq_impute)  ##992 124
  
  # element-wise difference
  seq_impute1=Beta2M(seq_impute)
  diff <- seq_impute1-(as.matrix(M))
  total.se <- sum(sum(diff^2,na.rm=T),na.rm=T)
  #total.se  ####379.908
  total.se/(n*m)  ####0.003088482
  
  results[i,2]=total.se/(n*m) 
  
}

write.csv(results,paste("0519results",outputName,".csv",sep=""),row.names=F)
png(paste("simu0519-",outputName,".png",sep=""),width=10,height=10,units="in",res=600)
plot(results[,1],results[,2],pch=".",cex=3,col="red",xlab="Lambda",ylab="Element-wise difference",main=paste(paste("ep1=",ep1,sep=""),paste("ep2=",ep2,sep=""),sep=","))
dev.off()


}



###########################################################################################
#######simulation
simul(0.1,0.1,0)
simul(0.2,0.2,1)
simul(0.05,0.05,2)
simul(0.2,0.05,3)
simul(0.05,0.2,4)

simul(0.1,0.05,5)
simul(0.05,0.1,6)
simul(0.2,0.1,7)
simul(0.1,0.2,8)

###########################################################################################
#######plot

ep1=c(0.1,0.2,0.05,0.2,0.05,0.1,0.05,0.2,0.1)
ep2=c(0.1,0.2,0.05,0.05,0.2,0.05,0.1,0.1,0.2)
png("simulation0602.png",width=10,height=10,units="in",res=600)
par(mfrow=c(3,3))
for (i in 1:9){
  j=i-1
results=read.csv(paste("0519results",j,".csv",sep=""))
plot(results[,1],results[,2],pch=".",cex=3,col="red",xlab="Lambda",ylab="Element-wise difference",main=paste(paste("ep1=",ep1[i],sep=""),paste("ep2=",ep2[i],sep=""),sep=","))

}

dev.off()
