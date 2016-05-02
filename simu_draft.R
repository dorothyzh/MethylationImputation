# M=log2(Beta/(1-Beta))

# library(MVN)
library(MASS)
library(wateRmelon)
library(MetImpute)

n = 1000
n.seq = 100
m = 200
m.array = 20 
Epsilon1=.1
Epsilon2=.1
var_a = 0.5 #variance 0.5.
s=1.05
t=0.05



nmz = n #number of CpGs, here e.g. 1000
M <- mvrnorm(m, mu = rep(0, nmz), Sigma = diag(nmz)*var_a)  ####200samples

M <- data.frame(M)
colnames(M)=as.character(1:n)

seq<-M[order(sample(1:n, n.seq)), ] ###sequencing
dim(seq) # 100 1000

array<-M[, order(sample(1:m, m.array))]  ###array    200中20 90% 
dim(array) #200  20


####sequencing
noise1=rnorm(dim(seq)[1]*dim(seq)[2], mean = 0, sd = 1) ##n observations
seq1=seq+Epsilon1*noise1   #####add random noise

par(mfrow=c(2,1))

seq2=m2beta(seq1)     #####convert M-value to beta value

####microarray
noise2=rnorm(dim(array)[1]*dim(array)[2], mean = 0, sd = 1) ##n observations
array1=array+Epsilon2*noise2

array2=m2beta(array1) ####m-value to beta value

array3=s*array2+t  ###make the scale of array the same with sequencing, 0.9<s<1.1, -0.1<b<0.1

#####see results
dat1=as.data.frame(seq2)
dim(dat1)
dat2=as.data.frame(array3)
dim(dat2)
qc_frac = 1-5/dim(dat2)[2]

seq_impute <- MetIm(sequence = dat1, microarray = dat2, lambda=0.3, cut=10, cvfold=0, use.mvalue = F, qc_frac = qc_frac)

dim(seq_impute)

actual=seq2
predicted=seq_impute
se(actual, predicted)  #####compare(M,M2) ####using element-wise squared error.

# element-wise difference
diff <- seq_impute-as.matrix(M)
total.se <- sum(sum(diff^2))

