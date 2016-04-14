
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

seq<-B[sample(1:n, n.seq), ]
array<-B[, sample(1:m, m.array)]
qc_frac = 1-5/dim(array)[2]
Seq <- MetIm(sequence = t(seq), microarray = t(array), lambda=0.3, cut=10, cvfold=0, use.mvalue = T, qc_frac = qc_frac)

Accuracy = compare(Seq, B)

# plot results


