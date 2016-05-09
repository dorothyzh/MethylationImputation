library(HMP)
library(ggplot2)
library(vegan)
library(micropower)

setwd("/Users/dzhi/Dropbox/microbiome_simulation")
data(saliva)
n=800
mius = c(0.0001, 0.001, 0.01, 0.05, 0.3)
# mius = c(0.0001, 0.3)

big_table = list()

# beta = 0.5; # beta of individual OTUs may have different effect!
depth = 100000
# depth = 10000
power_iter = 100;
iter=19
effect_resolution = 0.01;
effect_resolution = 0.05;
# effect_resolution = 0.3;
max_effect = 0.3;
# max_effect = 0.2;

fit.saliva <- dirmult(saliva)
gamma = sort(fit.saliva$gamma); #fit.saliva$gamma
gamma=exp(seq(log(50),log(0.0000000001), log(0.95)))
gamma = gamma[1:500]
gamma = sort(gamma)
alpha = 0.05;

for (cumsum_cutoff in mius) {
  ii=0;
  power_table = list();
  set.seed(1)
  start <- Sys.time ()
  
  for (beta in seq(0,max_effect,effect_resolution)) {
    pp=list()
    for (pi in 1:power_iter) {
      ii=ii+1
      Nrs <- rep(depth, n)
      otu <- Dirichlet.multinomial(Nrs, gamma)
      otu = cbind(matrix(0,dim(otu)[1],1), otu)
      # heatmap(otu)
      
      # k is the smallest number such that top k OTUs have total abundance above
      # cumsum_cutoff
      k=sum(cumsum(gamma)/sum(gamma)<cumsum_cutoff)+1
      # eps = nnorm
      x=scale(rowSums(otu[,c(1:k)]))
      y= x* beta + rnorm(n=n, m=0, sd=1)
      
      data = data.frame(cbind(y,otu))
      
      fit.adonis=adonis(otu ~ y, data=data, permutations=iter) # bray-curtis by default
      p=fit.adonis$aov.tab[1,6]
      pp=cbind(pp,p)
      
    }
    power = sum(pp<=alpha) / power_iter
    power_table = cbind(power_table, power)
  }
  
  Sys.time () - start
  
  final_table = data.frame(beta=seq(0,max_effect,effect_resolution), power=as.numeric(power_table))
  
  write.csv(final_table, file=sprintf("power.r%d.freq%.3f.table.csv", depth, cumsum_cutoff))
  
  big_table = cbind(big_table, final_table[,2]);
  
}

write.csv(big_table, file=sprintf("power.big.table.csv"))

# plot(seq(0,max_effect,effect_resolution),power_table,, main = "simulated power", xlab="beta", ylab="power")
