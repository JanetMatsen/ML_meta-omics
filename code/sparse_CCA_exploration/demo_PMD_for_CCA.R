# needed bioconductor
# https://www.bioconductor.org/install/
# library(BiocInstaller)  # https://support.bioconductor.org/p/58339/
# biocLite("impute")

library("PMA")
library("ggplot2")


# first, do CCA with type="standard"
# A simple simulated example
u <- matrix(c(rep(1,25),rep(0,75)),ncol=1)
v1 <- matrix(c(rep(1,50),rep(0,450)),ncol=1)
v2 <- matrix(c(rep(0,50),rep(1,50),rep(0,900)),ncol=1)
x <- u%*%t(v1) + matrix(rnorm(100*500),ncol=500)
z <- u%*%t(v2) + matrix(rnorm(100*1000),ncol=1000)
# Can run CCA with default settings, and can get e.g. 3 components
out <- CCA(x,z,typex="standard",typez="standard",K=3)
print(out,verbose=TRUE) # To get less output, just print(out)
# Or can use CCA.permute to choose optimal parameter values
perm.out <- CCA.permute(x,z,typex="standard",typez="standard",nperms=7)
print(perm.out)
plot(perm.out)
out <- CCA(x,z,typex="standard",typez="standard",K=1,penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init)
print(out)