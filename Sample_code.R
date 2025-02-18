setwd("D:/ycData/Dropbox/research/PROJECTS/EL/Bayes Variable Selection for Emp")
source("./Rcode/lib_v66.R") ## approximated normal proposal
require("nleqslv")
require("mvtnorm")
require("monomvn")
require(rmutil)

p = 10
n = 100

## generate data
x = matrix(rnorm(n*p),n)
b = matrix(0,p,1)
b[1:5] = 1:5
y = x %*% b + rnorm(n)*3

## center data
X = x
mu = apply(X,2,mean)
X = t(t(X) - mu)
Y = y-mean(y)

res2 = bayesemp_lp(X,Y,niter = 1000,burnin = 200,burnin2 = 500)
res2$beta


