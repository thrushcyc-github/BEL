setwd("D:/ycData/Dropbox/research/PROJECTS/EL/Bayes Variable Selection for Emp")
source("./Rcode/lib_v66.R") ## approximated normal proposal
require("nleqslv")
require("mvtnorm")
require("monomvn")
require(rmutil)

source("./Rcode/lib_v66.R") ## approximated normal proposal
source("./Rcode/lib_scad.R")

p = 10
n = 100

# data generation
x = matrix(rnorm(n*p),n)
b = matrix(0,p,1)
b[1:5] = 1:5
y = x %*% b + rnorm(n)*3

#############################
###
### BEL_lasso
###
################################

## cross validation 

err = matrix(0,5,50)
for(fold in 1:5){
  id = sample(n,round(n*0.7))
  X = x[id,]
  Y = y[id]
  res2 = bayesemp_lp(X,Y,niter = 1000,burnin = 200,burnin2 = 200)
  beta = res2$beta
  for(ii in 1:50){
    th = ii/100
    sel = which(abs(beta)>th)
    Xs = X[,sel]
    be = solve(t(Xs)%*%Xs)%*%(t(Xs)%*%Y)
    xnew = x[-id,sel]
    if(length(b)==1){ pred = xnew*be }else{ pred = xnew%*%be }
    err[fold,ii] = sum((y[-id]-pred)^2)
  }
}

# choose the cut-off
th = which.min(apply(err,2,mean))/100

## final run using the whole data and the cut-off
X = x
mu = apply(X,2,mean)
X = t(t(X) - mu)
Y = y-mean(y)
res2 = bayesemp_lp(X,Y,niter = 1000,burnin = 200,burnin2 = 500)
beta = res2$beta*(abs(res2$beta)>th)
