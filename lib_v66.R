prior_lp = function(beta,tau){
  pi_beta = - sum(beta^2/tau/2)
  return(sum(pi_beta))
}

#lik_lp <- function(U,lambda){
#  sum(log(1/(1+U%*%lambda)))
#}

lik_lp <- function(U,lambda){
  a = log(1/(1+U%*%lambda))
  if(sum(is.na(a))<3){return(sum(a,na.rm=TRUE))}else{return(sum(a))}
}

penlik_lp <- function(beta,tau,U,lambda){
  lik_lp(U,lambda) + prior_lp(beta,tau)
}

metmove_lp = function(U_old,lambda_old,beta_old,X,Y,tau,itern,p,Sigma_lp,mu_lp){
  out = list()
  out$accept = 0
  out$note = 0

  un = rnorm(1)
  beta_new = mvrnorm(mu=mu_lp,Sigma=Sigma_lp)
  U = X*as.vector(Y-X%*%beta_new)
  dslnex <- function(x) {
    y <- apply(U*as.vector(1/(1+U%*%x)),2,sum)
    y
  }
  xstart <- rep(0,p)
  lambda = nleqslv(xstart, dslnex, control=list(btol=.01))$x  
  
  gnew = dmvnorm(c(beta_new),mean=c(mu_lp),sigma=Sigma_lp,log=TRUE)
  gold = dmvnorm(c(beta_old),mean=c(mu_lp),sigma=Sigma_lp,log=TRUE)  
  fold = penlik_lp(beta_old,tau,U_old,lambda_old)
  fnew = penlik_lp(beta_new,tau,U,lambda) # on the logarithm scale the negative of the penalized likelihood.

  if(sum(abs(dslnex(lambda)))>.1) {fnew = -99999;print("note!!!");out$note=1}
  r = fnew + gold - fold - gnew
  if(!is.na(r)){
    u = runif(1)
    if(exp(r)>u){ # accept
      out$accept = 1
      out$beta = beta_new
      out$weight = fnew
      out$lambda = lambda
      out$U = U
    }
  }
  return(out)
}

bayesemp_lp = function(X,Y,niter = 10000,burnin = 5000,burnin2,start="OLS",b=NULL){
  p = ncol(X)
  n = length(Y)
  gamma = rexp(1,rate = 2)

  Weight = c()
  Beta = matrix(0,niter,p)
  Tau = matrix(0,niter,p)
  r = c()
  sig = c()
  gam = c()
  choice = c()
  note = 0

  gam[1] = gamma

  A = t(X)%*%X/n
  B = (t(X)%*%Y)/n
  beta_ols = solve(A)%*%(B)
  Eps = diag(as.vector(Y-X%*%beta_ols)^2)
  Tn = t(X)%*%Eps%*%X/n
  #Tn = t(X)%*%X/n*mean(diag(Eps))

  if(start=="OLS")  beta = beta_ols
  if(start=="origin") beta = rep(0,p)
  if(start=="provided") beta = b
  
  U = X*as.vector(Y-X%*%beta)
  dslnex <- function(x) {
    y <- apply(U*as.vector(1/(1+U%*%x)),2,sum)
    y
  }
  xstart <- rep(0,p)
  lambda = nleqslv(xstart, dslnex, control=list(btol=.01))$x
  

  Beta[1,] = beta
  Tau[1,] = rep(1,p)
  accept = 1
  
  Weight[1] = penlik_lp(beta,Tau[1,],U,lambda) # on the logarithm scale the negative of the penalized likelihood.

  U_old = U
  lambda_old = lambda
  
  for(iter in 2:niter){
    if(iter%%1000==0) print(iter)
    ## first copy everything to current. 
    ## if accept then update part of it
    
    Beta[iter,] = Beta[iter-1,]
    sig[iter] = sig[iter-1]
    Weight[iter] = Weight[iter-1]
    accept[iter] = 0
    
    uni = runif(1)
    
    # Step 1: update gamma
    
    gam[iter] = sqrt(2*p/sum(Tau[1:(iter-1),])*(iter-1))
    
    # step 2: sample tau
    
    Tau[iter,] = 1/rinvgauss(n=1,m=sqrt(gam[iter]^2/Beta[iter-1,]^2),s=gam[iter]^2)
    
    # step 3: sample beta
    beta_old = Beta[iter-1,]
    Sigma_lp = solve(n*t(A)%*%solve(Tn)%*%A + diag(1/Tau[iter,]))
    mu_lp = t(t(B)%*%solve(Tn)%*%A%*%Sigma_lp)*n

    
    out = metmove_lp(U_old,lambda_old,beta_old,X,Y,Tau[iter,],iter,p,Sigma_lp,mu_lp)
    note = note + out$note
    if(out$accept==1){
      Weight[iter] = out$weight
      Beta[iter,] = out$beta
      accept[iter] = 1
      U_old = out$U
      lambda_old = out$lambda
    }
  }
  
  reg = list()
  
  burn = burnin
  reg$beta = apply(Beta[-(1:burn),],2,mean)
  reg$accept = accept
  reg$tau = apply(Tau[-(1:burn),],2,mean)
  reg$gam = mean(gam[-(1:burn)])
  
  id2 = 1:((niter-burn)/10)*10+burn
  reg$beta2 = apply(Beta[id2,],2,mean)
  reg$tau2 = apply(Tau[id2,],2,mean)
  reg$gam2 = mean(gam[id2])
  
  reg$Beta = Beta
  reg$Weight = Weight
  reg$choice = choice
  reg$gamma = gam
  reg$Tau = Tau
  
  reg$note = note
  return(reg)
}

bayesemp_lp_fixgam = function(X,Y,niter = 10000,burnin = 5000,burnin2,start="OLS",b=NULL,gamma){
  p = ncol(X)
  n = length(Y)
  #gamma = rexp(1,rate = 2)
  
  Weight = c()
  Beta = matrix(0,niter,p)
  Tau = matrix(0,niter,p)
  r = c()
  sig = c()
  gam = rep(gamma,niter)
  choice = c()
  note = 0
  
  gam[1] = gamma
  
  A = t(X)%*%X/n
  B = (t(X)%*%Y)/n
  beta_ols = solve(A)%*%(B)
  Eps = diag(as.vector(Y-X%*%beta_ols)^2)
  Tn = t(X)%*%Eps%*%X/n
  #Tn = t(X)%*%X/n*mean(diag(Eps))
  
  if(start=="OLS")  beta = beta_ols
  if(start=="origin") beta = rep(0,p)
  if(start=="provided") beta = b
  
  U = X*as.vector(Y-X%*%beta)
  dslnex <- function(x) {
    y <- apply(U*as.vector(1/(1+U%*%x)),2,sum)
    y
  }
  xstart <- rep(0,p)
  lambda = nleqslv(xstart, dslnex, control=list(btol=.01))$x
  
  
  Beta[1,] = beta
  Tau[1,] = rep(1,p)
  accept = 1
  
  Weight[1] = penlik_lp(beta,Tau[1,],U,lambda) # on the logarithm scale the negative of the penalized likelihood.
  
  U_old = U
  lambda_old = lambda
  
  for(iter in 2:niter){
    if(iter%%1000==0) print(iter)
    ## first copy everything to current. 
    ## if accept then update part of it
    
    Beta[iter,] = Beta[iter-1,]
    sig[iter] = sig[iter-1]
    Weight[iter] = Weight[iter-1]
    accept[iter] = 0
    
    uni = runif(1)
    
    # Step 1: update gamma
    
    #gam[iter] = sqrt(2*p/sum(Tau[1:(iter-1),])*(iter-1))
    
    # step 2: sample tau
    
    Tau[iter,] = 1/rinvgauss(n=1,m=sqrt(gam[iter]^2/Beta[iter-1,]^2),s=gam[iter]^2)
    
    # step 3: sample beta
    beta_old = Beta[iter-1,]
    Sigma_lp = solve(n*t(A)%*%solve(Tn)%*%A + diag(1/Tau[iter,]))
    mu_lp = t(t(B)%*%solve(Tn)%*%A%*%Sigma_lp)*n
    
    
    out = metmove_lp(U_old,lambda_old,beta_old,X,Y,Tau[iter,],iter,p,Sigma_lp,mu_lp)
    note = note + out$note
    if(out$accept==1){
      Weight[iter] = out$weight
      Beta[iter,] = out$beta
      accept[iter] = 1
      U_old = out$U
      lambda_old = out$lambda
    }
  }
  
  reg = list()
  
  burn = burnin
  reg$beta = apply(Beta[-(1:burn),],2,mean)
  reg$accept = accept
  reg$tau = apply(Tau[-(1:burn),],2,mean)
  #reg$gam = mean(gam[-(1:burn)])
  
  id2 = 1:((niter-burn)/10)*10+burn
  reg$beta2 = apply(Beta[id2,],2,mean)
  reg$tau2 = apply(Tau[id2,],2,mean)
  #reg$gam2 = mean(gam[id2])
  
  reg$Beta = Beta
  reg$Weight = Weight
  reg$choice = choice
  #reg$gamma = gam
  reg$Tau = Tau
  
  reg$note = note
  return(reg)
}

