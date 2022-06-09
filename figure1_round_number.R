library(pqr)
library(MASS)
##############################计算beta eta程序
Pi.n = function(u,x,tau)
{
  weights=tau-ifelse(u<=tau,1,0)
  tmp=x*weights
  Pi.value=max(abs(apply(tmp,2,sum)))
  return(Pi.value)
}

penalty = function(x,tau,level, sim.n)
{
  n=dim(x)[1];
  U=matrix(runif(n*sim.n),n,sim.n) 
  values = apply(U,2,Pi.n,x=x,tau=tau)
  values = sort(values)
  return(values[ceiling(level*sim.n)])
}

ini.est = function(x,y,index,tau,level=0.85,lam=NULL,sim.n=1000){
  if(is.null(lam)){
    lambda=penalty(x[,-index],tau,level=level,sim.n=sim.n)
  }else{lambda=lam}
  p=dim(x)[2]
  lam.vec=numeric(p)
  lam.vec[-index]=lambda
  result = rq.fit.lasso(x,y,tau=tau,lambda=lam.vec)
  return(list(residual=result$residual,coefficient=result$coefficient))
}

kern=function(v){     # smoothing function
  return(pnorm(v,0,1)) 
}




###################
####################经验似然

llog <- function (z, eps)
{
  ans <- z
  lo <- (z < eps)
  ans[lo] <- log(eps) - 1.5 + 2 * z[lo]/eps - 0.5 * (z[lo]/eps)^2
  ans[!lo] <- log(z[!lo])
  ans
}

llogp <- function (z, eps)
{
  ans <- z
  lo <- (z < eps)
  ans[lo] <- 2/eps - z[lo]/eps^2
  ans[!lo] <- 1/z[!lo]
  ans
}

llogpp <- function (z, eps)
{
  ans <- z
  lo <- (z < eps)
  ans[lo] <- -1/eps^2
  ans[!lo] <- -1/z[!lo]^2
  ans
}

el.ee <- function (theta,gee,lam, maxit=25, gradtol=1e-7, svdtol=1e-9,
                   itertrace=F, TINY=sqrt(.Machine$double.xmin))
{
  mx=gee(theta);
  mx <- as.matrix(mx)
  n <- nrow(mx)
  p <- ncol(mx)
  
  if (n <= p)
    stop("Need more observations than length(mu) in el.test().")
  z <- mx
  scale <- mean(abs(z)) + TINY
  z <- z/scale
  if (!missing(lam)) {
    lam <- as.vector(lam)
    lam <- lam * scale
    if (logelr.ee(z, lam) > 0)
      lam <- rep(0, p)
  }
  if (missing(lam))
    lam <- rep(0, p)
  if (svdtol < TINY)
    svdtol <- TINY
  if (gradtol < TINY)
    gradtol <- TINY
  nwts <- c(3^-c(0:3), rep(0, 12))
  gwts <- 2^(-c(0:(length(nwts) - 1)))
  gwts <- (gwts^2 - nwts^2)^0.5
  gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
  nits <- 0
  gsize <- gradtol + 1
  while (nits < maxit && gsize > gradtol) {
    arg <- 1 + z %*% lam
    wts1 <- as.vector(llogp(arg, 1/(n)))
    wts2 <- as.vector(-llogpp(arg, 1/(n)))^0.5
    grad <- as.matrix(-z * wts1)
    grad <- as.vector(rowsum(grad, rep(1, nrow(grad))))
    gsize <- mean(abs(grad))
    hess <- z * wts2
    svdh <- svd(hess)
    if (min(svdh$d) < max(svdh$d) * svdtol)
      svdh$d <- svdh$d + max(svdh$d) * svdtol
    nstep <- svdh$v %*% (t(svdh$u)/svdh$d)
    nstep <- as.vector(nstep %*% matrix(wts1/wts2, n, 1))
    gstep <- -grad
    if (sum(nstep^2) < sum(gstep^2))
      gstep <- gstep * (sum(nstep^2)^0.5/sum(gstep^2)^0.5)
    ologelr <- -sum(llog(arg, 1/(n)))
    ninner <- 0
    for (i in 1:length(nwts)) {
      nlogelr <- logelr.ee(z, lam + nwts[i] * nstep + gwts[i] * gstep)
      if (nlogelr < ologelr) {
        lam <- lam + nwts[i] * nstep + gwts[i] * gstep
        ninner <- i
        break
      }
    }
    nits <- nits + 1
    if (ninner == 0)
      nits <- maxit
    if (itertrace)
      print(c(lam, nlogelr, gsize, ninner))
  }
  list(`-2LLR` = -2 * nlogelr, Pval = 1 - pchisq(-2 * nlogelr, df = p),
       lambda = lam/scale, grad = grad * scale, hess = t(hess) %*%
         hess * scale^2, wts = wts1, nits = nits)
  
  return(-2 * nlogelr)  
}

logelr.ee <- function (mx, lam)
{
  mx <- as.matrix(mx)
  n <- nrow(mx)
  p <- ncol(mx)
  if (n <= p)
    stop("Need more observations than variables in logelr.")
  z <- mx
  arg <- 1 + z %*% lam
  return(-sum(llog(arg, 1/n)))
}

#########################################################################################
n <- 150
M <- 10
N <- n*M
d <- 5  #the column of x
s <- 4   #z
p <- 250 #the column of z
P=d+p-1
Tau <- 0.5
Run <- 50
beta0 <- matrix(1,d,1)
eta0 <- c(3,1.5,1,2,numeric(p-s))


Sigma=diag(1,P,P)
rr=0.1;
for (i in 1:P){
  for (j in 1:P)
  {
    Sigma[i,j]=(rr^(abs(i-j)))/4
  }
}

beta_all=matrix(0,d,M)
eta_all=matrix(0,p,M)
H_hat=matrix(0,(d*M),p)

beta_hat1=matrix(0,d,M)
beta_hat8=matrix(0,d,M)
beta_hat_8=matrix(0,d,M)
beta_hat9=matrix(0,d,M)
beta_hat_9=matrix(0,d,M)
result=matrix(0,(p+d),M)


beta_pqr_1=matrix(0,(d*Run),1)
beta_pqr_2=matrix(0,(d*Run),1)
beta_pqr_3=matrix(0,(d*Run),1)
beta_pqr_4=matrix(0,(d*Run),1)
beta_pqr_5=matrix(0,(d*Run),1)

beta_pqr_1_all=matrix(0,(d*Run),1)
beta_pqr_2_all=matrix(0,(d*Run),1)
beta_pqr_3_all=matrix(0,(d*Run),1)
beta_pqr_4_all=matrix(0,(d*Run),1)
beta_pqr_5_all=matrix(0,(d*Run),1)


beta00=matrix(0,(d*Run),1)
X=matrix(0,N,P)
x1=matrix(0,n,d)
z1=matrix(0,n,(p-1))
y1=matrix(0,n,1)
eta1=matrix(0,p,1)
H1=matrix(0,d,p)

for (r in 1:Run){
  
  X=mvrnorm(N, rep(0,P), Sigma)
  x <- matrix(X[,1:d],N,d)
  z <- matrix(X[,(d+1):P],N,(p-1))
  y <- x%*%beta0 + cbind(1,z)%*%eta0 +rt(N,3)#*(1+0.2*x[,1])#rchisq(N,2)
  ######################################################################################average方法
  for (m in 1:M){
    result8 = ini.est(cbind(x[((m-1)*n+1):(m*n),],1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),],index=c(1:(d+1)),tau=Tau)  
    beta_all[,m]=result8$coefficient[1:d]
    eta_all[,m]=result8$coefficient[(d+1):(d+p)]                        #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
    fit = mvr(x[((m-1)*n+1):(m*n),],z[((m-1)*n+1):(m*n),])
    H_hat[((m-1)*d+1):(m*d),] = t(rbind(as.vector(fit$muhat),fit$Bhat))  #每一个机器算一个H
  }
  beta_hat=apply(beta_all,1,mean)
  eta_hat=apply(eta_all,1,mean)  
  beta00[((r-1)*d+1):(r*d),1]=beta_hat
  ######################################################################################################################################################
  H_csl=matrix(0,d,p)
  for (m in 1:M){
    H_csl=H_csl+H_hat[((m-1)*d+1):(m*d),]/M
  } 
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  res2=matrix(0,n,1)
  L_1=0
  for (i in 1:n)
  {
    res2[i,]=y[i,]-x[i,]%*%beta_hat-cbind(1,z)[i,]%*%eta_hat
  }
  h2=1.5*sd(c(res2))*n^(-1/3)
  for (i in 1:n){
    L_1=L_1+c( Tau-kern( (-y[i,]+x[i,]%*%beta_hat+cbind(1,z)[i,]%*%eta_hat)/h2 ) )*( matrix(x[i,])- H_csl%*%cbind(1,z)[i,] )
  }
  L_1=L_1/n
  
  beta_hat_4=beta_hat
  x1=x[1:n,]
  z1=z[1:n,]
  y1=matrix(y[1:n,])
  g_pqr<-function(beta)   
  {
    res=matrix(0,n,1)
    for (i in 1:n)
    {
      res[i,]=y1[i,]-x1[i,]%*%beta_hat_4-cbind(1,z1)[i,]%*%eta_hat
    }
    h=1.5*sd(c(res))*n^(-1/3)
    re=matrix(0,n,d)
    beta=as.matrix(beta)
    for (i in 1:n)
    {
      re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_1 )
    }
    return(re)
  }
  f_pqr=function(beta)
  {
    f=el.ee(beta,gee=g_pqr)
    return(f)    
  }
  beta_hat4_1=optim(beta_hat_4,f_pqr,method = "BFGS")$par
  beta_pqr_1[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #########################################################################################################################
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat4_1,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat_1=apply(eta_all,1,mean)  
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat4_1-cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat4_1+cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  res2=matrix(0,n,1)
  L_1=0
  for (i in 1:n)
  {
    res2[i,]=y[i,]-x[i,]%*%beta_hat4_1-cbind(1,z)[i,]%*%eta_hat_1
  }
  h2=1.5*sd(c(res2))*n^(-1/3)
  for (i in 1:n){
    L_1=L_1+c( Tau-kern( (-y[i,]+x[i,]%*%beta_hat4_1+cbind(1,z)[i,]%*%eta_hat_1)/h2 ) )*( matrix(x[i,])- H_csl%*%cbind(1,z)[i,] )
  }
  L_1=L_1/n
  
  g_pqr<-function(beta)   
  {
    res=matrix(0,n,1)
    for (i in 1:n)
    {
      res[i,]=y1[i,]-x1[i,]%*%beta_hat4_1-cbind(1,z1)[i,]%*%eta_hat_1
    }
    h=1.5*sd(c(res))*n^(-1/3)
    re=matrix(0,n,d)
    beta=as.matrix(beta)
    for (i in 1:n)
    {
      re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat_1)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_1 )
    }
    return(re)
  }
  f_pqr=function(beta)
  {
    f=el.ee(beta,gee=g_pqr)
    return(f)    
  }
  beta_hat4_1=optim(beta_hat4_1,f_pqr,method = "BFGS")$par
  beta_pqr_2[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #####################################################################################################
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat4_1,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat_1=apply(eta_all,1,mean)  
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat4_1-cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat4_1+cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  res2=matrix(0,n,1)
  L_1=0
  for (i in 1:n)
  {
    res2[i,]=y[i,]-x[i,]%*%beta_hat4_1-cbind(1,z)[i,]%*%eta_hat_1
  }
  h2=1.5*sd(c(res2))*n^(-1/3)
  for (i in 1:n){
    L_1=L_1+c( Tau-kern( (-y[i,]+x[i,]%*%beta_hat4_1+cbind(1,z)[i,]%*%eta_hat_1)/h2 ) )*( matrix(x[i,])- H_csl%*%cbind(1,z)[i,] )
  }
  L_1=L_1/n
  
  g_pqr<-function(beta)   
  {
    res=matrix(0,n,1)
    for (i in 1:n)
    {
      res[i,]=y1[i,]-x1[i,]%*%beta_hat4_1-cbind(1,z1)[i,]%*%eta_hat_1
    }
    h=1.5*sd(c(res))*n^(-1/3)
    re=matrix(0,n,d)
    beta=as.matrix(beta)
    for (i in 1:n)
    {
      re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat_1)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_1 )
    }
    return(re)
  }
  f_pqr=function(beta)
  {
    f=el.ee(beta,gee=g_pqr)
    return(f)    
  }
  beta_hat4_1=optim(beta_hat4_1,f_pqr,method = "BFGS")$par
  beta_pqr_3[((r-1)*d+1):(r*d),1]=beta_hat4_1
  ############################################################################
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat4_1,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat_1=apply(eta_all,1,mean)  
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat4_1-cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat4_1+cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  res2=matrix(0,n,1)
  L_1=0
  for (i in 1:n)
  {
    res2[i,]=y[i,]-x[i,]%*%beta_hat4_1-cbind(1,z)[i,]%*%eta_hat_1
  }
  h2=1.5*sd(c(res2))*n^(-1/3)
  for (i in 1:n){
    L_1=L_1+c( Tau-kern( (-y[i,]+x[i,]%*%beta_hat4_1+cbind(1,z)[i,]%*%eta_hat_1)/h2 ) )*( matrix(x[i,])- H_csl%*%cbind(1,z)[i,] )
  }
  L_1=L_1/n
  
  g_pqr<-function(beta)   
  {
    res=matrix(0,n,1)
    for (i in 1:n)
    {
      res[i,]=y1[i,]-x1[i,]%*%beta_hat4_1-cbind(1,z1)[i,]%*%eta_hat_1
    }
    h=1.5*sd(c(res))*n^(-1/3)
    re=matrix(0,n,d)
    beta=as.matrix(beta)
    for (i in 1:n)
    {
      re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat_1)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_1 )
    }
    return(re)
  }
  f_pqr=function(beta)
  {
    f=el.ee(beta,gee=g_pqr)
    return(f)    
  }
  beta_hat4_1=optim(beta_hat4_1,f_pqr,method = "BFGS")$par
  beta_pqr_4[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #################################################################################
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat4_1,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat_1=apply(eta_all,1,mean)  
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat4_1-cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat4_1+cbind(1,z)[((m-1)*n+i),]%*%eta_hat_1)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  res2=matrix(0,n,1)
  L_1=0
  for (i in 1:n)
  {
    res2[i,]=y[i,]-x[i,]%*%beta_hat4_1-cbind(1,z)[i,]%*%eta_hat_1
  }
  h2=1.5*sd(c(res2))*n^(-1/3)
  for (i in 1:n){
    L_1=L_1+c( Tau-kern( (-y[i,]+x[i,]%*%beta_hat4_1+cbind(1,z)[i,]%*%eta_hat_1)/h2 ) )*( matrix(x[i,])- H_csl%*%cbind(1,z)[i,] )
  }
  L_1=L_1/n
  
  g_pqr<-function(beta)   
  {
    res=matrix(0,n,1)
    for (i in 1:n)
    {
      res[i,]=y1[i,]-x1[i,]%*%beta_hat4_1-cbind(1,z1)[i,]%*%eta_hat_1
    }
    h=1.5*sd(c(res))*n^(-1/3)
    re=matrix(0,n,d)
    beta=as.matrix(beta)
    for (i in 1:n)
    {
      re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat_1)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_1 )
    }
    return(re)
  }
  f_pqr=function(beta)
  {
    f=el.ee(beta,gee=g_pqr)
    return(f)    
  }
  beta_hat4_1=optim(beta_hat4_1,f_pqr,method = "BFGS")$par
  beta_pqr_5[((r-1)*d+1):(r*d),1]=beta_hat4_1
  ##################################################################################################################
  
  beta_hat_9=beta_hat
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  for (m in 1:M){
    res3=matrix(0,n,1)
    L_m=0
    for (i in 1:n)
    {
      res3[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h3=1.5*sd(c(res3))*n^(-1/3)
    for (i in 1:n){
      L_m=L_m+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h3 ) )*( matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
      
    }
    L_m=L_m/n
    
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    g_pqr_ave<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_9-cbind(1,z1)[i,]%*%eta_hat
      }
      h=1.5*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_m )
      }
      return(re)
    }
    f_pqr_ave=function(beta)
    {
      f=el.ee(beta,gee=g_pqr_ave)
      return(f)    
    }
    beta_hat9[,m]=optim(beta_hat_9,f_pqr_ave,method = "BFGS")$par
  }
  beta_hat_9=apply(beta_hat9,1,mean)
  beta_pqr_1_all[((r-1)*d+1):(r*d),1]=beta_hat_9
  
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat_9,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat=apply(eta_all,1,mean) 
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  for (m in 1:M){
    res3=matrix(0,n,1)
    L_m=0
    for (i in 1:n)
    {
      res3[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h3=1.5*sd(c(res3))*n^(-1/3)
    for (i in 1:n){
      L_m=L_m+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h3 ) )*( matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
      
    }
    L_m=L_m/n
    
    
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    g_pqr_ave<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_9-cbind(1,z1)[i,]%*%eta_hat
      }
      h=1.5*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_m )
      }
      return(re)
    }
    f_pqr_ave=function(beta)
    {
      f=el.ee(beta,gee=g_pqr_ave)
      return(f)    
    }
    beta_hat9[,m]=optim(beta_hat_9,f_pqr_ave,method = "BFGS")$par
  }
  beta_hat_9=apply(beta_hat9,1,mean)
  beta_pqr_2_all[((r-1)*d+1):(r*d),1]=beta_hat_9
  
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat_9,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat=apply(eta_all,1,mean) 
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  for (m in 1:M){
    res3=matrix(0,n,1)
    L_m=0
    for (i in 1:n)
    {
      res3[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h3=1.5*sd(c(res3))*n^(-1/3)
    for (i in 1:n){
      L_m=L_m+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h3 ) )*( matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
      
    }
    L_m=L_m/n
    
    
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    g_pqr_ave<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_9-cbind(1,z1)[i,]%*%eta_hat
      }
      h=1.5*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_m )
      }
      return(re)
    }
    f_pqr_ave=function(beta)
    {
      f=el.ee(beta,gee=g_pqr_ave)
      return(f)    
    }
    beta_hat9[,m]=optim(beta_hat_9,f_pqr_ave,method = "BFGS")$par
  }
  beta_hat_9=apply(beta_hat9,1,mean)
  beta_pqr_3_all[((r-1)*d+1):(r*d),1]=beta_hat_9
  
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat_9,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat=apply(eta_all,1,mean) 
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  for (m in 1:M){
    res3=matrix(0,n,1)
    L_m=0
    for (i in 1:n)
    {
      res3[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h3=1.5*sd(c(res3))*n^(-1/3)
    for (i in 1:n){
      L_m=L_m+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h3 ) )*( matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
      
    }
    L_m=L_m/n
    
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    g_pqr_ave<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_9-cbind(1,z1)[i,]%*%eta_hat
      }
      h=1.5*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_m )
      }
      return(re)
    }
    f_pqr_ave=function(beta)
    {
      f=el.ee(beta,gee=g_pqr_ave)
      return(f)    
    }
    beta_hat9[,m]=optim(beta_hat_9,f_pqr_ave,method = "BFGS")$par
  }
  beta_hat_9=apply(beta_hat9,1,mean)
  beta_pqr_4_all[((r-1)*d+1):(r*d),1]=beta_hat_9
  
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%beta_hat_9,index=c(1),tau=Tau)  
    eta_all[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all[j,m])<0.2){
        eta_all[j,m]=0
      }
    }
  }
  eta_hat=apply(eta_all,1,mean) 
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-1/3)
    for (i in 1:n){
      L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    }
  }
  L_N=L_N/N
  
  for (m in 1:M){
    res3=matrix(0,n,1)
    L_m=0
    for (i in 1:n)
    {
      res3[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h3=1.5*sd(c(res3))*n^(-1/3)
    for (i in 1:n){
      L_m=L_m+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat)/h3 ) )*( matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
      
    }
    L_m=L_m/n
    
    
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    g_pqr_ave<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_9-cbind(1,z1)[i,]%*%eta_hat
      }
      h=1.5*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_m )
      }
      return(re)
    }
    f_pqr_ave=function(beta)
    {
      f=el.ee(beta,gee=g_pqr_ave)
      return(f)    
    }
    beta_hat9[,m]=optim(beta_hat_9,f_pqr_ave,method = "BFGS")$par
  }
  beta_hat_9=apply(beta_hat9,1,mean)
  beta_pqr_5_all[((r-1)*d+1):(r*d),1]=beta_hat_9
  print(r)
}
msepqr1_all=rep(0,Run)
madpqr1_all=rep(0,Run)
msepqr2_all=rep(0,Run)
madpqr2_all=rep(0,Run)
msepqr3_all=rep(0,Run)
madpqr3_all=rep(0,Run)
msepqr4_all=rep(0,Run)
madpqr4_all=rep(0,Run)
msepqr5_all=rep(0,Run)
madpqr5_all=rep(0,Run)

msepqr1_all_all=rep(0,Run)
madpqr1_all_all=rep(0,Run)
msepqr2_all_all=rep(0,Run)
madpqr2_all_all=rep(0,Run)
msepqr3_all_all=rep(0,Run)
madpqr3_all_all=rep(0,Run)
msepqr4_all_all=rep(0,Run)
madpqr4_all_all=rep(0,Run)
msepqr5_all_all=rep(0,Run)
madpqr5_all_all=rep(0,Run)
mse00_all=rep(0,Run)
mad00_all=rep(0,Run)



#gamma=matrix(0,d,1)
#gamma[1]=0.2
#beta0=beta0+gamma*qt(Tau,2)
for (r in 1:Run){
  mse00_all[r]=sum( (beta00[((r-1)*d+1):(r*d),1]-beta0)^2 )
  mad00_all[r]=sum( abs(beta00[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr1_all[r]=sum( (beta_pqr_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr1_all[r]=sum( abs(beta_pqr_1[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr2_all[r]=sum( (beta_pqr_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr2_all[r]=sum( abs(beta_pqr_2[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr3_all[r]=sum( (beta_pqr_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr3_all[r]=sum( abs(beta_pqr_3[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr4_all[r]=sum( (beta_pqr_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr4_all[r]=sum( abs(beta_pqr_4[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr5_all[r]=sum( (beta_pqr_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr5_all[r]=sum( abs(beta_pqr_5[((r-1)*d+1):(r*d),1]-beta0) )
  
  msepqr1_all_all[r]=sum( (beta_pqr_1_all[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr1_all_all[r]=sum( abs(beta_pqr_1_all[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr2_all_all[r]=sum( (beta_pqr_2_all[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr2_all_all[r]=sum( abs(beta_pqr_2_all[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr3_all_all[r]=sum( (beta_pqr_3_all[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr3_all_all[r]=sum( abs(beta_pqr_3_all[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr4_all_all[r]=sum( (beta_pqr_4_all[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr4_all_all[r]=sum( abs(beta_pqr_4_all[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr5_all_all[r]=sum( (beta_pqr_5_all[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr5_all_all[r]=sum( abs(beta_pqr_5_all[((r-1)*d+1):(r*d),1]-beta0) )
}


msepqr1=mean(msepqr1_all)
madpqr1=mean(madpqr1_all)
msepqr2=mean(msepqr2_all)
madpqr2=mean(madpqr2_all)
msepqr3=mean(msepqr3_all)
madpqr3=mean(madpqr3_all)
msepqr4=mean(msepqr4_all)
madpqr4=mean(madpqr4_all)
msepqr5=mean(msepqr5_all)
madpqr5=mean(madpqr5_all)

msepqr1_full=mean(msepqr1_all_all)
madpqr1_full=mean(madpqr1_all_all)
msepqr2_full=mean(msepqr2_all_all)
madpqr2_full=mean(madpqr2_all_all)
msepqr3_full=mean(msepqr3_all_all)
madpqr3_full=mean(madpqr3_all_all)
msepqr4_full=mean(msepqr4_all_all)
madpqr4_full=mean(madpqr4_all_all)
msepqr5_full=mean(msepqr5_all_all)
madpqr5_full=mean(madpqr5_all_all)
mse00=mean(mse00_all)
mad00=mean(mad00_all)
mse=c(
  mse00,
  msepqr1,
  msepqr2,
  msepqr3,
  msepqr4,
  msepqr5
)
mad=c(
  mad00,
  madpqr1,
  madpqr2,
  madpqr3,
  madpqr4,
  madpqr5
)
results=rbind(mse,mad)
rname = c("mse","mad")
cname = c("pqr0","pqr1","pqr2","pqr3","pqr4","pqr5")


mse_all=c(
  mse00,
  msepqr1_full,
  msepqr2_full,
  msepqr3_full,
  msepqr4_full,
  msepqr5_full
)
mad_all=c(
  mad00,
  madpqr1_full,
  madpqr2_full,
  madpqr3_full,
  madpqr4_full,
  madpqr5_full
)
results_all=rbind(mse_all,mad_all)
rname_all = c("mse","mad")
cname_all = c("pqr0","pqr1_all","pqr2_all","pqr3_all","pqr4_all","pqr5_all")

matrix(results,2,dimnames = list(rname,cname))
matrix(results_all,2,dimnames = list(rname_all,cname_all))
boxplot(mse00_all,msepqr1_all,msepqr2_all,msepqr1_all_all,msepqr2_all_all)
