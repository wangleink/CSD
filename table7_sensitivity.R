library(pqr)
library(MASS)
##############################计算beta eta程序
######################################################################################半惩罚算初值
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
#########################################################################################################
kern=function(v){     # smoothing function
  return(pnorm(v,0,1)) 
}
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

#########################################################################################方法比较
n <- 150
M <- 5
N <- n*M
d <- 5  #the column of x
s <- 4   #z
p <- 250 #the column of z
P=d+p-1
Tau <- 0.5
Run <- 20
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


beta_pqr_1_1=matrix(0,(d*Run),1)
beta_pqr_2_1=matrix(0,(d*Run),1)
beta_pqr_3_1=matrix(0,(d*Run),1)
beta_pqr_ave_1=matrix(0,(d*Run),1)
beta_sub_pqr_1=matrix(0,(d*Run),1)
beta_ave_pqr_1=matrix(0,(d*Run),1)
beta_full_pqr_1=matrix(0,(d*Run),1)

beta_pqr_1_2=matrix(0,(d*Run),1)
beta_pqr_2_2=matrix(0,(d*Run),1)
beta_pqr_3_2=matrix(0,(d*Run),1)
beta_pqr_ave_2=matrix(0,(d*Run),1)
beta_sub_pqr_2=matrix(0,(d*Run),1)
beta_ave_pqr_2=matrix(0,(d*Run),1)
beta_full_pqr_2=matrix(0,(d*Run),1)

beta_pqr_1_3=matrix(0,(d*Run),1)
beta_pqr_2_3=matrix(0,(d*Run),1)
beta_pqr_3_3=matrix(0,(d*Run),1)
beta_pqr_ave_3=matrix(0,(d*Run),1)
beta_sub_pqr_3=matrix(0,(d*Run),1)
beta_ave_pqr_3=matrix(0,(d*Run),1)
beta_full_pqr_3=matrix(0,(d*Run),1)

beta_pqr_1_4=matrix(0,(d*Run),1)
beta_pqr_2_4=matrix(0,(d*Run),1)
beta_pqr_3_4=matrix(0,(d*Run),1)
beta_pqr_ave_4=matrix(0,(d*Run),1)
beta_sub_pqr_4=matrix(0,(d*Run),1)
beta_ave_pqr_4=matrix(0,(d*Run),1)
beta_full_pqr_4=matrix(0,(d*Run),1)

beta_pqr_1_5=matrix(0,(d*Run),1)
beta_pqr_2_5=matrix(0,(d*Run),1)
beta_pqr_3_5=matrix(0,(d*Run),1)
beta_pqr_ave_5=matrix(0,(d*Run),1)
beta_sub_pqr_5=matrix(0,(d*Run),1)
beta_ave_pqr_5=matrix(0,(d*Run),1)
beta_full_pqr_5=matrix(0,(d*Run),1)

beta00=matrix(0,(d*Run),1)
beta11=matrix(0,(d*Run),1)
beta22=matrix(0,(d*Run),1)

X=matrix(0,N,P)
y_fullpqr=matrix(0,N,1)
x_fullpqr=matrix(0,N,d)
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
  
  for (m in 1:M){
    beta_hat_8[,m]=beta_all[,m]
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    eta1=matrix(eta_all[,m])
    H1=H_hat[((m-1)*d+1):(m*d),]
    g_sub<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_8[,m]-cbind(1,z1)[i,]%*%eta_all[,m]
      }
      h=1.5*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t(c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta1)/h ))*(matrix(x1[i,])-H1%*%cbind(1,z1)[i,]))
      }
      return(re)
    }
    f_sub=function(beta)
    {
      f=el.ee(beta,gee=g_sub)
      return(f)    
    }
    beta_hat8[,m]=optim(beta_hat_8[,m],f_sub,method = "BFGS")$par
  }
  beta_ave_pqr_1[((r-1)*d+1):(r*d),1]=apply(beta_hat8,1,mean)
  ####################################################################################sub方法
  beta_sub_pqr_1[((r-1)*d+1):(r*d),1]=beta_hat8[,1]
  ####################################################################################full方法
  fit6 = mvr(x,z)
  H_hat6 = t(rbind(as.vector(fit6$muhat),fit6$Bhat))  
  
  result6 = ini.est(cbind(x,1,z),y,index=c(1:(d+1)),tau=Tau)  
  beta_hat_6=result6$coefficient[1:d]
  eta_hat6=result6$coefficient[(d+1):(d+p)]      
  
  for (j in 1:p){
    if (abs(eta_hat6[j])<0.2){
      eta_hat6[j]=0
    }
  }
  #for (i in 1:N){
  #  y_fullpqr[i,]=y[i]-t(H_hat6%*%cbind(1,z)[i,])%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  #  x_fullpqr[i,]=x[i,]-t(H_hat6%*%cbind(1,z)[i,])
  #}
  #beta_hat6=rq(y_fullpqr~x_fullpqr, tau=Tau)$coefficient[2:(d+1)]
  res=matrix(0,N,1)
  for (i in 1:N)
  {
    res[i,]=y[i,]-x[i,]%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  }
  h_full=1.5*sd(c(res))*N^(-1/3)
  B=t(x)-H_hat6%*%t(cbind(1,z))
  g_full<-function(beta)   
  {
    re=matrix(0,N,d)
    beta=as.matrix(beta)
    for (i in 1:N)
    {
      re[i,]=t(c(Tau-kern( (-y[i,]+x[i,]%*%beta+cbind(1,z)[i,]%*%eta_hat6)/h_full ))*matrix(B[,i]) )
    }
    return(re)
  }
  f_full=function(beta)
  {
    f=el.ee(beta,gee=g_full)
    return(f)    
  }
  beta_hat6=optim(beta_hat_6,f_full,method = "BFGS")$par
  beta_full_pqr_1[((r-1)*d+1):(r*d),1]=beta_hat6
  ###########################################################################################csl方法
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
  beta_pqr_1_1[((r-1)*d+1):(r*d),1]=beta_hat4_1
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
  beta_pqr_2_1[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #####################################################################################################csl_all方法
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
  beta_pqr_3_1[((r-1)*d+1):(r*d),1]=beta_hat_9
  ###########################################################################################################
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
  beta_pqr_ave_1[((r-1)*d+1):(r*d),1]=beta_hat_9
###########################################################################################################################
  
  
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
  
  for (m in 1:M){
    beta_hat_8[,m]=beta_all[,m]
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    eta1=matrix(eta_all[,m])
    H1=H_hat[((m-1)*d+1):(m*d),]
    g_sub<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_8[,m]-cbind(1,z1)[i,]%*%eta_all[,m]
      }
      h=1.5*sd(c(res))*n^(-1/5)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t(c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta1)/h ))*(matrix(x1[i,])-H1%*%cbind(1,z1)[i,]))
      }
      return(re)
    }
    f_sub=function(beta)
    {
      f=el.ee(beta,gee=g_sub)
      return(f)    
    }
    beta_hat8[,m]=optim(beta_hat_8[,m],f_sub,method = "BFGS")$par
  }
  beta_ave_pqr_2[((r-1)*d+1):(r*d),1]=apply(beta_hat8,1,mean)
  ####################################################################################sub方法
  beta_sub_pqr_2[((r-1)*d+1):(r*d),1]=beta_hat8[,1]
  ####################################################################################full方法
  fit6 = mvr(x,z)
  H_hat6 = t(rbind(as.vector(fit6$muhat),fit6$Bhat))  
  
  result6 = ini.est(cbind(x,1,z),y,index=c(1:(d+1)),tau=Tau)  
  beta_hat_6=result6$coefficient[1:d]
  eta_hat6=result6$coefficient[(d+1):(d+p)]      
  
  for (j in 1:p){
    if (abs(eta_hat6[j])<0.2){
      eta_hat6[j]=0
    }
  }
  #for (i in 1:N){
  #  y_fullpqr[i,]=y[i]-t(H_hat6%*%cbind(1,z)[i,])%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  #  x_fullpqr[i,]=x[i,]-t(H_hat6%*%cbind(1,z)[i,])
  #}
  #beta_hat6=rq(y_fullpqr~x_fullpqr, tau=Tau)$coefficient[2:(d+1)]
  res=matrix(0,N,1)
  for (i in 1:N)
  {
    res[i,]=y[i,]-x[i,]%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  }
  h_full=1.5*sd(c(res))*N^(-1/5)
  B=t(x)-H_hat6%*%t(cbind(1,z))
  g_full<-function(beta)   
  {
    re=matrix(0,N,d)
    beta=as.matrix(beta)
    for (i in 1:N)
    {
      re[i,]=t(c(Tau-kern( (-y[i,]+x[i,]%*%beta+cbind(1,z)[i,]%*%eta_hat6)/h_full ))*matrix(B[,i]) )
    }
    return(re)
  }
  f_full=function(beta)
  {
    f=el.ee(beta,gee=g_full)
    return(f)    
  }
  beta_hat6=optim(beta_hat_6,f_full,method = "BFGS")$par
  beta_full_pqr_2[((r-1)*d+1):(r*d),1]=beta_hat6
  ###########################################################################################csl方法
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
    h1=1.5*sd(c(res1))*n^(-1/5)
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
  h2=1.5*sd(c(res2))*n^(-1/5)
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
    h=1.5*sd(c(res))*n^(-1/5)
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
  beta_pqr_1_2[((r-1)*d+1):(r*d),1]=beta_hat4_1
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
    h1=1.5*sd(c(res1))*n^(-1/5)
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
  h2=1.5*sd(c(res2))*n^(-1/5)
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
    h=1.5*sd(c(res))*n^(-1/5)
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
  beta_pqr_2_2[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #####################################################################################################csl_all方法
  beta_hat_9=beta_hat
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-1/5)
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
    h3=1.5*sd(c(res3))*n^(-1/5)
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
      h=1.5*sd(c(res))*n^(-1/5)
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
  beta_pqr_3_2[((r-1)*d+1):(r*d),1]=beta_hat_9
  ###########################################################################################################
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
    h1=1.5*sd(c(res1))*n^(-1/5)
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
    h3=1.5*sd(c(res3))*n^(-1/5)
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
      h=1.5*sd(c(res))*n^(-1/5)
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
  beta_pqr_ave_2[((r-1)*d+1):(r*d),1]=beta_hat_9
  ###############################################################################################################
  
  
  
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
  
  for (m in 1:M){
    beta_hat_8[,m]=beta_all[,m]
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    eta1=matrix(eta_all[,m])
    H1=H_hat[((m-1)*d+1):(m*d),]
    g_sub<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_8[,m]-cbind(1,z1)[i,]%*%eta_all[,m]
      }
      h=1.5*sd(c(res))*n^(-2/5)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t(c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta1)/h ))*(matrix(x1[i,])-H1%*%cbind(1,z1)[i,]))
      }
      return(re)
    }
    f_sub=function(beta)
    {
      f=el.ee(beta,gee=g_sub)
      return(f)    
    }
    beta_hat8[,m]=optim(beta_hat_8[,m],f_sub,method = "BFGS")$par
  }
  beta_ave_pqr_3[((r-1)*d+1):(r*d),1]=apply(beta_hat8,1,mean)
  ####################################################################################sub方法
  beta_sub_pqr_3[((r-1)*d+1):(r*d),1]=beta_hat8[,1]
  ####################################################################################full方法
  fit6 = mvr(x,z)
  H_hat6 = t(rbind(as.vector(fit6$muhat),fit6$Bhat))  
  
  result6 = ini.est(cbind(x,1,z),y,index=c(1:(d+1)),tau=Tau)  
  beta_hat_6=result6$coefficient[1:d]
  eta_hat6=result6$coefficient[(d+1):(d+p)]      
  
  for (j in 1:p){
    if (abs(eta_hat6[j])<0.2){
      eta_hat6[j]=0
    }
  }
  #for (i in 1:N){
  #  y_fullpqr[i,]=y[i]-t(H_hat6%*%cbind(1,z)[i,])%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  #  x_fullpqr[i,]=x[i,]-t(H_hat6%*%cbind(1,z)[i,])
  #}
  #beta_hat6=rq(y_fullpqr~x_fullpqr, tau=Tau)$coefficient[2:(d+1)]
  res=matrix(0,N,1)
  for (i in 1:N)
  {
    res[i,]=y[i,]-x[i,]%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  }
  h_full=1.5*sd(c(res))*N^(-2/5)
  B=t(x)-H_hat6%*%t(cbind(1,z))
  g_full<-function(beta)   
  {
    re=matrix(0,N,d)
    beta=as.matrix(beta)
    for (i in 1:N)
    {
      re[i,]=t(c(Tau-kern( (-y[i,]+x[i,]%*%beta+cbind(1,z)[i,]%*%eta_hat6)/h_full ))*matrix(B[,i]) )
    }
    return(re)
  }
  f_full=function(beta)
  {
    f=el.ee(beta,gee=g_full)
    return(f)    
  }
  beta_hat6=optim(beta_hat_6,f_full,method = "BFGS")$par
  beta_full_pqr_3[((r-1)*d+1):(r*d),1]=beta_hat6
  ###########################################################################################csl方法
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
    h1=1.5*sd(c(res1))*n^(-2/5)
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
  h2=1.5*sd(c(res2))*n^(-2/5)
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
    h=1.5*sd(c(res))*n^(-2/5)
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
  beta_pqr_1_3[((r-1)*d+1):(r*d),1]=beta_hat4_1
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
    h1=1.5*sd(c(res1))*n^(-2/5)
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
  h2=1.5*sd(c(res2))*n^(-2/5)
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
    h=1.5*sd(c(res))*n^(-2/5)
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
  beta_pqr_2_3[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #####################################################################################################csl_all方法
  beta_hat_9=beta_hat
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1.5*sd(c(res1))*n^(-2/5)
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
    h3=1.5*sd(c(res3))*n^(-2/5)
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
      h=1.5*sd(c(res))*n^(-2/5)
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
  beta_pqr_3_3[((r-1)*d+1):(r*d),1]=beta_hat_9
  ###########################################################################################################
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
    h1=1.5*sd(c(res1))*n^(-2/5)
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
    h3=1.5*sd(c(res3))*n^(-2/5)
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
      h=1.5*sd(c(res))*n^(-2/5)
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
  beta_pqr_ave_3[((r-1)*d+1):(r*d),1]=beta_hat_9
  #####################################################################################################
  
  
  
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
  
  for (m in 1:M){
    beta_hat_8[,m]=beta_all[,m]
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    eta1=matrix(eta_all[,m])
    H1=H_hat[((m-1)*d+1):(m*d),]
    g_sub<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_8[,m]-cbind(1,z1)[i,]%*%eta_all[,m]
      }
      h=1*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t(c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta1)/h ))*(matrix(x1[i,])-H1%*%cbind(1,z1)[i,]))
      }
      return(re)
    }
    f_sub=function(beta)
    {
      f=el.ee(beta,gee=g_sub)
      return(f)    
    }
    beta_hat8[,m]=optim(beta_hat_8[,m],f_sub,method = "BFGS")$par
  }
  beta_ave_pqr_4[((r-1)*d+1):(r*d),1]=apply(beta_hat8,1,mean)
  ####################################################################################sub方法
  beta_sub_pqr_4[((r-1)*d+1):(r*d),1]=beta_hat8[,1]
  ####################################################################################full方法
  fit6 = mvr(x,z)
  H_hat6 = t(rbind(as.vector(fit6$muhat),fit6$Bhat))  
  
  result6 = ini.est(cbind(x,1,z),y,index=c(1:(d+1)),tau=Tau)  
  beta_hat_6=result6$coefficient[1:d]
  eta_hat6=result6$coefficient[(d+1):(d+p)]      
  
  for (j in 1:p){
    if (abs(eta_hat6[j])<0.2){
      eta_hat6[j]=0
    }
  }
  #for (i in 1:N){
  #  y_fullpqr[i,]=y[i]-t(H_hat6%*%cbind(1,z)[i,])%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  #  x_fullpqr[i,]=x[i,]-t(H_hat6%*%cbind(1,z)[i,])
  #}
  #beta_hat6=rq(y_fullpqr~x_fullpqr, tau=Tau)$coefficient[2:(d+1)]
  res=matrix(0,N,1)
  for (i in 1:N)
  {
    res[i,]=y[i,]-x[i,]%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  }
  h_full=1*sd(c(res))*N^(-1/3)
  B=t(x)-H_hat6%*%t(cbind(1,z))
  g_full<-function(beta)   
  {
    re=matrix(0,N,d)
    beta=as.matrix(beta)
    for (i in 1:N)
    {
      re[i,]=t(c(Tau-kern( (-y[i,]+x[i,]%*%beta+cbind(1,z)[i,]%*%eta_hat6)/h_full ))*matrix(B[,i]) )
    }
    return(re)
  }
  f_full=function(beta)
  {
    f=el.ee(beta,gee=g_full)
    return(f)    
  }
  beta_hat6=optim(beta_hat_6,f_full,method = "BFGS")$par
  beta_full_pqr_4[((r-1)*d+1):(r*d),1]=beta_hat6
  ###########################################################################################csl方法
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
    h1=1*sd(c(res1))*n^(-1/3)
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
  h2=1*sd(c(res2))*n^(-1/3)
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
    h=1*sd(c(res))*n^(-1/3)
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
  beta_pqr_1_4[((r-1)*d+1):(r*d),1]=beta_hat4_1
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
    h1=1*sd(c(res1))*n^(-1/3)
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
  h2=1*sd(c(res2))*n^(-1/3)
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
    h=1*sd(c(res))*n^(-1/3)
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
  beta_pqr_2_4[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #####################################################################################################csl_all方法
  beta_hat_9=beta_hat
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=1*sd(c(res1))*n^(-1/3)
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
    h3=1*sd(c(res3))*n^(-1/3)
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
      h=1*sd(c(res))*n^(-1/3)
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
  beta_pqr_3_4[((r-1)*d+1):(r*d),1]=beta_hat_9
  ###########################################################################################################
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
    h1=1*sd(c(res1))*n^(-1/3)
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
    h3=1*sd(c(res3))*n^(-1/3)
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
      h=1*sd(c(res))*n^(-1/3)
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
  beta_pqr_ave_4[((r-1)*d+1):(r*d),1]=beta_hat_9
  ##############################################################################################################################
  
  
  
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
  
  for (m in 1:M){
    beta_hat_8[,m]=beta_all[,m]
    x1=x[((m-1)*n+1):(m*n),]
    z1=z[((m-1)*n+1):(m*n),]
    y1=matrix(y[((m-1)*n+1):(m*n),])
    eta1=matrix(eta_all[,m])
    H1=H_hat[((m-1)*d+1):(m*d),]
    g_sub<-function(beta)   
    {
      res=matrix(0,n,1)
      for (i in 1:n)
      {
        res[i,]=y1[i,]-x1[i,]%*%beta_hat_8[,m]-cbind(1,z1)[i,]%*%eta_all[,m]
      }
      h=2*sd(c(res))*n^(-1/3)
      re=matrix(0,n,d)
      beta=as.matrix(beta)
      for (i in 1:n)
      {
        re[i,]=t(c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta1)/h ))*(matrix(x1[i,])-H1%*%cbind(1,z1)[i,]))
      }
      return(re)
    }
    f_sub=function(beta)
    {
      f=el.ee(beta,gee=g_sub)
      return(f)    
    }
    beta_hat8[,m]=optim(beta_hat_8[,m],f_sub,method = "BFGS")$par
  }
  beta_ave_pqr_5[((r-1)*d+1):(r*d),1]=apply(beta_hat8,1,mean)
  ####################################################################################sub方法
  beta_sub_pqr_5[((r-1)*d+1):(r*d),1]=beta_hat8[,1]
  ####################################################################################full方法
  fit6 = mvr(x,z)
  H_hat6 = t(rbind(as.vector(fit6$muhat),fit6$Bhat))  
  
  result6 = ini.est(cbind(x,1,z),y,index=c(1:(d+1)),tau=Tau)  
  beta_hat_6=result6$coefficient[1:d]
  eta_hat6=result6$coefficient[(d+1):(d+p)]      
  
  for (j in 1:p){
    if (abs(eta_hat6[j])<0.2){
      eta_hat6[j]=0
    }
  }
  #for (i in 1:N){
  #  y_fullpqr[i,]=y[i]-t(H_hat6%*%cbind(1,z)[i,])%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  #  x_fullpqr[i,]=x[i,]-t(H_hat6%*%cbind(1,z)[i,])
  #}
  #beta_hat6=rq(y_fullpqr~x_fullpqr, tau=Tau)$coefficient[2:(d+1)]
  res=matrix(0,N,1)
  for (i in 1:N)
  {
    res[i,]=y[i,]-x[i,]%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  }
  h_full=2*sd(c(res))*N^(-1/3)
  B=t(x)-H_hat6%*%t(cbind(1,z))
  g_full<-function(beta)   
  {
    re=matrix(0,N,d)
    beta=as.matrix(beta)
    for (i in 1:N)
    {
      re[i,]=t(c(Tau-kern( (-y[i,]+x[i,]%*%beta+cbind(1,z)[i,]%*%eta_hat6)/h_full ))*matrix(B[,i]) )
    }
    return(re)
  }
  f_full=function(beta)
  {
    f=el.ee(beta,gee=g_full)
    return(f)    
  }
  beta_hat6=optim(beta_hat_6,f_full,method = "BFGS")$par
  beta_full_pqr_5[((r-1)*d+1):(r*d),1]=beta_hat6
  ###########################################################################################csl方法
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
    h1=2*sd(c(res1))*n^(-1/3)
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
  h2=2*sd(c(res2))*n^(-1/3)
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
    h=2*sd(c(res))*n^(-1/3)
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
  beta_pqr_1_5[((r-1)*d+1):(r*d),1]=beta_hat4_1
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
    h1=2*sd(c(res1))*n^(-1/3)
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
  h2=2*sd(c(res2))*n^(-1/3)
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
    h=2*sd(c(res))*n^(-1/3)
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
  beta_pqr_2_5[((r-1)*d+1):(r*d),1]=beta_hat4_1
  #####################################################################################################csl_all方法
  beta_hat_9=beta_hat
  
  res1=matrix(0,n,1)
  L_N=0
  for (m in 1:M){
    for (i in 1:n)
    {
      res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat
    }
    h1=2*sd(c(res1))*n^(-1/3)
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
    h3=2*sd(c(res3))*n^(-1/3)
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
      h=2*sd(c(res))*n^(-1/3)
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
  beta_pqr_3_5[((r-1)*d+1):(r*d),1]=beta_hat_9
  ###########################################################################################################
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
    h1=2*sd(c(res1))*n^(-1/3)
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
    h3=2*sd(c(res3))*n^(-1/3)
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
      h=2*sd(c(res))*n^(-1/3)
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
  beta_pqr_ave_5[((r-1)*d+1):(r*d),1]=beta_hat_9
  
  print(r)
}
msepqr1_all_1=rep(0,Run)
madpqr1_all_1=rep(0,Run)
msepqr2_all_1=rep(0,Run)
madpqr2_all_1=rep(0,Run)
msepqr3_all_1=rep(0,Run)
madpqr3_all_1=rep(0,Run)
msepqrave_all_1=rep(0,Run)
madpqrave_all_1=rep(0,Run)
msesubpqr_all_1=rep(0,Run)
madsubpqr_all_1=rep(0,Run)
mseavepqr_all_1=rep(0,Run)
madavepqr_all_1=rep(0,Run)
msefullpqr_all_1=rep(0,Run)
madfullpqr_all_1=rep(0,Run)
for (r in 1:Run){
  msepqr1_all_1[r]=sum( (beta_pqr_1_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr1_all_1[r]=sum( abs(beta_pqr_1_1[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr2_all_1[r]=sum( (beta_pqr_2_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr2_all_1[r]=sum( abs(beta_pqr_2_1[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr3_all_1[r]=sum( (beta_pqr_3_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr3_all_1[r]=sum( abs(beta_pqr_3_1[((r-1)*d+1):(r*d),1]-beta0) )
  msepqrave_all_1[r]=sum( (beta_pqr_ave_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqrave_all_1[r]=sum( abs(beta_pqr_ave_1[((r-1)*d+1):(r*d),1]-beta0) )
  msesubpqr_all_1[r]=sum( (beta_sub_pqr_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madsubpqr_all_1[r]=sum( abs(beta_sub_pqr_1[((r-1)*d+1):(r*d),1]-beta0) )
  mseavepqr_all_1[r]=sum( (beta_ave_pqr_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madavepqr_all_1[r]=sum( abs(beta_ave_pqr_1[((r-1)*d+1):(r*d),1]-beta0) )
  msefullpqr_all_1[r]=sum( (beta_full_pqr_1[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madfullpqr_all_1[r]=sum( abs(beta_full_pqr_1[((r-1)*d+1):(r*d),1]-beta0) )
}
msepqr1_1=mean(msepqr1_all_1)
madpqr1_1=mean(madpqr1_all_1)
msepqr2_1=mean(msepqr2_all_1)
madpqr2_1=mean(madpqr2_all_1)
msepqr3_1=mean(msepqr3_all_1)
madpqr3_1=mean(madpqr3_all_1)
msepqrave_1=mean(msepqrave_all_1)
madpqrave_1=mean(madpqrave_all_1)
msesubpqr_1=mean(msesubpqr_all_1)
madsubpqr_1=mean(madsubpqr_all_1)
mseavepqr_1=mean(mseavepqr_all_1)
madavepqr_1=mean(madavepqr_all_1)
msefullpqr_1=mean(msefullpqr_all_1)
madfullpqr_1=mean(madfullpqr_all_1)
mse_1=c(
  msesubpqr_1,
  mseavepqr_1,
  msefullpqr_1,
  msepqr1_1,
  msepqr2_1,
  msepqr3_1,
  msepqrave_1)
mad_1=c(
  madsubpqr_1,
  madavepqr_1,
  madfullpqr_1,
  madpqr1_1,
  madpqr2_1,
  madpqr3_1,
  madpqrave_1)

msepqr1_all_2=rep(0,Run)
madpqr1_all_2=rep(0,Run)
msepqr2_all_2=rep(0,Run)
madpqr2_all_2=rep(0,Run)
msepqr3_all_2=rep(0,Run)
madpqr3_all_2=rep(0,Run)
msepqrave_all_2=rep(0,Run)
madpqrave_all_2=rep(0,Run)
msesubpqr_all_2=rep(0,Run)
madsubpqr_all_2=rep(0,Run)
mseavepqr_all_2=rep(0,Run)
madavepqr_all_2=rep(0,Run)
msefullpqr_all_2=rep(0,Run)
madfullpqr_all_2=rep(0,Run)
for (r in 1:Run){
  msepqr1_all_2[r]=sum( (beta_pqr_1_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr1_all_2[r]=sum( abs(beta_pqr_1_2[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr2_all_2[r]=sum( (beta_pqr_2_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr2_all_2[r]=sum( abs(beta_pqr_2_2[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr3_all_2[r]=sum( (beta_pqr_3_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr3_all_2[r]=sum( abs(beta_pqr_3_2[((r-1)*d+1):(r*d),1]-beta0) )
  msepqrave_all_2[r]=sum( (beta_pqr_ave_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqrave_all_2[r]=sum( abs(beta_pqr_ave_2[((r-1)*d+1):(r*d),1]-beta0) )
  msesubpqr_all_2[r]=sum( (beta_sub_pqr_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madsubpqr_all_2[r]=sum( abs(beta_sub_pqr_2[((r-1)*d+1):(r*d),1]-beta0) )
  mseavepqr_all_2[r]=sum( (beta_ave_pqr_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madavepqr_all_2[r]=sum( abs(beta_ave_pqr_2[((r-1)*d+1):(r*d),1]-beta0) )
  msefullpqr_all_2[r]=sum( (beta_full_pqr_2[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madfullpqr_all_2[r]=sum( abs(beta_full_pqr_2[((r-1)*d+1):(r*d),1]-beta0) )
}
msepqr1_2=mean(msepqr1_all_2)
madpqr1_2=mean(madpqr1_all_2)
msepqr2_2=mean(msepqr2_all_2)
madpqr2_2=mean(madpqr2_all_2)
msepqr3_2=mean(msepqr3_all_2)
madpqr3_2=mean(madpqr3_all_2)
msepqrave_2=mean(msepqrave_all_2)
madpqrave_2=mean(madpqrave_all_2)
msesubpqr_2=mean(msesubpqr_all_2)
madsubpqr_2=mean(madsubpqr_all_2)
mseavepqr_2=mean(mseavepqr_all_2)
madavepqr_2=mean(madavepqr_all_2)
msefullpqr_2=mean(msefullpqr_all_2)
madfullpqr_2=mean(madfullpqr_all_2)
mse_2=c(
  msesubpqr_2,
  mseavepqr_2,
  msefullpqr_2,
  msepqr1_2,
  msepqr2_2,
  msepqr3_2,
  msepqrave_2)
mad_2=c(
  madsubpqr_2,
  madavepqr_2,
  madfullpqr_2,
  madpqr1_2,
  madpqr2_2,
  madpqr3_2,
  madpqrave_2)

msepqr1_all_3=rep(0,Run)
madpqr1_all_3=rep(0,Run)
msepqr2_all_3=rep(0,Run)
madpqr2_all_3=rep(0,Run)
msepqr3_all_3=rep(0,Run)
madpqr3_all_3=rep(0,Run)
msepqrave_all_3=rep(0,Run)
madpqrave_all_3=rep(0,Run)
msesubpqr_all_3=rep(0,Run)
madsubpqr_all_3=rep(0,Run)
mseavepqr_all_3=rep(0,Run)
madavepqr_all_3=rep(0,Run)
msefullpqr_all_3=rep(0,Run)
madfullpqr_all_3=rep(0,Run)
for (r in 1:Run){
  msepqr1_all_3[r]=sum( (beta_pqr_1_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr1_all_3[r]=sum( abs(beta_pqr_1_3[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr2_all_3[r]=sum( (beta_pqr_2_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr2_all_3[r]=sum( abs(beta_pqr_2_3[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr3_all_3[r]=sum( (beta_pqr_3_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr3_all_3[r]=sum( abs(beta_pqr_3_3[((r-1)*d+1):(r*d),1]-beta0) )
  msepqrave_all_3[r]=sum( (beta_pqr_ave_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqrave_all_3[r]=sum( abs(beta_pqr_ave_3[((r-1)*d+1):(r*d),1]-beta0) )
  msesubpqr_all_3[r]=sum( (beta_sub_pqr_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madsubpqr_all_3[r]=sum( abs(beta_sub_pqr_3[((r-1)*d+1):(r*d),1]-beta0) )
  mseavepqr_all_3[r]=sum( (beta_ave_pqr_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madavepqr_all_3[r]=sum( abs(beta_ave_pqr_3[((r-1)*d+1):(r*d),1]-beta0) )
  msefullpqr_all_3[r]=sum( (beta_full_pqr_3[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madfullpqr_all_3[r]=sum( abs(beta_full_pqr_3[((r-1)*d+1):(r*d),1]-beta0) )
}
msepqr1_3=mean(msepqr1_all_3)
madpqr1_3=mean(madpqr1_all_3)
msepqr2_3=mean(msepqr2_all_3)
madpqr2_3=mean(madpqr2_all_3)
msepqr3_3=mean(msepqr3_all_3)
madpqr3_3=mean(madpqr3_all_3)
msepqrave_3=mean(msepqrave_all_3)
madpqrave_3=mean(madpqrave_all_3)
msesubpqr_3=mean(msesubpqr_all_3)
madsubpqr_3=mean(madsubpqr_all_3)
mseavepqr_3=mean(mseavepqr_all_3)
madavepqr_3=mean(madavepqr_all_3)
msefullpqr_3=mean(msefullpqr_all_3)
madfullpqr_3=mean(madfullpqr_all_3)
mse_3=c(
  msesubpqr_3,
  mseavepqr_3,
  msefullpqr_3,
  msepqr1_3,
  msepqr2_3,
  msepqr3_3,
  msepqrave_3)
mad_3=c(
  madsubpqr_3,
  madavepqr_3,
  madfullpqr_3,
  madpqr1_3,
  madpqr2_3,
  madpqr3_3,
  madpqrave_3)

msepqr1_all_4=rep(0,Run)
madpqr1_all_4=rep(0,Run)
msepqr2_all_4=rep(0,Run)
madpqr2_all_4=rep(0,Run)
msepqr3_all_4=rep(0,Run)
madpqr3_all_4=rep(0,Run)
msepqrave_all_4=rep(0,Run)
madpqrave_all_4=rep(0,Run)
msesubpqr_all_4=rep(0,Run)
madsubpqr_all_4=rep(0,Run)
mseavepqr_all_4=rep(0,Run)
madavepqr_all_4=rep(0,Run)
msefullpqr_all_4=rep(0,Run)
madfullpqr_all_4=rep(0,Run)
for (r in 1:Run){
  msepqr1_all_4[r]=sum( (beta_pqr_1_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr1_all_4[r]=sum( abs(beta_pqr_1_4[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr2_all_4[r]=sum( (beta_pqr_2_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr2_all_4[r]=sum( abs(beta_pqr_2_4[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr3_all_4[r]=sum( (beta_pqr_3_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr3_all_4[r]=sum( abs(beta_pqr_3_4[((r-1)*d+1):(r*d),1]-beta0) )
  msepqrave_all_4[r]=sum( (beta_pqr_ave_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqrave_all_4[r]=sum( abs(beta_pqr_ave_4[((r-1)*d+1):(r*d),1]-beta0) )
  msesubpqr_all_4[r]=sum( (beta_sub_pqr_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madsubpqr_all_4[r]=sum( abs(beta_sub_pqr_4[((r-1)*d+1):(r*d),1]-beta0) )
  mseavepqr_all_4[r]=sum( (beta_ave_pqr_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madavepqr_all_4[r]=sum( abs(beta_ave_pqr_4[((r-1)*d+1):(r*d),1]-beta0) )
  msefullpqr_all_4[r]=sum( (beta_full_pqr_4[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madfullpqr_all_4[r]=sum( abs(beta_full_pqr_4[((r-1)*d+1):(r*d),1]-beta0) )
}
msepqr1_4=mean(msepqr1_all_4)
madpqr1_4=mean(madpqr1_all_4)
msepqr2_4=mean(msepqr2_all_4)
madpqr2_4=mean(madpqr2_all_4)
msepqr3_4=mean(msepqr3_all_4)
madpqr3_4=mean(madpqr3_all_4)
msepqrave_4=mean(msepqrave_all_4)
madpqrave_4=mean(madpqrave_all_4)
msesubpqr_4=mean(msesubpqr_all_4)
madsubpqr_4=mean(madsubpqr_all_4)
mseavepqr_4=mean(mseavepqr_all_4)
madavepqr_4=mean(madavepqr_all_4)
msefullpqr_4=mean(msefullpqr_all_4)
madfullpqr_4=mean(madfullpqr_all_4)
mse_4=c(
  msesubpqr_4,
  mseavepqr_4,
  msefullpqr_4,
  msepqr1_4,
  msepqr2_4,
  msepqr3_4,
  msepqrave_4)
mad_4=c(
  madsubpqr_4,
  madavepqr_4,
  madfullpqr_4,
  madpqr1_4,
  madpqr2_4,
  madpqr3_4,
  madpqrave_4)

msepqr1_all_5=rep(0,Run)
madpqr1_all_5=rep(0,Run)
msepqr2_all_5=rep(0,Run)
madpqr2_all_5=rep(0,Run)
msepqr3_all_5=rep(0,Run)
madpqr3_all_5=rep(0,Run)
msepqrave_all_5=rep(0,Run)
madpqrave_all_5=rep(0,Run)
msesubpqr_all_5=rep(0,Run)
madsubpqr_all_5=rep(0,Run)
mseavepqr_all_5=rep(0,Run)
madavepqr_all_5=rep(0,Run)
msefullpqr_all_5=rep(0,Run)
madfullpqr_all_5=rep(0,Run)
for (r in 1:Run){
  msepqr1_all_5[r]=sum( (beta_pqr_1_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr1_all_5[r]=sum( abs(beta_pqr_1_5[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr2_all_5[r]=sum( (beta_pqr_2_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr2_all_5[r]=sum( abs(beta_pqr_2_5[((r-1)*d+1):(r*d),1]-beta0) )
  msepqr3_all_5[r]=sum( (beta_pqr_3_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqr3_all_5[r]=sum( abs(beta_pqr_3_5[((r-1)*d+1):(r*d),1]-beta0) )
  msepqrave_all_5[r]=sum( (beta_pqr_ave_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madpqrave_all_5[r]=sum( abs(beta_pqr_ave_5[((r-1)*d+1):(r*d),1]-beta0) )
  msesubpqr_all_5[r]=sum( (beta_sub_pqr_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madsubpqr_all_5[r]=sum( abs(beta_sub_pqr_5[((r-1)*d+1):(r*d),1]-beta0) )
  mseavepqr_all_5[r]=sum( (beta_ave_pqr_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madavepqr_all_5[r]=sum( abs(beta_ave_pqr_5[((r-1)*d+1):(r*d),1]-beta0) )
  msefullpqr_all_5[r]=sum( (beta_full_pqr_5[((r-1)*d+1):(r*d),1]-beta0)^2 )
  madfullpqr_all_5[r]=sum( abs(beta_full_pqr_5[((r-1)*d+1):(r*d),1]-beta0) )
}
msepqr1_5=mean(msepqr1_all_5)
madpqr1_5=mean(madpqr1_all_5)
msepqr2_5=mean(msepqr2_all_5)
madpqr2_5=mean(madpqr2_all_5)
msepqr3_5=mean(msepqr3_all_5)
madpqr3_5=mean(madpqr3_all_5)
msepqrave_5=mean(msepqrave_all_5)
madpqrave_5=mean(madpqrave_all_5)
msesubpqr_5=mean(msesubpqr_all_5)
madsubpqr_5=mean(madsubpqr_all_5)
mseavepqr_5=mean(mseavepqr_all_5)
madavepqr_5=mean(madavepqr_all_5)
msefullpqr_5=mean(msefullpqr_all_5)
madfullpqr_5=mean(madfullpqr_all_5)
mse_5=c(
  msesubpqr_5,
  mseavepqr_5,
  msefullpqr_5,
  msepqr1_5,
  msepqr2_5,
  msepqr3_5,
  msepqrave_5)
mad_5=c(
  madsubpqr_5,
  madavepqr_5,
  madfullpqr_5,
  madpqr1_5,
  madpqr2_5,
  madpqr3_5,
  madpqrave_5)

results=rbind(mse_1,mad_1,mse_2,mad_2,mse_3,mad_3,mse_4,mad_4,mse_5,mad_5)
results=matrix(results,10)
