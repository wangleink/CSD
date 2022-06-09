library(Matrix); 
library(MASS); 
library(pqr)
library(flare)
data(eyedata)
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
###################
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
kern=function(v){     # smoothing function
  return(pnorm(v,0,1)) 
}
Kern<- function(x,h_1) {
  #re=1/h*exp(x/h)/(1+exp(x/h))^2;
  re=exp(-(x/h_1)^2/2)/(sqrt(2*pi)*h_1)
  #re=(-15/8*(x/h)^2+9/8)*(abs(x/h)<1)
  #re=3/4*(1-(x/h)^2)*(abs(x/h)<1)
  return(re)
}

data=apply(cbind(x,y),2,function(x){(x-mean(x))/(sd(x))})
N=120
n <- 24
M <- 5
N <- n*M
d <- 2  #the column of x
#s <- 2   #z
p <- 199 #the column of z
P=d+p-1
Tau <- 0.5
Run=1
r=1
beta_all=matrix(0,d,M)
eta_all=matrix(0,p,M)
H_hat=matrix(0,(d*M),p)

theta1=matrix(0,(d+p),M)
beta_hat1=matrix(0,d,M)
beta_hat8=matrix(0,d,M)
beta_hat_8=matrix(0,d,M)
beta_hat9=matrix(0,d,M)
beta_hat_9=matrix(0,d,M)
result=matrix(0,(p+d),M)

beta_pqr_1_all=matrix(0,d,Run)
beta_pqr_2_all=matrix(0,d,Run)
beta_pqr_3_all=matrix(0,d,Run)

beta_pqr_ave_all=matrix(0,d,Run)

beta_sub_all=matrix(0,d,Run)
beta_sub_pqr_all=matrix(0,d,Run)

beta_ave_all=matrix(0,d,Run)
beta_ave_pqr_all=matrix(0,d,Run)

beta_full_all=matrix(0,d,Run)
beta_full_pqr_all=matrix(0,d,Run)

y_fullpqr=matrix(0,N,1)
x_fullpqr=matrix(0,N,d)
x1=matrix(0,n,d)
z1=matrix(0,n,(p-1))
y1=matrix(0,n,1)
eta1=matrix(0,p,1)
H1=matrix(0,d,p)

X=data[(1:N),1:200]
y=matrix(data[(1:N),201])
x=X[(1:N),c(140,153)]
z=X[(1:N),-c(140,153)]
#y <- x%*%beta0 + cbind(1,z)%*%eta0 +rt(N,2)#rexp(N,1)rnorm(N,0,1)0.9*rnorm(N,-1,4)+0.1*rnorm(N,9,4) 
######################################################################################average方法
t1_ave=proc.time()
for (m in 1:M){
  result1 = rq.fit.lasso(cbind(1,x[((m-1)*n+1):(m*n),],z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),],tau=Tau,lambda = 1)  
  beta_hat1[,m]=result1$coefficient[2:(d+1)]
}
beta_ave_all[,r]=apply(beta_hat1,1,mean)
t2_ave=proc.time()


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

t1_avepqr=proc.time()
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
beta_ave_pqr_all[,r]=apply(t(na.omit(t(beta_hat8))),1,mean)
t2_avepqr=proc.time()
####################################################################################sub方法
t1_sub=proc.time()
m=1
result1 = rq.fit.lasso(cbind(1,x[((m-1)*n+1):(m*n),],z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),],tau=Tau,lambda = 1)  
beta_hat1[,m]=result1$coefficient[2:(d+1)]

beta_sub_all[,r]=beta_hat1[,1]
t2_sub=proc.time()



beta_hat_8[,m]=beta_all[,m]
x1=x[((m-1)*n+1):(m*n),]
z1=z[((m-1)*n+1):(m*n),]
y1=matrix(y[((m-1)*n+1):(m*n),])
eta1=matrix(eta_all[,m])
H1=H_hat[((m-1)*d+1):(m*d),]
t1_subpqr=proc.time()
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

beta_sub_pqr_all[,r]=beta_hat8[,1]
t2_subpqr=proc.time()
####################################################################################full方法
t1_full=proc.time()
beta_full_all[,r]=rq.fit.lasso(cbind(1,x,z),y,tau=Tau,lambda = 1)$coefficient[2:(d+1)]
t2_full=proc.time()
theta_full=rq.fit.lasso(cbind(1,x,z),y,tau=Tau,lambda = 1)$coefficient


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

t1_fullpqr=proc.time()
g_full<-function(beta)   
{
  res=matrix(0,N,1)
  for (i in 1:N)
  {
    res[i,]=y[i,]-x[i,]%*%beta_hat_6-cbind(1,z)[i,]%*%eta_hat6
  }
  h=1.5*sd(c(res))*N^(-1/3)
  re=matrix(0,N,d)
  beta=as.matrix(beta)
  for (i in 1:N)
  {
    re[i,]=t(c(Tau-kern( (-y[i,]+x[i,]%*%beta+cbind(1,z)[i,]%*%eta_hat6)/h ))*(matrix(x[i,])-H_hat6%*%cbind(1,z)[i,]) )
    # re[i,]=t( c( Tau-ifelse((y[i,]-x[i,]%*%beta-cbind(1,z)[i,]%*%eta_hat6)<0,1,0) )*(matrix(x[i,])-H_hat6%*%cbind(1,z)[i,]) )
    
  }
  return(re)
}
f_full=function(beta)
{
  f=el.ee(beta,gee=g_full)
  return(f)    
}
beta_hat6=optim(beta_hat_6,f_full,method = "BFGS")$par#rq(y_fullpqr~x_fullpqr, tau=Tau)$coefficient[2:(d+1)]
beta_full_pqr_all[,r]=beta_hat6
t2_fullpqr=proc.time()
######################################################################################################################################################
H_csl=matrix(0,d,p)
for (m in 1:M){
  H_csl=H_csl+H_hat[((m-1)*d+1):(m*d),]/M
} 
t1_cslpqr=proc.time()
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
beta_pqr_1_all[,r]=beta_hat4_1
t2_cslpqr=proc.time()
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
t3_cslpqr=proc.time()
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
beta_pqr_2_all[,r]=beta_hat4_1
t4_cslpqr=proc.time()
#####################################################################################################
t1_gelpqr=proc.time()
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
beta_pqr_3_all[,r]=beta_hat_9
t2_gelpqr=proc.time()
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
eta_hat_ave=apply(eta_all,1,mean) 
t3_gelpqr=proc.time()
res1=matrix(0,n,1)
L_N=0
for (m in 1:M){
  for (i in 1:n)
  {
    res1[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat_ave
  }
  h1=1.5*sd(c(res1))*n^(-1/3)
  for (i in 1:n){
    L_N=L_N+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat_ave)/h1 ) )*(matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
  }
}
L_N=L_N/N

for (m in 1:M){
  res3=matrix(0,n,1)
  L_m=0
  for (i in 1:n)
  {
    res3[i,]=y[((m-1)*n+i),]-x[((m-1)*n+i),]%*%beta_hat_9-cbind(1,z)[((m-1)*n+i),]%*%eta_hat_ave
  }
  h3=1.5*sd(c(res3))*n^(-1/3)
  for (i in 1:n){
    L_m=L_m+c( Tau-kern( (-y[((m-1)*n+i),]+x[((m-1)*n+i),]%*%beta_hat_9+cbind(1,z)[((m-1)*n+i),]%*%eta_hat_ave)/h3 ) )*( matrix(x[((m-1)*n+i),])-H_csl%*%cbind(1,z)[((m-1)*n+i),])
    
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
      res[i,]=y1[i,]-x1[i,]%*%beta_hat_9-cbind(1,z1)[i,]%*%eta_hat_ave
    }
    h=1.5*sd(c(res))*n^(-1/3)
    re=matrix(0,n,d)
    beta=as.matrix(beta)
    for (i in 1:n)
    {
      re[i,]=t( c(Tau-kern( (-y1[i,]+x1[i,]%*%beta+cbind(1,z1)[i,]%*%eta_hat_ave)/h ))*(matrix(x1[i,])- H_csl%*%cbind(1,z1)[i,]) +L_N-L_m )
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
beta_pqr_ave_all[,r]=beta_hat_9
t4_gelpqr=proc.time()
beta_sub_all 
beta_sub_pqr_all
beta_ave_all
beta_ave_pqr_all
beta_full_all
beta_full_pqr_all
beta_pqr_1_all
beta_pqr_2_all
beta_pqr_3_all
beta_pqr_ave_all
(t2_sub-t1_sub)[3]
(t2_ave-t1_ave)[3]
(t2_full-t1_full)[3]
(t2_subpqr-t1_subpqr)[3]
(t2_avepqr-t1_avepqr)[3]
(t2_fullpqr-t1_fullpqr)[3]
(t2_cslpqr-t1_cslpqr)[3]
(t4_cslpqr-t3_cslpqr)[3]
(t2_gelpqr-t1_gelpqr)[3]
(t4_gelpqr-t3_gelpqr)[3]
beta_sub_pqr=beta_sub_pqr_all
beta_ave_pqr=beta_ave_pqr_all
beta_full_pqr=beta_full_pqr_all
beta_pqr_1=beta_pqr_1_all
beta_pqr_2=beta_pqr_2_all
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
eta_all_ave=matrix(0,p,M)

up_pqr_1=matrix(0,1,Run)
down_pqr_1=matrix(0,1,Run)
up_pqr_2=matrix(0,1,Run)
down_pqr_2=matrix(0,1,Run)
up_pqr_ave=matrix(0,1,Run)
down_pqr_ave=matrix(0,1,Run)
up_full_pqr=matrix(0,1,Run)
down_full_pqr=matrix(0,1,Run)
up_sub_pqr=matrix(0,1,Run)
down_sub_pqr=matrix(0,1,Run)
up_ave_pqr=matrix(0,1,Run)
down_ave_pqr=matrix(0,1,Run)
Run=1
r=1
v=matrix(c(0,1),d,1)
error=matrix(0,n,1)
for (r in 1:Run){
  for (m in 1:M){
    result8 = ini.est(cbind(1,z[((m-1)*n+1):(m*n),]),y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%c(beta_ave_pqr[((r-1)*d+1):(r*d),1]),index=c(1),tau=Tau)  
    eta_all_ave[,m]=result8$coefficient                       #每一个机器算一次beta eta
    for (j in 1:p){
      if (abs(eta_all_ave[j,m])<0.2){
        eta_all_ave[j,m]=0
      }
    }
  }
  eta_ave_pqr=apply(eta_all,1,mean)
  var=0
  for (m in 1:M){
    D=matrix(0,d,d)
    for (i in 1:n){
      D=D+( (matrix(x[((m-1)*n+i),])-H_hat[((m-1)*d+1):(m*d),]%*%matrix(cbind(1,z)[((m-1)*n+i),])) %*% t(matrix(x[((m-1)*n+i),])-H_hat[((m-1)*d+1):(m*d),]%*%matrix(cbind(1,z)[((m-1)*n+i),])) ) 
    }
    D=Tau*(1-Tau)*(D/n)
    
    Q=matrix(0,d,d)
    error=matrix(y[(((m-1)*n+1):(m*n)),])-x[(((m-1)*n+1):(m*n)),]%*%c(beta_ave_pqr[((r-1)*d+1):(r*d),1])-cbind(1,z)[(((m-1)*n+1):(m*n)),]%*%eta_all_ave[,m]
    h_err=0.9*n^(-1/5)*min(sd(error),IQR(error)/1.34)
    f_err=mean(Kern(error,h_err))
    for (i in 1:n){
      Q=Q+ f_err*( (matrix(x[((m-1)*n+i),])- H_hat[((m-1)*d+1):(m*d),]%*%matrix(cbind(1,z)[((m-1)*n+i),]))%*%t(matrix(x[((m-1)*n+i),])) )
    }
    Q=Q/n
    var=var+t(v)%*%solve(Q)%*%D%*%solve(t(Q))%*%v
  }
  var=var/M
  
  up_ave_pqr[,r]=(t(v)%*%beta_ave_pqr[((r-1)*d+1):(r*d),1])+sqrt(var/N)*qnorm(0.975)
  down_ave_pqr[,r]=(t(v)%*%beta_ave_pqr[((r-1)*d+1):(r*d),1])-sqrt(var/N)*qnorm(0.975)
  var_ave=sqrt(var/N)
  ####################################################################################sub方法
  
  
  result8 = ini.est(cbind(1,z[(1:n),]),y[(1:n),]-x[(1:n),]%*%c(beta_sub_pqr[((r-1)*d+1):(r*d),1]),index=c(1),tau=Tau)  
  eta_sub_pqr=result8$coefficient                       #每一个机器算一次beta eta
  for (j in 1:p){
    if (abs(eta_sub_pqr[j])<0.2){
      eta_sub_pqr[j]=0
    }
  }
  
  D=matrix(0,d,d)
  for (i in 1:n){
    D=D+( (matrix(x[i,])- H_hat[1:d,]%*%matrix(cbind(1,z)[i,])) %*% t(matrix(x[i,])- H_hat[1:d,]%*%matrix(cbind(1,z)[i,])) ) 
  }
  D=Tau*(1-Tau)*(D/n)
  
  Q=matrix(0,d,d)
  error=matrix(y[1:n,])-x[1:n,]%*%c(beta_sub_pqr[((r-1)*d+1):(r*d),1])-cbind(1,z)[1:n,]%*%eta_sub_pqr
  h_err=0.9*n^(-1/5)*min(sd(error),IQR(error)/1.34)
  f_err=mean(Kern(error,h_err))
  for (i in 1:n){
    Q=Q+ f_err*( (matrix(x[i,])- H_hat[1:d,]%*%matrix(cbind(1,z)[i,]))%*%t(matrix(x[i,])) )
  }
  Q=Q/n
  var=t(v)%*%solve(Q)%*%D%*%solve(t(Q))%*%v
  up_sub_pqr[,r]=(t(v)%*%beta_sub_pqr[((r-1)*d+1):(r*d),1])+sqrt(var/n)*qnorm(0.975)
  down_sub_pqr[,r]=(t(v)%*%beta_sub_pqr[((r-1)*d+1):(r*d),1])-sqrt(var/n)*qnorm(0.975)
  var_sub=sqrt(var/n)
  ####################################################################################full方法
  result8 = ini.est(cbind(1,z[(1:N),]),y[(1:N),]-x[(1:N),]%*%c(beta_full_pqr[((r-1)*d+1):(r*d),1]),index=c(1),tau=Tau)  
  eta_full_pqr=result8$coefficient                       #每一个机器算一次beta eta
  for (j in 1:p){
    if (abs(eta_full_pqr[j])<0.2){
      eta_full_pqr[j]=0
    }
  }
  D=matrix(0,d,d)
  for (i in 1:N){
    D=D+( (matrix(x[i,])- H_hat6%*%matrix(cbind(1,z)[i,])) %*% t(matrix(x[i,])- H_hat6%*%matrix(cbind(1,z)[i,])) ) 
  }
  D=Tau*(1-Tau)*(D/N)
  
  Q=matrix(0,d,d)
  error=matrix(y-x%*%c(beta_full_pqr[((r-1)*d+1):(r*d),1])-cbind(1,z)%*%eta_full_pqr)[1:N,1]
  h_err=0.9*N^(-1/5)*min(sd(error),IQR(error)/1.34)
  f_err=mean(Kern(error,h_err))
  for (i in 1:N){
    Q=Q+ f_err*( (matrix(x[i,])- H_hat6%*%matrix(cbind(1,z)[i,]))%*%t(matrix(x[i,])) )
  }
  Q=Q/N
  var=t(v)%*%solve(Q)%*%D%*%solve(t(Q))%*%v
  up_full_pqr[,r]=(t(v)%*%beta_full_pqr[((r-1)*d+1):(r*d),1])+sqrt(var/N)*qnorm(0.975)
  down_full_pqr[,r]=(t(v)%*%beta_full_pqr[((r-1)*d+1):(r*d),1])-sqrt(var/N)*qnorm(0.975)
  var_full=sqrt(var/N)
  ######################################################################################################################################################
  var=0
  for (m in 1:M){
    D=matrix(0,d,d)
    for (i in 1:n){
      D=D+( (matrix(x[((m-1)*n+i),])-H_csl%*%matrix(cbind(1,z)[((m-1)*n+i),])) %*% t(matrix(x[((m-1)*n+i),])-H_csl%*%matrix(cbind(1,z)[((m-1)*n+i),])) ) 
    }
    D=Tau*(1-Tau)*(D/n)
    
    error=y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%c(beta_pqr_1[((r-1)*d+1):(r*d),1])-cbind(1,z)[((m-1)*n+1):(m*n),]%*%eta_all[,m]
    h_err=0.9*n^(-1/5)*min(sd(error),IQR(error)/1.34)
    f_err=mean(Kern(error,h_err))
    
    Q=matrix(0,d,d)
    for (i in 1:n){
      Q=Q+ ( (matrix(x[((m-1)*n+i),])-H_csl%*%matrix(cbind(1,z)[((m-1)*n+i),]))%*%t(matrix(x[((m-1)*n+i),])) )
    }
    Q=f_err*Q/n
    var=var+t(v)%*%solve(Q)%*%D%*%solve(t(Q))%*%v
  }
  var=var/M
  up_pqr_1[,r]=(t(v)%*%beta_pqr_1[((r-1)*d+1):(r*d),1])+sqrt(var/N)*qnorm(0.975)
  down_pqr_1[,r]=(t(v)%*%beta_pqr_1[((r-1)*d+1):(r*d),1])-sqrt(var/N)*qnorm(0.975)
  var_pqr1=sqrt(var/N)
  #################################################################################################
  var=0
  for (m in 1:M){
    D=matrix(0,d,d)
    for (i in 1:n){
      D=D+( (matrix(x[((m-1)*n+i),])-H_csl%*%matrix(cbind(1,z)[((m-1)*n+i),])) %*% t(matrix(x[((m-1)*n+i),])-H_csl%*%matrix(cbind(1,z)[((m-1)*n+i),])) ) 
    }
    D=Tau*(1-Tau)*(D/n)
    
    error=y[((m-1)*n+1):(m*n),]-x[((m-1)*n+1):(m*n),]%*%c(beta_pqr_2[((r-1)*d+1):(r*d),1])-cbind(1,z)[((m-1)*n+1):(m*n),]%*%eta_all[,m]
    h_err=0.9*n^(-1/5)*min(sd(error),IQR(error)/1.34)
    f_err=mean(Kern(error,h_err))
    
    Q=matrix(0,d,d)
    for (i in 1:n){
      Q=Q+ ( (matrix(x[((m-1)*n+i),])-H_csl%*%matrix(cbind(1,z)[((m-1)*n+i),]))%*%t(matrix(x[((m-1)*n+i),])) )
    }
    Q=f_err*Q/n
    var=var+t(v)%*%solve(Q)%*%D%*%solve(t(Q))%*%v
  }
  var=var/M
  up_pqr_2[,r]=(t(v)%*%beta_pqr_2[((r-1)*d+1):(r*d),1])+sqrt(var/N)*qnorm(0.975)
  down_pqr_2[,r]=(t(v)%*%beta_pqr_2[((r-1)*d+1):(r*d),1])-sqrt(var/N)*qnorm(0.975)
  var_pqr2=sqrt(var/N)
  print(r)
}
beta_sub_all 
beta_sub_pqr_all
beta_ave_all
beta_ave_pqr_all
beta_full_all
beta_full_pqr_all
beta_pqr_1_all
beta_pqr_2_all
beta_pqr_3_all
beta_pqr_ave_all
var_sub
var_ave
var_full
var_pqr1
var_pqr2

down_sub_pqr
up_sub_pqr
down_ave_pqr
up_ave_pqr
down_full_pqr
up_full_pqr
down_pqr_1
up_pqr_1
down_pqr_2
up_pqr_2  
 
c(mean(up_pqr_1-down_pqr_1),
  mean(up_pqr_2-down_pqr_2),
  mean(up_full_pqr-down_full_pqr),
  mean(up_ave_pqr-down_ave_pqr),
  mean(up_sub_pqr-down_sub_pqr))
