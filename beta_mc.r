# Beta - Simulação Monte Carlo

beta_mc <- function(x, alpha=5, beta=3, n=50, R=100){
  
inicio<-Sys.time()

set.seed(369)

betafit<-function(x)
{
  n<-length(x)
  
  #MM
  x_bar<-mean(x)
  sigma2_hat<-(1/(n-1))*sum((x-x_bar)^2)
  
  alpha_hat<-x_bar*(((x_bar*(1-x_bar))/sigma2_hat)-1)
  beta_hat<-(1-x_bar)*(((x_bar*(1-x_bar))/sigma2_hat)-1)
  
  par_ini<-c(alpha_hat,beta_hat)
  
  #MV
  loglik<-function(par)
  {
    alpha<-par[1]
    beta<-par[2]
    
    suppressWarnings(sum(log(dbeta(x,shape1=alpha,shape2=beta))))
    # suppressWarnings(-n*lbeta(alpha,beta) + (alpha-1)*sum(log(x)) + (beta-1)*sum(log(1-x)))
  }
  
  escore<-function(par)
  {
    alpha<-par[1]
    beta<-par[2]
    
    Ualpha<-n*digamma(alpha+beta)-n*digamma(alpha)+sum(log(x))
    Ubeta<-n*digamma(alpha+beta)-n*digamma(beta)+sum(log(1-x))
    
    c(Ualpha,Ubeta) 
  }
  
  max<-optim(par_ini, loglik, escore, method="BFGS", control=list(fnscale=-1))
  
  z<-c()
  z$mm<-c(alpha_hat,beta_hat)
  z$mv<-max$par
  z$conv<-max$convergence
  
  return(z)
}  

result<-c(); i=falhas=0
while(i < R)
{ 
  x<-rbeta(n,shape1=alpha,shape2=beta)
  fit<-betafit(x)
  
  if(fit$conv == 0){
    result<-rbind(result,c(fit$mm,fit$mv))
    i=i+1
  } 
  else falhas=falhas+1
  
  # if( (100*i/R)%%25==0) print(c((100*i/R),"% --------------- R ---------------"), quote = F)
}
result
medias<-colMeans(result)
real_par<-rep(c(alpha,beta),2)

#vies
vies<-medias-real_par

#VR
vr<- (vies/real_par)*100

#erro padrao
sd <- apply(result,2,sd) 

#EQM
eqm <- apply(result,2,var)+vies^2 

# final
all<-cbind(real_par,medias,vies,vr,sd,eqm)

cat("\n")
print(paste("Contador Falhas:",falhas),q=F)
print(paste("n =",n),q=F)
rownames(all)<-c("alpha MM","beta  MM","alpha MV","beta  MV")
print(round(all,3))

fim<-Sys.time()
print(fim-inicio)

}