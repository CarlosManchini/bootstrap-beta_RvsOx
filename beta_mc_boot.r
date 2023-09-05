inicio<-Sys.time()

alpha<-5
beta<-3

n<-30

R<-1000
B<-500

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
    # (-n*lbeta(alpha,beta) + (alpha-1)*sum(log(x)) + (beta-1)*sum(log(1-x)))
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

result=resultBmm=resultBmv=c(); i=falhas=j=0
while(i < R)
{ 
  x<-rbeta(n,shape1=alpha,shape2=beta)
  
  fit<-betafit(x)
  
  if(fit$conv == 0){
    result<-rbind(result,c(fit$mm,fit$mv))
    i=i+1
  } 
  else falhas=falhas+1
  
  for(j in 1:B){
    x_boot_mm <- rbeta(n,shape1=fit$mm[1],shape2=fit$mm[2])                                              # reamostra bootstrap paramétrica
    fit_boot_mm<-betafit(x_boot_mm)                                                                    # faz com que a reamostra passe pelas mesmas funções
    resultBmm<-rbind(resultBmm,c(fit_boot_mm$mm,fit_boot_mm$mv))
    
    x_boot_mv <- rbeta(n,shape1=fit$mv[1],shape2=fit$mv[2])
    fit_boot_mv<-betafit(x_boot_mv)
    resultBmv<-rbind(resultBmv,c(fit_boot_mv$mm,fit_boot_mv$mv))
    
  }
  if( (100*i/R)%%20==0) print(c((100*i/R),"% --------------- R ---------------"), quote = F)
  
}

#MM corrigido por MM
# boot_emmA_mm<- 2*mean(result[,1]) - mean(resultBmm[,1]) #ou theta_hat - (theta_boot - theta_hat)
# boot_emmB_mm<- 2*mean(result[,2]) - mean(resultBmm[,2]) 

#MM corrigido por MV
boot_emmA_mv<- 2*mean(result[,1]) - mean(resultBmv[,1]) #take care bro
boot_emmB_mv<- 2*mean(result[,2]) - mean(resultBmv[,2])

#MV corrigido por MM
# boot_emvA_mm<- 2*mean(result[,3]) - mean(resultBmv[,1])
# boot_emvB_mm<- 2*mean(result[,4]) - mean(resultBmv[,2])

#MV corrigido por MV
boot_emvA_mv<- 2*mean(result[,3]) - mean(resultBmv[,3])
boot_emvB_mv<- 2*mean(result[,4]) - mean(resultBmv[,4])

# medias_boot_mm <-c(boot_emmA_mm,boot_emmB_mm, boot_emvA_mm,boot_emvB_mm)
medias_boot_mv <-c(boot_emmA_mv,boot_emmB_mv, boot_emvA_mv,boot_emvB_mv)

# result
medias<-c(colMeans(result),medias_boot_mv)#,medias_boot_mm)
real_par<-rep(c(alpha,beta),4) #6

#vies
vies<-medias-real_par

#VR
vr<- (vies/real_par)*100

#erro padrao
sd <- apply(result,2,sd) 

#EQM
eqm <- apply(result,2,var)+vies^2 

# resultados
all<-cbind(real_par,medias,vies,vr,sd,eqm)

cat("\n")
print(paste("Contador Falhas:",falhas),q=F)
print(paste("n =",n),q=F)
rownames(all)<-c("alpha MM","beta  MM","alpha MV","beta  MV",
                 "alpha MM*_mv","beta  MM*_mv","alpha MV*_mv","beta  MV*_mv")
                  # "alpha MM*_mm","beta  MM*_mm","alpha MV*_mm","beta  MV*_mm"
print(round(all,3))

fim<-Sys.time()
print(fim-inicio)

