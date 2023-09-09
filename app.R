library(shiny)
library(shinythemes)


# Define UI for random distribution app ----
ui <- fluidPage(
  # shinythemes::themeSelector(),
  theme = shinytheme("lumen"),
  
  withMathJax(),
  
  # App title ----
  titlePanel(h3("Beta Distribution")),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      
      
      # Input: Select the random distribution type ----
      div("Probability density function:",style="text-indent:20px;font-size:125%;"),
      
      withMathJax(helpText("$$f(x;\\alpha,\\beta)=\\frac{\\Gamma(\\alpha+\\beta)}{\\Gamma(\\alpha)\\Gamma(\\beta)}  x^{\\alpha-1}  (1-x)^{\\beta-1}$$",style="text-indent:20px;font-size:100%;")),
      withMathJax(helpText("$$x \\in (0,1) \\qquad \\alpha>0, \\quad \\beta>0$$",style="text-indent:20px;font-size:100%;")),
      
      
      # div(textOutput("fdp"),style="text-indent:20px;font-size:125%;"),
      
      # br() element to introduce extra vertical spacing ----
      br(),
      
      
      
      sliderInput("sliderAlpha", 
                  "\\(\\alpha\\) parameter value:", 
                  min = 1, 
                  max = 15, 
                  value = 2, 
                  animate = animationOptions(interval = 300)),
      
      sliderInput("sliderBeta", 
                  "\\(\\beta\\) parameter value:", 
                  min = 1, 
                  max = 15, 
                  value = 5, 
                  animate = animationOptions(interval = 200))
      
      
    ),
    
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  tabPanel("Density", plotOutput("densityplot")),
                  tabPanel("Monte Carlo Simulation", 
                           br(),
                    sidebarPanel(
                      numericInput("n_mc","Sample size",value=30, min=3, max=4500),
                      numericInput("R_mc","Monte Carlo replications",value = 500, min = 2, max = 50000),
                      tags$hr(),
                      checkboxInput("boot_mc","Bootstrap bias correction"),
                      helpText(h6("Note: It can take a while")),
                      conditionalPanel(condition = "input.boot_mc",
                                       numericInput("B_mc_boot","Bootstrap replications",value = 10, min = 2, max = 1000)
                                       )
                      # numericInput("B_mc_boot","Bootstrap replications",value = 100, min = 2, max = 1000)
                    ), 
                    titlePanel(h4("Results of Monte Carlo simulation for the point estimation of the beta distribution parameters",align="center")),
                    verbatimTextOutput("beta_mc"), verbatimTextOutput("beta_mc_boot")
                  )
      )
      
    )
  )
)

server <- function(input, output) {
  
  
  output$densityplot <- renderPlot({
    alpha <- input$sliderAlpha
    beta<-input$sliderBeta
    
    xrange <- seq(from=0.001 , to=1, by=0.009)
    xprob <- dbeta(xrange , alpha , beta)
    mainlabel <- expression (paste ( "Probability Density Beta ( " ,alpha,", ",beta," )", sep = " " ) )
    plot(xrange , xprob , type = "l" , main = mainlabel , cex=.9, xlab="x" , ylab = "density",lty=1,lwd=3,col="darkblue")
    legend("topright",cex=1.7, bty="n",
           legend = c(sapply(paste("",alpha,"   "), function(x) as.expression(substitute(alpha==B,list(B=as.name(x))))),
                      sapply(paste(beta), function(x) as.expression(substitute(beta==C,list(C=as.name(x)))))) )
    
    })
  
  output$beta_mc <- renderPrint({
    R <- input$R_mc
    if(input$boot_mc != 1){ 
    
    alpha <- input$sliderAlpha
    beta<-input$sliderBeta
    R <- input$R_mc
    n <- input$n_mc
    
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
        
        # suppressWarnings(sum(log(dbeta(x,shape1=alpha,shape2=beta))))
        (-n*lbeta(alpha,beta) + (alpha-1)*sum(log(x)) + (beta-1)*sum(log(1-x)))
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
    means<-colMeans(result)
    true_par<-rep(c(alpha,beta),2)
    
    #vies
    bias<-means-true_par
    
    #VR
    RB<- (bias/true_par)*100
    
    #erro padrao
    SE <- apply(result,2,sd) 
    
    #EQM
    MSE <- apply(result,2,var)+bias^2 
    
    # final
    all<-cbind(true_par,means,bias,RB,SE,MSE)
    
    # cat("\n")
    # print(paste("Contador Falhas:",falhas),q=F)
    # cat("\n","Results of Monte Carlo simulation for the point estimation of the beta distribution parameters","\n")
    cat("n =",n,"\n")
    cat("R =",R,"\n")
    rownames(all)<-c("alpha MM","beta  MM","alpha MV","beta  MV")
    print(round(all,3))
    cat("\n")
      } #else { cat("Loading bootstrap...")}
  })
  
  
  output$beta_mc_boot <- renderPrint({
    if(input$boot_mc){ 
      
      alpha <- input$sliderAlpha
      beta<-input$sliderBeta
      R <- input$R_mc
      n <- input$n_mc
      B= input$B_mc_boot
    
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
        
        # suppressWarnings(sum(log(dbeta(x,shape1=alpha,shape2=beta))))
        (-n*lbeta(alpha,beta) + (alpha-1)*sum(log(x)) + (beta-1)*sum(log(1-x)))
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
    
    withProgress(message="Loading bootstrap...", value=0,{
      
    
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
        
        x_boot_mv <- rbeta(n,shape1=fit$mv[1],shape2=fit$mv[2])
        fit_boot_mv<-betafit(x_boot_mv)
        resultBmv<-rbind(resultBmv,c(fit_boot_mv$mm,fit_boot_mv$mv))
        
      }
      # if( (100*i/R)%%20==0) print(c((100*i/R),"% --------------- R ---------------"), quote = F)
      incProgress(1/R, detail = paste(round(100*i/R,1),"%"))
      Sys.sleep(0.1)
      
    }
    
    })
    
    
    #MM corrigido por MV
    boot_emmA_mv<- 2*mean(result[,1]) - mean(resultBmv[,1]) #take care bro
    boot_emmB_mv<- 2*mean(result[,2]) - mean(resultBmv[,2])
    
    
    #MV corrigido por MV
    boot_emvA_mv<- 2*mean(result[,3]) - mean(resultBmv[,3])
    boot_emvB_mv<- 2*mean(result[,4]) - mean(resultBmv[,4])
    
    # medias_boot_mm <-c(boot_emmA_mm,boot_emmB_mm, boot_emvA_mm,boot_emvB_mm)
    medias_boot_mv <-c(boot_emmA_mv,boot_emmB_mv, boot_emvA_mv,boot_emvB_mv)
    
    means<-c(colMeans(result),medias_boot_mv)
    true<-rep(c(alpha,beta),4)
    
    #vies
    bias<-means-true
    
    #VR
    RB<- (bias/true)*100
    
    #erro padrao
    SE <- apply(result,2,sd) 
    
    #EQM
    MSE <- apply(result,2,var)+bias^2 
    
    # final
    all<-cbind(true,means,bias,RB,SE,MSE)
    
    # cat("\n")
    # print(paste("Contador Falhas:",falhas),q=F)
    cat("n =",n,"\n")
    cat("R =",R,"\n")
    cat("B =",B)
    cat("\n")
    rownames(all)<-c("alpha MM","beta  MM","alpha MV","beta  MV",
                     "alpha MM*","beta  MM*","alpha MV*","beta  MV*")
    # "alpha MM*_mm","beta  MM*_mm","alpha MV*_mm","beta  MV*_mm"
    print(round(all,3))
    }
  })
  
  
}

# Create Shiny app ----
shinyApp(ui, server)





