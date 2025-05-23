#
# Shiny app build by Lisa Buche - April 2023 

#This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Install the packages below if you do not have these already installed.
#install.packages('tidyverse')
#install.packages('medicaldata')
#install.packages('shiny')
#install.packages('mlbench')
#install.packages('plotly')
#install.packages('ggforce')
library(tidyverse)
library(medicaldata)
library(shiny)
library(mlbench)
library(shiny)
library(plotly)
library(ggplot2)
library(ggpubr)
library(rmarkdown)
library(knitr)
library(pander)
library(ggforce)
#install.packages('rsconnect')
library(rsconnect)

setwd("/Users/lisabuche/Documents/Projects/Density-dependentFunctions/shiny")
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Neighbours density-dependent effect on Fecundity - Pairwise example"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      h3("Values of parameters:"),width=3,
      sliderInput("Ainit",
                  label =withMathJax("Initial effect of j on i when j = 0: \\(\\alpha_{0,i,j}\\)"),
                  min = -1,
                  max = 1,
                  value = -0.1,
                  step=0.01),
      sliderInput("Aslope",
                  label =withMathJax("Decay of interactions from positive or low to higher competition: \\(\\alpha_{ij}:\\)"),
                  min = -1,
                  max = 1,
                  value = -0.8,
                  step=0.01),
      sliderInput("No",
                  label =withMathJax("Optimal density of j that maximises fecundity of i: \\(N_{o,j}\\)"),
                  min = 0,
                  max = 10,
                  value = 5,
                  step=1),
      sliderInput("Ni",
                  label =withMathJax("Density of \\(N_{i}\\)"),
                  min = 0,
                  max = 10,
                  value = 0,
                  step=1),
      #sliderInput("m",
      #           label ="value of m:",
      #           min = 0,
      #          max = 1,
      #          value = 0.5,
      #         step=0.01),
      sliderInput("C",
                  label ="Stretching parameter C:",
                  min = -1,
                  max = 1,
                  value = -0.2,
                  step=0.01),
      sliderInput("RangeNj", 
                  label = withMathJax("Range of Neighbour density \\(N_{j}:\\)"),
                  min = -10, max = 50, value = c(0, 20),
                  step = 1),
      sliderInput("lambda",
                  label ="Lambda value (intrinsic fecundity):",
                  min = 0,
                  max = 10,
                  value = 5,
                  step=1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(   h3("Community dynamics:"),width=9,
                 br(),
                 fluidRow(
                   splitLayout(cellWidths = c("50%", "50%"), 
                               column(5,
                                      uiOutput("ComDyn")),
                               column(4,
                                      uiOutput("Conditions"),
                                      radioButtons("functiontype", h4("Select function type to display fecundity:"),
                                                   choices = list("Constant function" = 1, "Linear function" = 2,
                                                                  "Exponential function" = 3,"Sigmoidal function" = 4),
                                                   selected = 1)
                               )
                   )
                 ),
                 br(),
                 fluidRow(
                   splitLayout(cellWidths = c("50%", "50%"), 
                               tags$b("Graphs of the 4 functions:"),
                               tags$b("Graphs of fecundity of species i in fucntion of Nj:"))
                 ),
                 fluidRow(
                   splitLayout(cellWidths = c("50%", "50%"), 
                               plotOutput("plotfunctions"), 
                               plotOutput("Fecundityplot"))
                 ),
                 br()
                 
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  alpha_function2 <- function(Amin, Aslope,N,No){
    alpha = Amin + Aslope*(N-No)
    return(alpha)
  }
  half.signoidal.function <- function(Amin, Aslope,C,N,No){
    alpha = Amin+C*(1-exp(Aslope*(N-No)))
    return(alpha)
  }
  signoidal.function <- function(Amin, Aslope,C,N,No){
    e = exp(Aslope*(N-No)) # c is stretching the graph horizontally 
    a = C*(1-e) #stretching the graph vertically
    d = Amin
    alpha = (a/(1 + e)) + d
    
    return(alpha)
  }
  functions.plot <- function(Amin, Aslope, C ,N,No){
    
    f1 <- ggplot(data.frame(N= N, alpha = Amin), aes(y=alpha, x= N))+
      geom_line(aes(colour = after_stat(y < 0)))+      
      guides(color="none") + labs(title="function 1") + 
      geom_hline(yintercept=0,linetype="dashed") +
      xlab("Neighbour density of j") + ylab("Resulting effect of j on i") +
      theme_bw() +
      if(Amin <0){
        scale_colour_manual(
          values = c("red")) 
      }else{   scale_colour_manual(
        values = c("limegreen")) }
    
    
    f2 <- ggplot(data.frame(N= N, alpha = alpha_function2(Amin,Aslope,N,No)), 
                 aes(y=alpha, x= N)) +
      geom_line(aes(colour = after_stat(y < 0)))+     
      guides(color="none") + labs(title="function 2") + 
      geom_hline(yintercept=0,linetype="dashed") +
      xlab("Neighbour density of j") + ylab("Resulting effect of j on i")+
      theme_bw() +
      if(all(alpha_function2(Amin, Aslope,N,No) < 0)){
        scale_colour_manual(
          values = c("red")) 
      }else{scale_colour_manual(
        values = c("limegreen","red")) }
    
    
    f3 <- ggplot(data.frame(N= N, alpha =  half.signoidal.function(Amin, Aslope,C,N,No)), 
                 aes(y=alpha, x= N))+
      #geom_smooth(alpha=0.8,color="grey") +
      labs(title="function 3") + 
      geom_hline(yintercept=0,linetype="dashed") +
      geom_point(alpha=0.5,shape=20,aes(colour = after_stat(y < 0)))+    
      guides(color="none") +
      xlab("Neighbour density of j") + ylab("Resulting effect of j on i")+
      theme_bw()+
      if(all( half.signoidal.function(Amin, Aslope,C,N,No) < 0)){
        scale_colour_manual(
          values = c("red")) 
      }else{scale_colour_manual(
        values = c("limegreen","red")) }
    
    f4 <- ggplot(data.frame(N= N , alpha = signoidal.function(Amin, Aslope,C ,N,No)), 
                 aes(y=alpha, x= N))+
      #geom_smooth(alpha=0.8,color="grey") + 
      labs(title="function 4") + 
      geom_point(alpha=0.5,shape=20,aes(colour = after_stat(y < 0)))+     
      guides(color="none") +
      geom_hline(yintercept=0,linetype="dashed") +
      xlab("Neighbour density of j") + ylab("per capita effect of j on i")+
      theme_bw()+
      if(all(signoidal.function(Amin, Aslope,C ,N,No) < 0)){
        scale_colour_manual(
          values = c("red")) 
      }else{scale_colour_manual(
        values = c("limegreen","red")) }
    
    
    
    ggarrange(f1,f2,f3,f4,ncol=2,nrow=2)
    
  }
  fecundity.plot <- function(function.alpha,lambda,
                             Amin, Aslope,C,Nj,Ni,Ni0,Nj0){
    if(function.alpha == 1){
      aii <- 0.2
      aij <- Amin
    }
    if(function.alpha == 2){
      aii <- 0.2 + 0.2*(Ni-Ni0)
      aij <- alpha_function2(Amin, Aslope,Nj,Nj0)
    }
    if(function.alpha == 3){
      aii <- half.signoidal.function(-0.2, -0.2,C ,Ni,No =Ni0)
      aij <- half.signoidal.function(Amin, Aslope,C ,Nj,No =Nj0)
    }
    if(function.alpha == 4){
      aii <- signoidal.function(-0.2, -0.2, C ,Ni,No =Ni0)
      aij <- signoidal.function(Amin, Aslope,C ,Nj, No =Nj0)
    }
    
    
    Fi <- exp(lambda + aii*Ni + aij*Nj)
    return(log(Fi))
  }
  
  output$ComDyn <- renderUI({
    x <- input$Ainit
    y <- input$Aslope
    Ni <- input$Ni
    withMathJax(
      paste0("Equations of community dynamics"),
      paste0("Fecundity distribution: \\( F_i = \\lambda_i * e^{\\alpha_{ii}N_i} * e^{\\alpha_{ij}N_j} \\)"),
      br(),
      paste0("Constant function: \\({\\alpha}_{ij} = \\alpha_{ 0,i,j} \\)"),
      br(),
      paste0("Linear function: \\({\\alpha}_{ij} = \\alpha_{ 0,i,j} + \\alpha_{ij} * (N_j - N_{o,j}) \\)"),
      br(),
      paste0("Exponential function: \\({\\alpha}_{ij} =  \\alpha_{ 0,i,j} + C*(1-e^{\\alpha_{ij}({N}_{j}-N_{o,j})} \\)"),
      br(),
      br(),
      paste0("Sigmoidal function: \\({\\alpha}_{ij} = \\alpha_{ 0,i,j} + \\dfrac{C*(1-e^{\\alpha_{ij}({N}_j- N_{0,j})})}{1+e^{\\alpha_{ij}({N}_j - N_{o,j})}} \\)"),
      br(),
      paste0("Interaction of j on i when neighbourhood density \\({N}_j = 0\\), is equal to (\\(\\alpha_{ 0,i,j}\\))  :", x),
      br(),
      paste0("Interaction of j on i increases per capita by (\\(\\alpha_{J}\\)) :", y),
      br(),
      paste0("Neighbours density of species i (conspecific density) is kept constant (\\(N_{i}\\)) :",Ni),
      br(),
      paste0("Intrapecfic parameters are kept constant, \\( \\alpha_{0,i,i} = -0.2 ; \\alpha_{N_{i}} = -0.2 \\)")
    )
  })
  
  
  output$plotfunctions <- renderPlot({
    Ainit =input$Ainit
    Aslope = input$Aslope
    C = input$C
    No = input$No
    Ni = input$Ni
    functions.plot(Amin =Ainit, Aslope =Aslope,
                   C=C ,N=seq(input$RangeNj[1],
                              input$RangeNj[2],0.25),
                   No=No)
  }, res = 96)
  
  output$Fecundityplot<- renderPlot({
    Amin =input$Ainit
    Aslope = input$Aslope
    C = input$C
    Nj0= input$No
    Ni0= 0  
    Ni = input$Ni
    Nj = seq(input$RangeNj[1],
             input$RangeNj[2],1)
    lambda = input$lambda
    function.alpha= input$functiontype
    
    ggplot(data.frame(Nj = Nj,
                      fecundity =  fecundity.plot(function.alpha=function.alpha,lambda=lambda,
                                                  Amin=Amin, Aslope=Aslope, 
                                                  C =C,Nj=Nj,Ni=Ni,Ni0=Ni0,Nj0 = Nj0)), 
           aes(y=fecundity, x= Nj))+
      #geom_smooth(alpha=0.8,color="grey") + 
      labs(title="Fecundity distribution") + 
      geom_point() +
      xlab("Neighbour density of species j") + ylab("Fecundity of species i (ln)")+
      theme_bw()
  }, res = 96)
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
