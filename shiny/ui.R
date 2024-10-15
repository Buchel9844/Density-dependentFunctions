#
# Shiny app build by Lisa Buche - April 2023 

#This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
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
library(rsconnect)

# Define UI for application that draws a histogram
fluidPage(
  
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
                  step=0.1),
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