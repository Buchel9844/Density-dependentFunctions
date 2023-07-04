Project to modify pairwise interaction effect to be a function of neighbourhood density. 
Lisa Buche - June 2023

A Rmarkdown is available to understand the story in Summary. Rmd

A shiny App is available to display the 4 different functions of neighbourhood density-dependency. To play with the shiny app - run the "shiny/functions.R" code in a Rstudio session. 

Here is the Rscripts list and what they do.
Code to generate simulated 2 species communities:
GenerateSimData_wrapper.R (to generate the simulated community)
if Stouffer= T - Following Daniel Stouffer model, 2023 - GenerateContinuousData-Stouffer.R; ; GenerateSimData_Stouffer.R (generate the surface response experiment used as our observations dataset)
if Ricker = T - Following a Ricker model - GenerateSimData_Ricker.R
if Malyon= F - Following Malyon Bimler model, 2023 - GenerateSimData_Malyon.R

Code to run Bayesian Model: 
ModelFit.R (to run the model for each function for the Stouffer model)
DensityFunct_APM_Final.stan (stan code for Bayesian model)
ModelFit_Ricker.R (to run the model for each function for the Ricker model)

stan_modelcheck_rem.R (function to check posterior distributions)

Code to make figures and compute likelihood for each model: 
Figure. R 

Code to project population through multiple generations: 
PopProjection.R 
PopProjection_toolbox.R 

