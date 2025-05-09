Read me file for the article entitled 
"Neighbor density-dependent facilitation promotes coexistence and internal oscillation"
by
Lisa Buche, Lauren Shoemaker, Lauren Hallett, Peter Vesk, Oscar Godoy and Margie Mayfield 
Published in 
""

Buche L., Shoemaker L. et al (2025), Buchel9844/Density-dependentFunctions: Initial release (V.0). Zenodo. DOI: 10.5281/zenodo.14618455  (version V.1 will be released upon acceptance).
Contact details: Lisa Buche (buchel9844@gmail.com or lbuche@student.unimelb.edu.au)

Statement of authorship: All co-authors designed the study. Margaret Mayfield, Oscar Godoy, Lauren Hallett and Lauren Shoemaker obtained the necessary funding. Lisa Buche built the framework and performed the simulated and empirical analysis with substantial input from Lauren Shoemaker, Oscar Godoy, and Peter Vesk. Members of the Mayfield lab, including Lisa Buche, collected the data. Lisa Buche wrote the manuscript with substantial input from Lauren Shoemaker and Margaret Mayfield. All co-authors contributed to the review of the manuscripts. 

Abstract: 
The ability of species to form diverse communities is not fully understood. Species are known to interact in various ways with their neighborhood. Despite this, common models of species coexistence assume that per capita interactions are constant and competitive, even as the environment changes. In this study, we investigate how neighbor density-dependent variation in the strength and sign of species interactions changes species and community dynamics. We show that by including these sources of variation, predictions of ecological dynamics are significantly improved compared to outcomes of typical models that hold interaction strengths constant. We compared how well models based on different functions of neighbor density and identity did in describing population trajectories (i.e., persistence over time) and community dynamics (i.e., temporal stability, synchrony and degree of oscillation) in simulated two-species communities and a real diverse annual plant system. In our simulated communities, we found the highest level of coexistence between species pairs when species interactions varied from competitive to facilitative according to neighbor density (i.e., following a sigmoid function). Introducing within-guild facilitation through a nonlinear bounded function allowed populations, both simulated and empirical, to avoid extinction or runaway growth. In fact, nonlinear bounded functions (i.e., exponential and sigmoid functions) predicted population trends over time within the range of abundances observed over the last 10 years. With the sigmoid function, the simulated communities of two species displayed a higher probability of synchrony and oscillation than other functional forms. These simulated communities did not always show temporal stability but were predicted to coexist. Overall, varying species interactions lead to realistic ecological trajectories and community dynamics when bounded by asymptotes based on neighbor density. These findings are important for advancing our understanding of how diverse communities are sustained and for operationalizing ecological theory in the study of the real world.

Authorship of data: Members of the Mayfield lab. We want to thank John Dwyer, Claire Wainwright, Maia Raymundo, Trace Martyn, Victoria Reynolds, Catherine Bowler, Aubrie James, Abigail Pastore, Manuel Sevenello, and Courtney Taylor for the data collected in the Perenjori region. 

Authorship of code: R Code was written by Lisa Buche.


Details of R script and their function:

1.ComputeSim.R - Simulated the communities within the parameters' ranges based on the Ricker model and the functional forms which are written as fucntion in  "1.1PopProjection_tollbox.R". This script runs a loop to simulated the time series of the 1500 communities. The script also gives multiple ways to visualise the simulated communities and create Fig 2 and Fig S1 in the Appendix (some of which are not included in the manuscript).  


2.Comcomposition.R - Extract the population trajectory of the simulated communities. It determines the results to create Fig 5 in the maintext and Fig S2 in the Appendix. 

3.StabilityAssesmentR - Extract the community dynamics of the simulated communities based on the function of the script "3.1.TimeSerie_tollbox.R". It determines the results to create Fig 7 in the maintext and Fig S3 in the Appendix. 

3.1.TimeSerie_tollbox.R was modfied from the script from Ives et al, 2010 - Analysis of ecological time series with ARMA(p,q) models. Anthony R. Ives, Karen C. Abbott, Nicolas L. Ziebarth. 01 March 2010 https://doi.org/10.1890/09-0442.1


4.SensitivityAnalysis.R - Perform the GLM that test how the paramters influence the communities stability when they are predicted to coexist. It creates the results and figure of Box 1 and Fig S4-9 in the appendix.

5.NatureDataFit_Perenjory.R - Bayesian fit of the case study for LARO based on the Stan code in " DensityFunct_Final.stan". This script fit the competition data and projects it based on the family density projection. It creates the data for Fig 3-4, 6 in the maintext and Figure S11 in the appendix.

6.CheckModelBehavior.R - Extracts the models' convergence, Posterior predictive distribution of LARO’s fecundity (Fig S12) and explanatory power of the final fitted models, which were evaluated according to Root mean squared deviance and the leave-one-out approximation (see Appendix section "Model behaviour check" and section "Model comparison"). 

6.1. Mean.abundance.sp.R - Compute the mean abundances observed in the neighbourhood of a specific focal for a specific year of individuals within and across trophic levels. 

