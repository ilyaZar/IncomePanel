################################################################################
################ Simultaneous Bayesian Modelling of a Panel of #################
############### Parametric Income Distributions for Grouped Data ###############
############# Tobias Eckernkemper, Bastian Gribisch, Ilya Zarubin ##############
################################################################################
#
#
#
#
#
General To-Do
- Reorder Code in project folder as well as github repo
  --> change locally, then pull to gitHub
  --> must inlcude re-ordering of current code along the following steps 
  1.) data generation
  2.) MCMC with R i.e. some *.R file
  3.) MCMC with STAN i.e. some *.Stan file
  4.) A main file running 2.) and 3) with variuos settings
  5.) Add general *.R file for convergence diagnostics including some shiny app
  command argument: function(..., shiny = TRUE), where shiny = TRUE 
  runs shinyStan, but when false, ... runs nicely output diagnostics
- Devise Project according to 7 steps:
#
#
#
#
#