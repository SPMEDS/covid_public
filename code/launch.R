#####################################################

scen = 'IT1'

tmonths <- 18  # total number of months  for the simulation

per_deb <- 1  # start period for graphs
per_fin <- 18 # end period for graphs  

datebreak = "2 months"
fobs = "IT_observed"

source('code/covid_model.R')


#####################################################

scen = "US1"

tmonths <- 18  # total years for the simulation

per_deb <- 1  # start period for graphs
per_fin <- 18 # end  period for graphs  

datebreak = "2 months"
fobs = "US_observed"


source('code/covid_model.R')

#####################################################

scen = 'KOR1'

tmonths <- 18  # total number of months  for the simulation

per_deb <- 1  # start period for graphs
per_fin <- 18 # end period for graphs  

datebreak = "2 months"
fobs = "KOR_observed"

source('code/covid_model.R')

