# Covid-model

The covid vaccine project contains the code enabling to produce the results
presented in the article entitled "Potential impact of introducing vaccines against COVID-19 under supply and uptake constraints in France: a modelling study" article 

The project is organized as follows:
- code
  + model_function.R: the SEIR transmission model itself, with utility functions,
                   to load inputs or plot graphs for example, and general
                   parameters
  + compare_incid.R: the code to generate graph comparing different scenarios.
  + graph_hosp.R: the code to generate hospitalisations graph for a range of vaccine efficacy.
- input : this folder contains all the input files used for generating the different scenarios.
- result: The folder to store result files corresponding to each scenario
- figure: The folder to store graph allowing comparing scenario results
              

## Installation & execution

The code was developed and tested with R version 3.6.2 (2019-12-12) through
R Studio. 

### Dependencies
Two libraries need to be present on your R system:
- deSolve
- tidyverse
- scales


### Execution

1. Open the RStudio project
2. Launch the R notebook "Results_vaccination"
3. You will find newly produced results in a Results directory sub-folder and Figures in figure sub-folder

## Disclaimers

- The code is provided to ensure transparency and reproductability of graphs and
  figures exposed in the article. No specific support for someone willing to reuse this
  model in another context will be provided.
- As with any mathematical model, misconfiguring inputs can lead to
  misinterpretation and/or meaningless outputs. The development team only
  can only endorse outputs it has itself generated.
- We have developed and checked thoroughly this code. Nevertheless, an error is
  possible. We would be grateful if you consider contacting us, should you find
  any errors.

## Acknowledgments
The code presented here expands on a code initially developped for the following article https://dx.doi.org/10.1016%2Fj.vaccine.2020.10.034

-> Gibub version