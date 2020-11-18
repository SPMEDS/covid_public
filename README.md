# Covid-model

The _covid-model_ project contains the code enabling to reproduce the results
presented in the article entitled "Critical uncertainties impacting accelerated industry response to 
COVID-19: modeling insights Age-structured" article 

The project is organized as follows:
- Code
  + covid_model.R: the SEIR transmission model itself, with utility functions,
                   to load inputs or plot graphs for example, and general
                   parameters
  + launch.R: code to launch the computation and creation of figures for three
              use cases: the mid-term medium thresholds scenario for Italy, United
              States and South Korea.
- Input : this folder contains the configuration settings specific to each
          scenario.
  + IT1: folder for Italy
  + KOR1: folder for South-Korea
  + US1 : folder for the United States of America
- Observed: reported cases and deaths, from European Center for Disease
            Prevention and Control (https://qap.ecdc.europa.eu/public/extensions/COVID-19/COVID-19.html),
            and adjusted cases per day (method described in the Supplementary material)
- Results: folder with one sub-folder per country including both csv files containing 
           model outputs by day and age group and png files being the output graphs.
           The files provided correspond to the three scenarios described above.
              
It is possible to generate the other graphs in the manuscript by adapting
parameters in launch.R script and input files based on information
detailed in the Supplementary material.

## Installation & execution

The code was developed and tested with R version 3.6.2 (2019-12-12) through
R Studio.

### Dependencies
Two libraries need to be present on your R system:
- deSolve
- tidyverse

### Execution

1. Ensure the working directory for R is "~/your_own_path/covid-model"
2. Launch the execution of part or whole launch.R
3. You will find newly produced results in a Results directory sub-folder(s)

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
The development of this code benefitted from code implemented for a malaria transmission model at the University of Cape Town, South Afric
https://github.com/sheetalsilal/METCAP