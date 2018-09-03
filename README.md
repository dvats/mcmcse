# mcmcse
An R package for computing Monte Carlo standard
errors (MCSE) in Markov chain Monte Carlo (MCMC) settings. MCSE
computation for expectation and quantile estimators is
supported as well as multivariate estimations. The package also provides 
functions for computing effective sample size and for plotting
Monte Carlo estimates versus sample size.


# Installation
mcmcse can be downloaded directly into R through the the `devtools` package:
```{r}
install.packages("devtools")
library(devtools)
devtools::install_github("statvats/mcmcse")
```
# Citation
Please copy the following into a  `.bib` file.

`
@Manual{fleg:hugh:2017,
title = {mcmcse: Monte Carlo Standard Errors for MCMC},
author = {Flegal, James M and Hughes, John and Vats, Dootika and  Dai, Ning},
year = {2017},
address = {Riverside, CA and Minneapolis, MN},
note = {R package version 1.3-2},
}
`

