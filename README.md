# GDFMM
R package for Hierarchical Mixture of Finite Mixtures

# R/
This folder contains all .R files. Do not modify GDFMM-package.R and zzz.R (used to load and unload the dynamic C++ library created from source files). In GDFMM_export.R there are all functions implemented for the library, so those functions that ca be called from R. All functions are documented with Roxygen 

# src/
This folder contains all c++ and Rcpp files, the ones to be compiled. There are three types of files:
1. Makevars and Makevars.win are makefiles for linux and window, respectively.
2. **GDFMM_exports.cpp** is the file to export C++ functions in R environment with Rcpp. All functions are documented with Roxygen

# man/
This folder contains the Rmarkdown files for the documentation. Do not edit by hand.

## Install the package
Open R and run
```R
devtools::install_github("alessandrocolombi/GDFMM","HMFM") 
```
DEVO RICONTROLLARE L'ULTIMO COMANDO. SECONDO ME C'è UN MODO PER INSTALLARE UN BRANCH SPECIFICO E NON IL MASTER, MA NON SO SE SIA QUELLO






