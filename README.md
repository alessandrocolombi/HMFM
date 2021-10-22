# GDFMM
R package for Group Dependent Finite Mixture Models

# R/
This folder contains all .R files. Do not modify GDFMM-package.R and zzz.R (used to load and unload the dynamic c++ library created from source files). Add you R functions to GDFMM_export.R. Functions has to be documented using Roxygen. Roxygen commands starts with #'. Add #' @export at the end to include the function in the NAMESPACE, those are the functions made available to the user, the one you use in the package. RcppExports.R can not be modified by hand.
