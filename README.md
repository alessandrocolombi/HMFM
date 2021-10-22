# GDFMM
R package for Group Dependent Finite Mixture Models

# R/
This folder contains all .R files. Do not modify GDFMM-package.R and zzz.R (used to load and unload the dynamic c++ library created from source files). Add you R functions to GDFMM_export.R. Functions has to be documented using Roxygen. Roxygen commands starts with #'. Add #' @export at the end to include the function in the NAMESPACE, those are the functions made available to the user, the one you use in the package. RcppExports.R can not be modified by hand.

# src/
This folder contains all c++ and Rcpp files, the ones to be compiled. There are three types of files:
1. Makevars and Makevars.win are sort of makefiles for linux and window, respectively.
2. GDFMM_exports.cpp is the one where you write the functions that has to be made available. Use Roxygen to document those functions, in this case Roxygen comments start with //' and not with #' as in .R file. There are you ways to make a function available. You may only add // [[Rcpp::export]], this means that the function can be called from R (in particular it can be called from the R function you write in GDFMM_export.R!) but it is not added to the NAMESPACE file. As a consequence, in principle the user may call those functions, but since they are not included in the NAMESPACE nobody knows that those functions exist and are never called. They are kind of hidden from the user. Otherwise you can add both //' @export and // [[Rcpp::export]] to export the function and include it in the NAMESPACE, so that it can be called from the user. See EigenTest() and TestGSL() as examples (even if they are not in this file, the idea is the same). RcppExports.cpp file is automatically generated and can not be modified by hands.
3. Test_Eigen.cpp and Test_GSL.cpp are just examples to be runned. You may think of them as the main function one usually call when running c++ code. They take no input and do something. 
4. All other files are you c++ code. You can add the header files (**Remark** use .h extension, not .hpp) and the source files (use .cpp extension) that define your classes. Try to avoid Rcpp function calls inside those functions and comment them as you normally do in c++.

# man/
This folder contains the Rmarkdown files for the documentation. Do not edit by hand.

**The workflow**
TODO
