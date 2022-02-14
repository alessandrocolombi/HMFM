# GDFMM
R package for Group Dependent Finite Mixture Models

# R/
This folder contains all .R files. Do not modify GDFMM-package.R and zzz.R (used to load and unload the dynamic C++ library created from source files). In GDFMM_export.R there are all functions implemented for the library, so those functions that ca be called from R. All functions are documented with Roxygen 

# src/
This folder contains all c++ and Rcpp files, the ones to be compiled. There are three types of files:
1. Makevars and Makevars.win are sort of makefiles for linux and window, respectively.
2. **GDFMM_exports.cpp** is the file to export C++ functions in R environment with Rcpp. All functions are documented with Roxygen
3. **Test_Eigen.cpp** and **Test_GSL.cpp** are just example functions that test if Eigen and GSL libraries are installed correctly. Moreover this functions thest if user-defined classes to sample from GSL library works correctly 
4. **include_headers.h** contains all the `#include<lib>` for the C++ internal libraries. Custom classes, i.e `#include "MyClass.h" `, are not included in this file.
5. **recurrent_traits.h** contains all the typedefs for values or data structures that are recurrently called
6. **utils.h** contains custom free functions, i.e those operations that are not class methods.
7. **FullConditional.h** is a pure virtual class from wich are derived all classes for Full Conditionals, i.e. files **FC_*.h** and **Partition.h**
8. **GS_data.h** and out_data.h are data structure to contain data needed in a Gibbs Sample iteration and for the output of GDFMM_sampler, respectively.
9. **GibbsSampler.h** class for a Gibbs Sampler object. Gives the possibility to specify if the partition is fixed and the prior distribution for the components
10. **RcppExports.cpp** automatically generated, do not edit by hand

# man/
This folder contains the Rmarkdown files for the documentation. Do not edit by hand.

# tests/
This folder contains R scripts to test the GDFMM algorithm on simulated datasets or on the Mathscore dataset (Hoff 2009)

# Compile, document, build and run
The first thing to be done is to clone this repository in a local directory. Open Git or Git Bash and set the directory where you want the folder to be cloned into. The run 
```shell
$ git clone https://github.com/alessandrocolombi/GDFMM.git
```
or run 
```shell
$ git clone git@github.com:alessandrocolombi/GDFMM.git
```

To compile the package from Rstudio you have to open the Rproject GDFMM.Rproj.
Once you did it to create the library and use all functions implemented you need to run <br/>
`ctrl+shift+D` and then `ctrl+shift+B` . <br/>
The first command compile the c++ code (it creates the GDFMM.so library actually) and update the documentation. The second one simulate the building of the R package.

## Install the package
When the work is done, or whenever you want to install packages that have been published only on github and on the CRAN, you can install the package. The difference with respect to the previous section is that this step is not done when you are inside the GDFMM.Rproj. Therefore it is not local to the project but you make the package available in all your R system.
There are several ways to install a github package. The procedures presented here are valid in general, not only for GDFMM.
1. clone the repo as explained in [# Compile, document, build and run]. Open R and set the repo as working directory (select the folder where there is the .Rproj file). Run.
```R
devtools::install()
```
2. Sometimes one does not want to clone the repo of the package. If so, just open R and run
```R
devtools::install_github("author_github_name/package_name") #for example devtools::install_github("alessandrocolombi/GDFMM")
```
The drawback for those cases, is that the repo has to be public or at least you have access to it.<br/>
3. Another possibility, is that the creator of the package builds the tar.gz file of the package, send it to the user which can install it. The tar.gz file can be created by running 
```shell
$ cd package_repo
$ R CMD build package_name
```
from command line. (This works in Linux, never tested that on Windows). Once that the file is sent to the user, the latter just need to open R, set the working directory to match the file location and run
```R
install.package("package_name.tar.gz", repos = NULL )
```







