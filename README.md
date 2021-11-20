# GDFMM
R package for Group Dependent Finite Mixture Models

# R/
This folder contains all .R files. Do not modify GDFMM-package.R and zzz.R (used to load and unload the dynamic c++ library created from source files). Add you R functions to GDFMM_export.R. Functions has to be documented using Roxygen. Roxygen commands starts with #'. Add #' @export at the end to include the function in the NAMESPACE, those are the functions made available to the user, the one you use in the package. RcppExports.R can not be modified by hand.

# src/
This folder contains all c++ and Rcpp files, the ones to be compiled. There are three types of files:
1. Makevars and Makevars.win are sort of makefiles for linux and window, respectively.
2. GDFMM_exports.cpp is the one where you write the functions that has to be made available. Use Roxygen to document those functions, in this case Roxygen comments start with //' and not with #' as in .R file. There are you ways to make a function available. You may only add // [[Rcpp::export]], this means that the function can be called from R (in particular it can be called from the R function you write in GDFMM_export.R!) but it is not added to the NAMESPACE file. As a consequence, in principle the user may call those functions, but since they are not included in the NAMESPACE nobody knows that those functions exist and are never called. They are kind of hidden from the user. Otherwise you can add both //' @export and // [[Rcpp::export]] to export the function and include it in the NAMESPACE, so that it can be called from the user. See EigenTest() and TestGSL() as examples (even if they are not in this file, the idea is the same). RcppExports.cpp file is automatically generated and can not be modified by hands.
3. Test_Eigen.cpp and Test_GSL.cpp are just examples to be runned. You may think of them as the main function one usually call when running c++ code. They take no input and do something. 
4. Put in include_headers.h all the `#include<lib>` for the c++ internal libraries. Do not include here custom classes, i.e `#include "MyClass.h" ` .
5. recurrent_traits.h contains all the typedefs. Use custom and significant names for your types, use typedef as much as you can.
6. Use utils.h to define custom free functions, i.e those operations that are not class methods.
7. Add all the c++ files you want. You can add the header files (**Remark** use .h extension, not .hpp) and the source files (use .cpp extension) that define your classes. Try to avoid Rcpp function calls inside those functions and comment them as you normally do in c++.

# man/
This folder contains the Rmarkdown files for the documentation. Do not edit by hand.

# Getting started
The first thing to be done is to clone this repository in a local directory. Open Git or Git Bash and set the directory where you want the folder to be cloned into. The run 
```shell
$ git clone https://github.com/alessandrocolombi/GDFMM.git
```
or run 
```shell
$ git clone git@github.com:alessandrocolombi/GDFMM.git
```
if you want to use the SSH protocall (prefered but requires having an ssh key). From now on, you have a local git repo of the package. Run the GDFMM.Rproj file to open the R project in Rstudio. Next to the environment button you should see two additional buttons, Build and Git. 

# How to use Git
There are three important actions in Git. Pull, commit and Push. The fourth main action is the Merge, the one we want to avoid because it is needed to solve conflicts. There are two levels of repositories: your local repo and the shared one, that is the one you see in the GitHub page.
1. **pull:** using this action you update your local repo by "downloading" everything that is available in the shared repo. Fox example, one has modified a function. By pulling you update your local file. **It is good habit to pull everytime you start working**, most of the conflicts that rises are due to the fact that you forget to pull.<br/>
How to pull? You the blu arrow you see in Rstudio in the Git window. The corresponding command is 
```shell
$ git pull
```
2. **commit:** the idea of the commit action is to "save safely" your local changes, but do not share them with the others. For example, you change something in file.R and save the file (i.e you save the R file the usual way). By doing so, the file has been saved but it is local to your pc, if something happens, you may lose it. By committing a file, you save it in the git (but not in github) file manager system. Your are able to recover it whenever you want. But not only, if some times later your are not satisfied by your change, you can recover the older version of the file. In other words, by committing files it is very difficult to mess up you work. That is why it is considered a "safe save". <br/> How to commit? When you create a new file or modify an existing file, the name of the file appears in the Git window in R. Check the box next to the files you want to commit and then click the commit button. Note that every commit has to be associated to a meaningful message, so that your are able to remeber what was the stage of your work at that time. The corresponding commands are
```shell
$ git add --all
$ git commit -m "message"
```
By doing so, you commit all the modified files. If you want to commit only specific files, you can substitute --all with the name of the files.
3. **push:** this is the dual action of pull. You push the previously committed files to the shared repo in github. Note that you always have to commit the files before pushing them. When you push files, everyone else is able to "download" them using the pull command.<br/>
How to push? You the green arrow you see in Rstudio in the Git window. The corresponding command is 
```shell
$ git push
```
4. **status:** Using the following command
```shell
$ git status
``` 
you can analyze the differences between the status of your local repo and the shared one, also called master. If the master is ahead, you need to pull. If your repo is ahead, you can push. The status keeps track of the modified files, check it frequently. **Always check the status before using pull, commit and push!**<br/>
5. **merge:** if you are not careful, conflicts may happen. The most common case is when you modify some files before pulling. Or when you push without a previous pull. In those cases, you need to solve the conflict by editing files (usually this is done by hand).<br/>

A me git non piace (perché non lo so usare bene). Io di solito mi tengo i file salvati normalmente e quando raggiungo uno stato che mi piace faccio commit di tutto e immediatamente dopo il push. Non è una sana abitudine perché si perdono appunto i benefici del commit che ho scritto sopra, però così evito di fare casini. Spero siate più bravi di me, usate git come preferite ma usatelo, non tenete i file salvati in qualche vostra cartella condivisa tipo OneDrive. Siete in 5 e ognuno modificherà parti e file diverse, se state attenti e lo usate facendo operazioni semplici e con attenzione, git è uno strumento molto potente che vi impedisce di fare casini ed è integrato molto bene con il sistema di pacchetti di R. 

# The workflow

## R package code organization
What are the benefits of an `R` package with compiled code? This is a nice solution to get the flexibility and simplicity of the `R` language to manipulate data, table and to visualize plots, without giving up efficient code for the computationally intense operations.<br/>
The idea is to write in c++ the classes needed to run the sampler. It may useful to use write c++ code also to compute posterior expectations and posterior manipulations, this  may depend on the specific problem. Everything else can be written in R. Use the src/ folder to add all the c++ files where you define the classes needed to define the sampler. You may think to this folder as an independent c++ library (actually, a dynamic c++ library is created). Let us assume that you created a complex c++ class called GDFMM_sampling_strategy. It is unlikely that the inputs needed by such a class are simple built-in types but it probably needs to get other c++ custom classes. In this example, you need to pass the classes where you defined the likelihoos, the priors, the clustering mechanism and so on.
However, an package is effective if the user can run the sampler in the simplest possible way. Usually one has just a matrix (or data.frame) where data are stored and some input parameters (number of iterations, hyperparameters, initial values). This mean that it is not straightforward to write an R function that gets such inputs and is able to call the run method define inside the GDFMM_sampling_strategy class. We need an intermediate step. This is where the Rcpp sytax kicks in.
In the src/ folder, more precisely in the GDFMM_exports.cpp file, we define a custom Rcpp function that is a bridge from R data to c++ types. An Rcpp function is a c++ function (you are using c++, not R) that uses the Rcpp library which defines many types that are very similar to the R syntax (for example Rcpp::List). Such a language is created to let R and c++ to comunicate. The workflow is the following:
1. Create an R function that takes the data. Pre-process data with simple operations that can be done in R.
2. Inside such a function, call the Rcpp function that takes the pre-processed data from R and create all the c++ object and classes you need. Within this function, you have all you need to create a GDFMM_sampling_strategy, run its main method, collect the result and return them in R.<br/>

The workflow presented above is valid for those functions that require complex or custom c++ objects as input. As explained in [src/](#src/), you may also have c++ functions that takes "simple" c++ objects as input. Those objects that can automatically converted from R language to c++ thaks to the Rcpp package, such as: int, double, strings, lists and matrices. If so, it is enough to implement the Rcpp function without creating the R counterpart. Finally, remember that most of the data manipulation operations can be implemented in R, without speeding up the code using c++.  <br/>

A couple of examples are reported.

## Compile, document, build and run
Whenever you want to test or check your code, you need to compile it. As we are dealing with an R package, this step involves several substeps. Fortunately, most of the operations happens under the hood. Most of the times, you compile your code when you are inside the GDFMM.Rproj. If so, run  <br/>
`ctrl+shift+D` and then `ctrl+shift+B` . <br/>
The first command compile the c++ code (it creates the GDFMM.so library actually) and update the documentation. The second one simulate the building of the R package. This second part requires a deeper knowledge of R packages to be properly understood, in practice what it does is to compile the code again and install the package. Of course, if no files is modifies after running `ctrl+shift+D`, there is nothing to be compiled and this second step is much faster. Once it is completed, you have a test version of your package installed, you may check it in the package window. It works as all other R packages, you can check the available functions, the documentation and run all the exported functions (those whose name appears in the NAMESPACE at least). You can test whatever you want.<br/>
**Remark:** the c++ compiler you have in R in not as friendly as the one you have in a pure c++ environments. Error messages are ugly and long. Use the usual rule to trust and fix only the first one that is printed. Moreover, R is very bad in failure handling. That is, runtime errors such as segmentation fault or similar do not give you meaningful messages but results in R abort (that is why R session aborts are so frequent!) If so, be patient when debugginig.

## Install the package
When the work is done, or whenever you want to install packages that have been published only on github and on the CRAN, you can install the package. The difference with respect to the previous section is that this step is not done when you are inside the GDFMM.Rproj. Therefore it is not local to the project but you make the package available in all your R system.
There are several ways to install a github package. The procedures presented here are valid in general, not only for GDFMM.
1. clone the repo as explained in [# Getting started](# Getting started). Open R and set the repo as working directory (select the folder where there is the .Rproj file). Run
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
from command line. (Questo funziona in linux, non ho mai creato il tar.gz in window). Once that the file is sent to the user, the latter just need to open R, set the working directory to match the file location and run
```R
install.package("package_name.tar.gz", repos = NULL )
```
# Guidelines
1. Siete liberi di sviluppare il codice come preferite. Se quando fate esperienza con R imparate un modo diverso per chiamare le funzioni ben venga. L'importante è che per l'user sia tutto il più semplice possibile. Serve flessibilità nel passare i dati in forme diverse e pochi input da decidere.
2. Commentate per bene tutto il codice! Non solo all'inizio della funzione, idealmente per ogni scope serve una breve descrizione che dice cosa state facendo.
3. Documentate bene le funzioni che esportate.
4. Usate nomi significativi, usate tante typedef e definitele nel file che vi ho indicato sopra.
5. Le regole di buon coding dicono di usare sostantivi come nomi degli oggetti e verbi per le funzioni. Se riuscite, meglio.







