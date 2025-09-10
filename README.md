README
================

<!-- Knit this README.Rmd to generate README.md -->

# CPxplorer

<!-- badges: start -->

<img src="inst/CPxplorer_Logo.png" align="right" height="150px" />

<!-- badges: end -->

### Installation

<!-- You can install the released version of CPxplorer from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("CPxplorer") -->
<!-- ``` -->

To install this R package directly from Github:

``` r
install.packages("devtools")
devtools::install_github("WBS-TW/CPxplorer")
```

After installation, attach the package in RStudio by:

``` r
library(CPxplorer)
```

These functions will be available:  
`CPions()`: opens a shiny app to generate ions of PCAs and analogues.  
`CPquant()`: opens a shiny app to analyze and quantify output from
Skyline.

Server versions is available at shinyapp.io (this is a free tier account
and therefore subjected to monthly limit for the server calculations):  
CPions: <https://wbs-tw.shinyapps.io/CPions/>  
CPquant: <https://wbs-tw.shinyapps.io/CPquant/>

The apps can also be opened directly in a web browser of the local
computer from these sites (no need to install the R package and verified
to work with Chrome):  
HOWEVER: Shinylive is still experimental and there are still some bugs

<https://wbs-tw.github.io/CPions_Shinylive/>  
<https://wbs-tw.github.io/CPquant_Shinylive/>

SOME KNOWN BUGS in Shinylive (but works in the R package and server
versions):  
- CPions_Shinylive: *currently only export to excel the pages in the
panel (not the entire table)*  
- CPquant_Shinylive: *currently Shinylive cannot export all results to
excel, and the user is referred to the other versions for this
functionality*
