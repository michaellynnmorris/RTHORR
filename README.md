
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RTHORR

<!-- badges: start -->

<!-- badges: end -->

The goal of RTHORR is to make it simple to run a text file of
correlation matrices - or raw data - through the randomization test of
hypothesized order relations (RTHOR; Hubert & Arabie, 1987) and make
comparisons of fit between pairs of matrices (or raw dataframes). This
program is especially valuable in the evaluation of circumplex models of
data as found in color perception, vocational interests, and
interpersonal behavior.

Hubert, L., & Arabie, P. (1987). Evaluating order hypotheses within
proximity matrices. Psychological Bulletin, 102(1), 172–178.
<https://doi.org/10.1037/0033-2909.102.1.172>

## Installation

You can install this development version of RTHORR from
[github](https://github.com/) with:

``` r
devtools::install_github("michaellynnmorris/RTHORR")
```

## Examples

``` r
library(RTHORR)
#> Loading required package: gdata
#> gdata: read.xls support for 'XLS' (Excel 97-2004) files ENABLED.
#> 
#> gdata: read.xls support for 'XLSX' (Excel 2007+) files ENABLED.
#> 
#> Attaching package: 'gdata'
#> The following object is masked from 'package:stats':
#> 
#>     nobs
#> The following object is masked from 'package:utils':
#> 
#>     object.size
#> The following object is masked from 'package:base':
#> 
#>     startsWith

#run randall on input.txt
#outputs a single data frame with RTHOR results
randall_output <- RTHORR::randall(n=6,
                                  nmat=3,
                                  ord = "circular6",
                                  input=system.file("extdata", "input.txt", package = "RTHORR"),
                                  description = c("sample_one", "sample_two", "sample_three"))


#run randmf on input.txt
#outputs a list with two dataframes (RTHOR results and comparisons)
randmf_output <- RTHORR::randmf(n=6,
                                nmat=3,
                                ord = "circular6",
                                input=system.file("extdata", "input.txt", package = "RTHORR"))


#view an example of correlation matrix input
# file.edit(system.file("extdata", "input.txt", package = "RTHORR"))
```
