# SparseBiplots
Sparse Biplots are modern variant of HJ-Biplot which attempts to find sparse loadings, i.e., weight vectors with only a few nonzero values. Introducing a penalty parameter helps improve the interpretation of the model because the axis are formed as a linear combination of only some of the original variables.

## Installation
* Install *SparseBiplots* from **CRAN**:
```R
install.packages("SparseBiplots")
```

* Install *SparseBiplots* from **GitHub** using [devtools](https://github.com/r-lib/devtools):
```R
install.packages("devtools")
devtools::install_github("mitzicubillamontilla/SparseBiplots")
```

## Usage

### Functions

* `HJBiplot`: 
* `Ridge_HJBiplot`: 
* `LASSO_HJBiplot`: 
* `ElasticNet_HJBiplot`: 

### Parameters

### Values

## Example
```R
library("SparseBiplots")
HJBiplot(mtcars, transform_data = 'scale', ind_name  = TRUE)
```
