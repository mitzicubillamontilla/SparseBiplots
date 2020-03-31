# SparseBiplots
Sparse Biplots are modern variant of HJ-Biplot which attempts to find sparse loadings, i.e., weight vectors with only a few nonzero values. Introducing a penalty parameter helps improve the interpretation of the model because the axis are formed as a linear combination of only some of the original variables.

This package provides functions to perform the HJ-Biplot (Galindo, 1986) and modifications introducing Ridge, LASSO and Elastic Net penalty:

* `HJBiplot`: performs the HJ-Biplot .
* `Ridge_HJBiplot`: performs the HJ-Biplot introducing Ridge penalty.
* `LASSO_HJBiplot`: performs the HJ-Biplot introducing LASSO penalty.
* `ElasticNet_HJBiplot`: performs the HJ-Biplot introducing Elastic Net penalty. For this use the [spca](https://github.com/erichson/spca) package.

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

### Parameters

* `X`: data frame which provides the data to be analyzed. All the variables must be numeric.
* `Transform_data`: a value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. For default is "scale".
* `ind_name`: if it is TRUE it prints the name for each row of X. If it is FALSE (default) does not print the names.
* `vec_name`: 

### Values

## Example
```R
library("SparseBiplots")
HJBiplot(mtcars, transform_data = 'scale', ind_name  = TRUE)
```
