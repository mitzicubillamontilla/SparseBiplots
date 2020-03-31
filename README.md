# SparseBiplots
*SparseBiplots* library is a modern variant of HJ-Biplot which attempts to find sparse loadings, i.e., weight vectors with only a few nonzero values. 

## Installation
* Install *SparseBiplots* from **CRAN**:
```R
install.packages("SparseBiplots")
```

* Install *SparseBiplots* from **GitHub**:
```R
# install.packages("devtools")
devtools::install_github("mitzicubillamontilla/SparseBiplots")
```

## Usage

## Example
```R
library("SparseBiplots")
HJBiplot(mtcars, transform_data = 'scale', ind_name  = TRUE)
```
