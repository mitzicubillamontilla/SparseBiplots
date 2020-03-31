# SparseBiplots
*SparseBiplots* library is a modern variant of HJ-Biplot which attempts to find sparse loadings, i.e., weight vectors with only a few nonzero values. 

## Installation
* Install *SparseBiplots* from **CRAN**:
```{r}
install.packages("SparseBiplots")
```

* Install *SparseBiplots* from **GitHub**:
```{r}
# install.packages("devtools")
devtools::install_github("mitzicubillamontilla/SparseBiplots")
```

## Usage

## Example
```{r}
library("SparseBiplots")
HJBiplot(mtcars, transform_data = 'scale', ind_name  = TRUE)
```
