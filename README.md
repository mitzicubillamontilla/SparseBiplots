# SparseBiplots
Este repositorio continen las funciones que se encuentran en el paquete de R *SparseBiplot*

## Installation
* Install `SparseBiplots` from CRAN
```{r}
install.packages("SparseBiplots")
```

* Install `SparseBiplots` from GitHub:
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
