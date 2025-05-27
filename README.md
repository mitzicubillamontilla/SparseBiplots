# SparseBiplots

[![CRAN status](https://www.r-pkg.org/badges/version/SparseBiplots)](https://CRAN.R-project.org/package=SparseBiplots)

## Overview

Sparse Biplots are modern variant of HJ-Biplot which attempts to find sparse loadings, i.e., weight vectors with only a few nonzero values. Introducing a penalty parameter helps improve the interpretation of the model because the axis are formed as a linear combination of only some of the original variables.

This package provides functions to perform the HJ-Biplot (Galindo, 1986) and modifications introducing Ridge, LASSO and Elastic Net penalty:

* `HJBiplot`: performs the HJ-Biplot .
* `Ridge_HJBiplot`: performs the HJ-Biplot introducing Ridge penalty.
* `LASSO_HJBiplot`: performs the HJ-Biplot introducing LASSO penalty.
* `ElasticNet_HJBiplot`: performs the HJ-Biplot introducing Elastic Net penalty. For this use the [spca](https://github.com/erichson/spca) package.
* `Plot_Biplot`: create the plot of the results obtained with any of the above functions. It generates elegant data visualization using [ggplot2](https://github.com/tidyverse/ggplot2) with little code

## Installation

* Install *SparseBiplots* from **CRAN**:
```{r Install from CRAN}
install.packages("SparseBiplots")
```

* Install *SparseBiplots* from **GitHub** using [devtools](https://github.com/r-lib/devtools):
```{r Install from GitHub}
install.packages("devtools")
devtools::install_github("mitzicubillamontilla/SparseBiplots")
```

## Usage

### Parameters

* `X`: data frame which provides the data to be analyzed. All the variables must be numeric.
* `Transform.Data`: character indicating whether the columns have to be scaled or centered. Allowed values are: `"center"` and `"scale"` (default).
* `lambda`: tuning parameter for each penalty. In Elastic Net it refers to LASSO penalty.
* `àlpha`: tuning parameter of the Ridge shrinkage (just allowed with Elastic Net penalty).
* `operator`: operator used to solve the norm L1 (just allowed with LASSO penalty). Allowed values are: `'Hard-Thresholding'` and `'Soft-Thresholding'`. 
* `ind_name`: logical value which indicates if prints the name for each row of X (FALSE by default).
* `vec_name`: logical value which indicates if prints the name for each column of X (TRUE by default).

### Values

* `loadings`: loadings matrix.
* `coord_ind`: row coordinates matrix (individuals).
* `coord_var`: column coordinates matrix (variables).
* `eigenvalues`: approximated eigenvalues vector.
* `explvar`: vector containing the proportion of variance explained by the k axis obtained.
* `n_ceros`: matrix which indicates the number of loadings equal to cero in each axis (output exclusive just for LASSO and Elastic net penalties)

### Plot_Biplot parameters
* `X`: list containing the output of one of the functions of the package. 
* `axis`: vector with lenght 2 which contains the axis ploted in x and y axis.
* `hide`: vector specifying the elements to be hidden on the plot. Default value is “none”. Other allowed values are “ind” and “var”.
* `labels`: it indicates the label for points. If it is "auto" the labels are the row names of the coordinates of individuals. If it isn't auto it would be a vector containing the labels.
* `ind.shape`: points shapes. It can be a number to indicate the shape of all the points or a factor to indicate different shapes.
* `ind.color`: points colors. It can be a character indicating the color of all the points or a factor to use different colors.
* `ind.size`: numeric value indicating the size of points.
* `ind.label`: logical name indicating if prints the row names. 
* `ind.label.size`: numeric value indicating the size of the label of points.
* `var.col`: character indicating the color of the arrows.
* `var.size`: size of arrow.
* `var.label`: logical name indicating if prints the column names. 
* `var.label.size`: numeric value indicating the size of the labels of variables.
* `var.label.angle`: logical value indicating if print the vector names with orentation of the angle of the vector. 


## Example

```{r Fit HJ-Biplot}
library("SparseBiplots")
hj.biplot <- HJBiplot(mtcars)
```
It fit the HJ-Biplot on the R data `mtcars` returning a list that contains the results of the fiting (`hj_mtcars`). To create the plot use the function `Plot_Biplot` over the list obtained. 

```{r Visualize HJ-Biplot}
Plot_Biplot(hj.biplot, ind.label = TRUE)
```
<img src="https://github.com/mitzicubillamontilla/SparseBiplots/blob/master/plots/HJBiplot_example.png" width="750">

## References

* [Cubilla-Montilla, M. I. (2019). Contribuciones al análisis biplot basadas en soluciones factoriales disjuntas y en soluciones sparse (Doctoral dissertation, Universidad de Salamanca).](https://gredos.usal.es/handle/10366/140389)
* [Erichson, N. B., et al. (2018). Sparse principal component analysis via variable projection.](https://arxiv.org/abs/1804.00341)
* [Galindo, M. P. (1986). Una Alternativa de Representación Simultanea: HJ-Biplot. Qüestiió, 10:13-23.](http://diarium.usal.es/pgalindo/files/2012/07/0article-HJ-1986.pdf)
