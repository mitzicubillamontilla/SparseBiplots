#' @title Elastic Net HJ Biplot
#'
#' @description This function is a generalization of the Ridge regularization method and the LASSO penalty. Realizes the representation of the SPARSE HJ Biplot through a combination of LASSO and Ridge, on the data matrix. This means that with this function you can eliminate weak variables completely as with the LASSO regularization or contract them to zero as in Ridge.
#'
#' @usage ElasticNet_HJBiplot(X, Lambda = 1e-04, Alpha = 1e-04, Transform.Data = 'scale')
#'
#' @param X array_like; \cr
#'     A data frame with the information to be analyzed
#'
#' @param Transform.Data character; \cr
#'    A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. Default is "scale".
#'
#' @param Lambda  float; \cr
#'     Tuning parameter of the LASSO penalty. Higher values lead to sparser components.
#'
#' @param Alpha  float; \cr
#'     Tuning parameter of the Ridge shrinkage
#'
#' @details Algorithm used to perform automatic selection of variables and continuous contraction simultaneously. With this method, the model obtained is simpler and more interpretable. It is a particularly useful method when the number of variables is much greater than the number of observations.
#'
#' @return \code{ElasticNet_HJBiplot} returns a list containing the following components:
#' \item{loadings}{  array_like; \cr
#'           penalized loadings, the loadings of the sparse principal components.
#'           }
#'
#'\item{n_ceros}{  array_like; \cr
#'           number of loadings equal to cero in each component.
#'           }
#'
#' \item{coord_ind}{  array_like; \cr
#'           matrix with the coordinates of individuals.
#'           }
#'
#' \item{coord_var}{  array_like; \cr
#'           matrix with the coordinates of variables.
#'           }
#'
#' \item{eigenvalues}{  array_like; \cr
#'           vector with the eigenvalues penalized.
#'           }
#'
#' \item{explvar}{  array_like; \cr
#'           an vector containing the proportion of variance explained by the first 1, 2,.,k sparse principal components obtained.
#'           }
#'
#' @author Mitzi Cubilla-Montilla, Carlos Torres-Cubilla, Ana Belen Nieto Librero and Purificacion Galindo Villardon
#'
#' @references
#' \itemize{
#'  \item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.
#'  \item Erichson, N. B., Zheng, P., Manohar, K., Brunton, S. L., Kutz, J. N., & Aravkin, A. Y. (2018). Sparse principal component analysis via variable projection. arXiv preprint arXiv:1804.00341.
#'  \item Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.
#' }
#'
#' @seealso \code{\link{spca}}, \code{\link{Plot_Biplot}}
#'
#' @examples
#'  ElasticNet_HJBiplot(mtcars, Lambda = 0.2, Alpha = 0.1)
#'
#' @import sparsepca
#'
#' @export

ElasticNet_HJBiplot <- function(X, Lambda = 1e-04, Alpha = 1e-04, Transform.Data = 'scale') {

  # List of objects that the function returns
  hj_elasticnet <-
    list(
      eigenvalues = NULL,
      explvar = NULL,
      loadings = NULL,
      n_ceros = NULL,
      coord_ind = NULL,
      coord_var = NULL
      )

  # Sample's tags
  ind_tag <- rownames(X)

  # Variable's tags
  vec_tag <- colnames(X)


  #### 1. Transform data ####
  if (Transform.Data == 'center') {
    X <-
      scale(
        as.matrix(X),
        center = TRUE,
        scale = FALSE
      )
  }

  if (Transform.Data == 'scale') {
    X <-
      scale(
        as.matrix(X),
        center = TRUE,
        scale = TRUE
      )
  }

  if (Transform.Data == 'none') {
    X <- as.matrix(X)
  }


  #### 2. SVD decomposition ####
  svd <- svd(X)
  U <- svd$u
  d <- svd$d
  D <- diag(d)
  V <- svd$v


  #### 3. Components calculated ####
  PCs <- vector()
  for (i in 1:dim(V)[2]){
    npc <- vector()
    npc <- paste("Dim",i)
    PCs <- cbind(PCs, npc)
  }


  #### 4. Introduce sparsity ####
  V_L <- V
  spca <-
    spca(
      X,
      k=dim(X)[2],
      alpha = Alpha,
      beta = Lambda,
      center = FALSE,
      scale = FALSE,
      verbose = FALSE
      )
  V_L[, c(1:dim(spca$loadings)[2])] <- spca$loadings
  rownames(V_L) <- colnames(X)


  #### 5. Sparsity magnitude ####
  n_ceros <- vector()
  for (i in 1:ncol(V_L)){
    n_ceros <- cbind(n_ceros, sum((V_L[,i] == 0) * 1))
  }
  colnames(n_ceros) <- PCs[, 1:dim(n_ceros)[2]]


  ##### 6. Output ####

  #### >Eigenvalues ####
  hj_elasticnet$eigenvalues <- spca$eigenvalues
  names(hj_elasticnet$eigenvalues) <- PCs

  #### > Explained variance ####
  hj_elasticnet$explvar <-
    round(
    hj_elasticnet$eigenvalues / sum(hj_elasticnet$eigenvalues),
    digits = 4
    ) * 100

  #### >Loagings ####
  hj_elasticnet$loadings <- V_L
  row.names(hj_elasticnet$loadings) <- vec_tag #update row
  colnames(hj_elasticnet$loadings) <- PCs #update col

  #### > Sparsity magnitude ####
  hj_elasticnet$n_ceros <- n_ceros

  #### >Row coordinates ####
  hj_elasticnet$coord_ind <- X %*% hj_elasticnet$loadings

  #### >Column coordinates ####
  hj_elasticnet$coord_var <- t(D %*% t(hj_elasticnet$loadings))
  colnames(hj_elasticnet$coord_var) <- PCs[, 1:dim(hj_elasticnet$coord_var)[2]]


  hj_elasticnet

}
