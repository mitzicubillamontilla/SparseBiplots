#' @title Ridge HJ Biplot 
#'
#' @description This function performs the representation of the HJ Biplot applying the Ridge regularization, on the original data matrix, implementing the norm L2.
#'
#' @usage Ridge_HJBiplot (X, Lambda, Transform.Data = 'scale')
#'
#' @param X array_like; \cr
#'     A data frame which provides the data to be analyzed. All the variables must be numeric.
#'
#' @param Transform.Data character; \cr
#'     A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. Default is "scale".
#'
#' @param Lambda  float; \cr
#'     Tuning parameter for the Ridge penalty
#'
#' @details Algorithm used to contract the loads of the main components towards zero, but without achieving the nullity of any. If the penalty parameter is less than or equal to 1e-4 the result is like Galindo's HJ Biplot (1986).
#'
#' @return \code{Ridge_HJBiplot} returns a list containing the following components:
#' \item{eigenvalues}{  array_like; \cr
#'           vector with the eigenvalues penalized.
#'           }
#'
#' \item{explvar}{  array_like; \cr
#'           an vector containing the proportion of variance explained by the first 1, 2,.,k sparse principal components obtained.
#'           }
#'
#' \item{loadings}{  array_like; \cr
#'           penalized loadings, the loadings of the sparse principal components.
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
#' @author Mitzi Cubilla-Montilla, Carlos Torres-Cubilla, Ana Belen Nieto Librero and Purificacion Galindo Villardon
#'
#' @references
#' \itemize{
#'  \item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.
#'  \item Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1), 55-67.
#'  \item Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse principal component analysis. Journal of computational and graphical statistics, 15(2), 265-286.
#' }
#'
#' @seealso \code{\link{Plot_Biplot}}
#'
#' @examples
#'  Ridge_HJBiplot(mtcars, Lambda = 0.2)
#'
#' @import stats
#'
#' @export

Ridge_HJBiplot <- function(X, Lambda, Transform.Data = 'scale'){

  # List of objects that the function returns
  hj_ridge <-
    list(
      eigenvalues = NULL,
      explvar = NULL,
      loadings = NULL,
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
  V <- svd$v/(1+Lambda)

  #### 3. Components calculated ####
  PCs <- vector()
  for (i in 1:dim(V)[2]){
    npc = vector()
    npc = paste("Dim",i)
    PCs = cbind(PCs, npc)
  }

  ##### 4. Output ####

  #### >Eigenvalues ####
  hj_ridge$eigenvalues <- eigen(cor(X))$values
  names(hj_ridge$eigenvalues) <- PCs

  #### > Explained variance ####
  hj_ridge$explvar <-
    round(
      hj_ridge$eigenvalues / sum(hj_ridge$eigenvalues),
      digits = 4
      ) * 100
  names(hj_ridge$explvar) <- PCs

  #### >Loagings ####
  hj_ridge$loadings <- V
  row.names(hj_ridge$loadings) <- vec_tag
  colnames(hj_ridge$loadings) <- PCs


  #### >Row coordinates ####
  hj_ridge$coord_ind = X%*%V
  colnames(hj_ridge$coord_ind) = PCs

  #### >Column coordinates ####
  hj_ridge$coord_var = t(D%*%t(V))
  row.names(hj_ridge$coord_var) = vec_tag
  colnames(hj_ridge$coord_var) = PCs


  hj_ridge

}
