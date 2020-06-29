#' @title LASSO HJ Biplot 
#'
#' @description This function performs the representation of the SPARSE HJ Biplot applying the LASSO regularization, on the original data matrix, implementing the norm L1.
#'
#' @usage LASSO_HJBiplot(X, Lambda, Transform.Data = 'scale', Operator = 'Hard-Thresholding')
#'
#' @param X array_like; \cr
#'      A data frame which provides the data to be analyzed. All the variables must be numeric.
#'
#' @param Transform.Data character; \cr
#'     A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. Default is "scale".
#'
#' @param Lambda  float; \cr
#'     Tuning parameter for the LASSO penalty
#'
#' @param Operator character; \cr
#'     The operator used to solve the norm L1. Allowed values are "Soft-Thresholding" and "Hard-Thresholding".
#'
#' @details Algorithm that performs a procedure of contraction and selection of variables. LASSO imposes a penalty that causes the charges of some components to be reduced to zero. By producing zero loadings for some components and not zero for others, the Lasso technique performs selection of variables. As the value of the penalty approaches one, the loadings approach zero.
#'
#' @return \code{LASSO_HJBiplot} returns a list containing the following components:
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
#'  \item Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society: Series B (Methodological), 58(1), 267-288.
#'  \item Tibshirani, R. (2011). Regression shrinkage and selection via the lasso: a retrospective. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(3), 273-282.
#' }
#'
#' @seealso \code{\link{Plot_Biplot}}
#'
#' @examples
#'  LASSO_HJBiplot(mtcars, Lambda = 0.2, Operator = 'Hard-Thresholding')
#'
#' @import stats
#'
#' @export
LASSO_HJBiplot <- function(X, Lambda, Transform.Data = 'scale', Operator = 'Hard-Thresholding') {

  # List of objects that the function returns
  hj_lasso <-
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
  if (Operator == 'Soft-Thresholding'){
    V_LASSO <- abs(V) - Lambda
    V_LASSO[V_LASSO <= 0] <- 0
    V_LASSO <- sign(V) * V_LASSO
  }

  if (Operator == 'Hard-Thresholding'){
    V_LASSO <- V
    V_LASSO[abs(V_LASSO) <= Lambda] <- 0
  }

  #### 5. Sparsity magnitude ####
  n_ceros <- vector()
  for (i in 1:ncol(V_LASSO)){
    n_ceros = c(n_ceros, sum((V_LASSO[,i] == 0) * 1))
  }
  names(n_ceros) <- PCs[, 1:length(n_ceros)]

  ##### 6. Output ####

  #### >Loagings ####
  hj_lasso$loadings <- V_LASSO
  row.names(hj_lasso$loadings) <- vec_tag
  colnames(hj_lasso$loadings) <- PCs

  #### >Sparsity magnitude ####
  hj_lasso$n_ceros <- n_ceros

  #### >Row coordinates ####
  hj_lasso$coord_ind <- X %*% hj_lasso$loadings

  #### >Column coordinates ####
  hj_lasso$coord_var <- t(D %*% t(hj_lasso$loadings)) #
  colnames(hj_lasso$coord_var) <- PCs[, 1:dim(hj_lasso$coord_var)[2]]

  #### >Eigenvalues ####

  ## Calculate variance
  QR <- qr(hj_lasso$coord_ind) #Descomposicion QR
  R <- qr.R( QR )
  Variance <- abs(diag(R))
  Variance <- as.vector(Variance)

  ## Eigenvalues
  hj_lasso$eigenvalues <-
    eigen(
      cor(
        svd$u %*% diag(Variance) %*% t(V_LASSO)
        )
      )$values
  names(hj_lasso$eigenvalues) <- PCs

  #### > Explained variance ####
  hj_lasso$explvar <-
    round(
      hj_lasso$eigenvalues / sum(hj_lasso$eigenvalues),
      digits = 4
    ) * 100


  hj_lasso

}
