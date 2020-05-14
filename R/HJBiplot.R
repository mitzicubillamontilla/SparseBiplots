#' @title HJ Biplot
#'
#' @description This function performs the representation of HJ Biplot (Galindo, 1986).
#'
#' @usage HJBiplot (X, Transform.Data = 'scale')
#'
#' @param X array_like; \cr
#'     A data frame which provides the data to be analyzed. All the variables must be numeric.
#'
#' @param Transform.Data character; \cr
#'     A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. Default is "scale".
#'
#'
#' @details Algorithm used to construct the HJ Biplot. The Biplot is obtained as result of the configuration of markers for individuals and markers for variables in a reference system defined by the factorial axes resulting from the Decomposition in Singular Values (DVS).
#'
#' @return \code{HJBiplot} returns a list containing the following components:
#' \item{eigenvalues}{  array_like; \cr
#'           vector with the eigenvalues.
#'           }
#'
#' \item{explvar}{  array_like; \cr
#'           an vector containing the proportion of variance explained by the first 1, 2,.,k principal components obtained.
#'           }
#'
#' \item{loadings}{  array_like; \cr
#'           the loadings of the principal components.
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
#'  \item Gabriel, K. R. (1971). The Biplot graphic display of matrices with applications to principal components analysis. Biometrika, 58(3), 453-467.
#'  \item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.
#' }
#'
#' @seealso \code{\link{Plot_Biplot}}
#'
#' @examples
#'  HJBiplot(mtcars)
#'
#' @import stats
#'
#' @export

HJBiplot <- function(X, Transform.Data = 'scale'){

  # List of objects that the function returns
  hjb <-
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
  V <- svd$v


  #### 3. Components calculated ####
  PCs <- vector()
  for (i in 1:dim(V)[2]){
    npc <- vector()
    npc <- paste("Dim", i)
    PCs <- cbind(PCs, npc)
  }


  ##### 4. Output ####

  #### >Eigenvalues ####
  hjb$eigenvalues <- eigen(cor(X))$values
  names(hjb$eigenvalues) <- PCs

  #### > Explained variance ####
  hjb$explvar <-
    round(
      hjb$eigenvalues / sum(hjb$eigenvalues),
      digits = 4
    ) * 100
  names(hjb$explvar) <- PCs

  #### >Loagings ####
  hjb$loadings <- V
  row.names(hjb$loadings) <- vec_tag
  colnames(hjb$loadings) <- PCs

  #### >Row coordinates ####
  hjb$coord_ind <- X %*% V
  colnames(hjb$coord_ind) <- PCs

  #### >Column coordinates ####
  hjb$coord_var <- t(D %*% t(V))
  row.names(hjb$coord_var) <- vec_tag
  colnames(hjb$coord_var) <- PCs


  hjb

  }
