#' @title Ridge HJ Biplot
#'
#'
#' @description This function performs the representation of the SPARSE HJ Biplot applying the Ridge regularization, on the original data matrix, implementing the norm L2.
#'
#' @usage Ridge_HJBiplot (X, lambda, transform_data = 'scale', ind_name=FALSE,
#'     vec_name = TRUE)
#'
#' @param X array_like; \cr
#'     A data frame which provides the data to be analyzed. All the variables must be numeric.
#'
#' @param lambda  float; \cr
#'     Tuning parameter for the Ridge penalty
#'
#' @param transform_data character; \cr
#'     A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. For default it is "scale".
#'
#' @param ind_name bool; \cr
#'     If it is TRUE it prints the name for each row of X. If it is FALSE (default) does not print the names.
#'
#' @param vec_name bool; \cr
#'     If it is TRUE (default) it prints the name for each column of X. If it FALSE does not print the names.
#'
#' @details Algorithm used to contract the loads of the main components towards zero, but without achieving the nullity of any. If the penalty parameter is less than or equal to 1e-4 the result is like Galindo's HJ Biplot (1986).
#'
#' @return \code{Ridge_HJBiplot} returns a list containing the following components:
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
#' \item{eigenvalues}{  array_like; \cr
#'           vector with the eigenvalues penalized.
#'           }
#'
#' \item{explvar}{  array_like; \cr
#'           an vector containing the proportion of variance explained by the first 1, 2,.,k sparse principal components obtained.
#'           }
#'
#' @author Mitzi Cubilla-Montilla, Carlos Torres, Ana Belen Nieto Librero and Purificacion Galindo Villardon
#'
#' @references
#' \itemize{
#'  \item Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1), 55-67.
#'  \item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.
#'  \item Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse principal component analysis. Journal of computational and graphical statistics, 15(2), 265-286.
#' }
#'
#' @examples
#'  data(mtcars)
#'  Ridge_HJBiplot(mtcars, 0.2, transform_data = 'scale', ind_name = TRUE)
#'
#' @importFrom graphics abline arrows plot text
#'
#' @export

Ridge_HJBiplot = function(X, lambda, transform_data = 'scale', ind_name = FALSE,
                          vec_name = TRUE){

  # List of objects that the function returns
  hj_ridge = list(loadings = NULL,
                  coord_ind = NULL,
                  coord_var = NULL,
                  eigenvalues = NULL,
                  explvar = NULL)

  # Sample's tags
  ind_tag = rownames(X)

  # Variable's tags
  vec_tag = colnames(X)

  # Transform the data
  if (transform_data == 'center') {
    X = scale(as.matrix(X), center = TRUE)
  }

  if (transform_data == 'scale') {
    X = scale(as.matrix(X), scale = TRUE)
  }

  # SVD decomposition
  svd = svd(X)
  U = svd$u
  d = svd$d
  D = diag(d)
  V = svd$v/(1+lambda)

  # Components' names
  PCs = vector()
  for (i in 1:dim(V)[2]){
    npc = vector()
    npc = paste(c("PC",i), collapse = "")
    PCs = cbind(PCs, npc)
  }

  # Objects returned by the function
  hj_ridge$loadings = V
  row.names(hj_ridge$loadings) = vec_tag #update row
  colnames(hj_ridge$loadings) = PCs #update col
  hj_ridge$coord_ind = X%*%V
  colnames(hj_ridge$coord_ind) = PCs #update col
  hj_ridge$coord_var = t(D%*%t(V))
  row.names(hj_ridge$coord_var) = vec_tag #update row
  colnames(hj_ridge$coord_var) = PCs #update col
  QR = qr(hj_ridge$coord_ind) #Descomposicion QR
  R = qr.R( QR )
  hj_ridge$eigenvalues = round(abs(diag(R)),digits=4) #AUTOVALORES
  vari=hj_ridge$eigenvalues^2
  hj_ridge$explvar = round(vari/sum(vari), digits = 4)*100 #VARIANZA EXPLICADA

  # Limits of the plot
  xmin = min(hj_ridge$coord_ind[,1], hj_ridge$coord_var[,1]) - 0.3
  xmax = max(hj_ridge$coord_ind[,1], hj_ridge$coord_var[,1]) + 0.3
  ymin = min(hj_ridge$coord_ind[,2], hj_ridge$coord_var[,2]) - 0.3
  ymax = max(hj_ridge$coord_ind[,2], hj_ridge$coord_var[,2]) + 0.3

  # x and y labels
  var1 = paste(c("(", hj_ridge$explvar[1], "%", ")"), collapse = "")
  axis1 = paste(PCs[1], var1)
  var2 = paste(c("(", hj_ridge$explvar[2], "%", ")"), collapse = "")
  axis2 = paste(PCs[2], var2)

  # Plot
  plot(hj_ridge$coord_ind, col = "green4",
       xlab = axis1, ylab = axis2, pch = 20,
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  arrows(0, 0, hj_ridge$coord_var[,1], hj_ridge$coord_var[,2],
         col="blue", length = 0.1, angle = 20, lwd=2)
  abline(h = 0, v = 0, col="gray30", lwd=1, lty = 2)

  # Print sample's tags (if required)
  if (ind_name == TRUE){
    i=1
    for (tag in ind_tag){
      text(hj_ridge$coord_ind[i, 1]+0.15, hj_ridge$coord_ind[i, 2]+0.15,
           labels = ind_tag[i], col = "green4", cex = 0.65,font=2)
      i=i+1
    }
  }

  # Print varables' tags (if required)
  if(vec_name == TRUE){
    i=1
    for(tag in vec_tag){
      x = hj_ridge$coord_var[i, 1]
      y = hj_ridge$coord_var[i, 2]
      if (x > 0){
        text(x, y, labels = vec_tag[i], col = "blue",
             cex = 0.75, pos = 4,font=2)
      } else{
        text(x, y, labels = vec_tag[i], col = "blue",
             cex = 0.75, pos = 2, font=2)
      }
      i=i+1
    }
  }

  hj_ridge

}
