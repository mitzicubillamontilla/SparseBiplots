#' @title HJ Biplot
#'
#' @description This function performs the representation of HJ Biplot (Galindo, 1986).
#'
#' @usage HJBiplot (X, transform_data = 'scale', ind_name=FALSE,
#'     vec_name = TRUE)
#'
#' @param X array_like; \cr
#'     A data frame which provides the data to be analyzed. All the variables must be numeric.
#'
#' @param transform_data character; \cr
#'     A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. For default is "scale".
#'
#' @param ind_name bool; \cr
#'     If it is TRUE it prints the name for each row of X. If it is FALSE (default) does not print the names.
#'
#' @param vec_name bool; \cr
#'     If it is TRUE (default) it prints the name for each column of X. If it FALSE does not print the names.
#'
#' @details Algorithm used to construct the HJ Biplot. The Biplot is obtained as result of the configuration of markers for individuals and markers for variables in a reference system defined by the factorial axes resulting from the Decomposition in Singular Values (DVS).
#'
#' @return \code{HJBiplot} returns a list containing the following components:
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
#' \item{eigenvalues}{  array_like; \cr
#'           vector with the eigenvalues.
#'           }
#'
#' \item{explvar}{  array_like; \cr
#'           an vector containing the proportion of variance explained by the first 1, 2,.,k principal components obtained.
#'           }
#'
#' @author Mitzi Cubilla-Montilla, Carlos Torres, Ana Belen Nieto Librero and Purificacion Galindo Villardon
#'
#' @references
#' \itemize{
#'  \item Gabriel, K. R. (1971). The Biplot graphic display of matrices with applications to principal components analysis. Biometrika, 58(3), 453-467.
#'  \item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.
#' }
#'
#' @examples
#'  data(mtcars)
#'  HJBiplot(mtcars, transform_data = 'scale', ind_name  = TRUE)
#'
#' @importFrom graphics abline arrows plot text
#'
#' @export

HJBiplot = function(X, transform_data = 'scale', ind_name = FALSE, vec_name = TRUE){

  # List of objects that the function returns
  hjb = list(loadings = NULL,
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
  V = svd$v

  # Components' names
  PCs = vector()
  for (i in 1:dim(V)[2]){
    npc = vector()
    npc = paste(c("PC",i), collapse = "")
    PCs = cbind(PCs, npc)
  }

  # Objects returned by the function
  hjb$loadings = V
  row.names(hjb$loadings) = vec_tag #update row
  colnames(hjb$loadings) = PCs #update col
  hjb$coord_ind = X%*%V
  colnames(hjb$coord_ind) = PCs #update col
  hjb$coord_var = t(D%*%t(V))
  row.names(hjb$coord_var) = vec_tag #update
  colnames(hjb$coord_var) = PCs #update col
  hjb$eigenvalues = t(as.matrix(d))
  colnames(hjb$eigenvalues) = PCs #update col
  hjb$explvar = t(as.matrix(round(d^2/sum(d^2), digits = 4)*100))
  colnames(hjb$explvar) = PCs #update col

  # Limits of the plot
  xmin = min(hjb$coord_ind[,1], hjb$coord_var[,1]) - 0.5
  xmax = max(hjb$coord_ind[,1], hjb$coord_var[,1]) + 0.5
  ymin = min(hjb$coord_ind[,2], hjb$coord_var[,2]) - 0.5
  ymax = max(hjb$coord_ind[,2], hjb$coord_var[,2]) + 0.5

  # x and y labels
  var1 = paste(c("(", hjb$explvar[1], "%", ")"), collapse = "")
  eje1 = paste(PCs[1], var1)
  var2 = paste(c("(", hjb$explvar[2], "%", ")"), collapse = "")
  eje2 = paste(PCs[2], var2)

  # Plot
  plot(hjb$coord_ind,
       col ="blue",
       xlab = eje1,
       ylab = eje2,
       pch = 19,
       xlim = c(xmin, xmax),
       ylim = c(ymin, ymax))
  arrows(0, 0, hjb$coord_var[,1], hjb$coord_var[,2],
         col="green4", length = 0.1, angle = 20,lwd=2)
  abline(h = 0, v = 0, col="gray30", lwd=1, lty = 2)

  # Print sample's tags (if required)
  if (ind_name == TRUE){
    i=1
    for (tag in ind_tag){
      text(hjb$coord_ind[i, 1]+0.2, hjb$coord_ind[i, 2]+0.2,
           labels = ind_tag[i], col ="blue", cex = 0.65, font=2)
      i=i+1
    }
  }

  # Print varables' tags (if required)
  if(vec_name == TRUE){
    i=1
    for(tag in vec_tag){
      x = hjb$coord_var[i, 1]
      y = hjb$coord_var[i, 2]
      if (x > 0){
        text(x, y, labels = vec_tag[i], col = "green4",
             cex = 0.75, pos = 4, font=2)
        } else{
          text(x, y, labels = vec_tag[i], col = "green4",
               cex = 0.75, pos = 2, font=2)
          }
      i=i+1
    }
  }

  hjb

}
