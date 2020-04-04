#' @title LASSO HJ Biplot
#'
#' @description This function performs the representation of the SPARSE HJ Biplot applying the LASSO regularization, on the original data matrix, implementing the norm L1.
#'
#' @usage LASSO_HJBiplot(X, lambda, transform_data = 'scale',
#'     operator = 'Hard-Thresholding',
#'     ind_name=FALSE, vec_name = TRUE)
#'
#' @param X array_like; \cr
#'      A data frame which provides the data to be analyzed. All the variables must be numeric.
#'
#' @param lambda  float; \cr
#'     Tuning parameter for the LASSO penalty
#'
#' @param transform_data character; \cr
#'     A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" that removes the columns means and "scale" that removes the columns means and divide by its standard deviation. For default is "scale".
#'
#' @param operator character; \cr
#'     The operator used to solve the norm L1.
#'
#' @param ind_name bool; \cr
#'     If it is TRUE it prints the name for each row of X. If it is FALSE (default) does not print the names.
#'
#' @param vec_name bool; \cr
#'     If it is TRUE (default) it prints the name for each column of X. If it is FALSE does not print the names.
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
#' @author Mitzi Cubilla-Montilla, Carlos Torres, Ana Belen Nieto Librero and Purificacion Galindo Villardon
#'
#' @references
#' \itemize{
#'  \item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.
#'  \item Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society: Series B (Methodological), 58(1), 267-288.
#'  \item Tibshirani, R. (2011). Regression shrinkage and selection via the lasso: a retrospective. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(3), 273-282.
#' }
#'
#' @examples
#'  data(mtcars)
#'  LASSO_HJBiplot(mtcars, 0.2, transform_data = 'scale', operator = 'Hard-Thresholding',
#'  ind_name = TRUE)
#'
#' @importFrom graphics abline arrows plot text
#'
#' @export
LASSO_HJBiplot = function(X, lambda, transform_data = 'scale',
                          operator = 'Hard-Thresholding',
                          ind_name = FALSE, vec_name = TRUE) {

  #Los valores permitidos para el par?metro "transform_data" son: 'scale' y 'center'

  # List of objects that the function returns
  hj_lasso = list(loadings = NULL,
                  n_ceros = NULL,
                  coord_ind = NULL,
                  coord_var = NULL,
                  eigenvalues = NULL,
                  explvar = NULL)

  # Sample's tags
  ind_tag = rownames(X)

  # Variable's tags
  vec_tag = colnames(X)

  # Transfomr the data
  if (transform_data == 'scale') {
    X = scale(as.matrix(X), scale = TRUE)
  }
  if (transform_data == 'center') {
    X = scale(as.matrix(X), center = TRUE)
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

  # Penalize the weights by the penalization
  if (operator == 'Soft-Thresholding'){
    V_LASSO = abs(V)-lambda
    V_LASSO[V_LASSO<=0] = 0
    V_LASSO = sign(V) * V_LASSO
  }
  if (operator == 'Hard-Thresholding'){
    V_LASSO = V
    V_LASSO[abs(V_LASSO)<=lambda]=0
  }

  #Number of null weights by each component
  n_ceros = vector()
  for (i in 1:ncol(V_LASSO)){
    n_ceros = cbind(n_ceros, sum((V_LASSO[,i] == 0) * 1))
  }
  colnames(n_ceros) = PCs[, 1:dim(n_ceros)[2]]

  # Objects returned by the function
  hj_lasso$loadings = V_LASSO #CARGAS
  row.names(hj_lasso$loadings) = vec_tag #update row
  colnames(hj_lasso$loadings) = PCs #update col
  hj_lasso$n_ceros = n_ceros #N CEROS
  hj_lasso$coord_ind = X%*%hj_lasso$loadings #INDIVIDUOS
  hj_lasso$coord_var = t(D%*%t(hj_lasso$loadings)) #VARIABLES
  colnames(hj_lasso$coord_var) = PCs[, 1:dim(hj_lasso$coord_var)[2]] #update
  QR = qr(hj_lasso$coord_ind) #Descomposicion QR
  R = qr.R( QR )
  hj_lasso$eigenvalues = round(abs(diag(R)),digits=4) #AUTOVALORES
  vari=hj_lasso$eigenvalues^2
  hj_lasso$explvar = round(vari/sum(vari), digits = 4)*100 #VARIANZA EXPLICADA

  # Limits of the plot
  xmin = min(hj_lasso$coord_ind[,1], hj_lasso$coord_var[,1]) - 0.5
  xmax = max(hj_lasso$coord_ind[,1], hj_lasso$coord_var[,1]) + 0.5
  ymin = min(hj_lasso$coord_ind[,2], hj_lasso$coord_var[,2]) - 0.5
  ymax = max(hj_lasso$coord_ind[,2], hj_lasso$coord_var[,2]) + 0.5

  # x and y labels
  var1 = paste(c("(", hj_lasso$explvar[1], "%", ")"), collapse = "")
  axis1 = paste(PCs[1], var1)
  var2 = paste(c("(", hj_lasso$explvar[2], "%", ")"), collapse = "")
  axis2 = paste(PCs[2], var2)

  # Plot
  plot(hj_lasso$coord_ind, col = "navy",
       xlab = axis1, ylab = axis2, pch = 20,
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  arrows(0, 0, hj_lasso$coord_var[,1], hj_lasso$coord_var[,2],
         length = 0.1, angle = 20, col="red")
  abline(h = 0, v = 0, col="gray30", lwd=1, lty = 2)

  # Print sample's tags (if required)
  if (ind_name == TRUE){
    i=1
    for (tag in ind_tag){
      text(hj_lasso$coord_ind[i, 1]+0.2, hj_lasso$coord_ind[i, 2]+0.2,
           labels = ind_tag[i], col = "navy", cex = 0.75,font=2, offset = 0.01)
      i=i+1
    }
  }

  # Print variables' tags (if required)
  if(vec_name == TRUE){
    i=1
    for(tag in vec_tag){
      x = hj_lasso$coord_var[i, 1]
      y = hj_lasso$coord_var[i, 2]
      if (x > 0){
        text(x, y, labels = vec_tag[i], col = "red",
             cex = 0.75, pos = 4, offset = 0.1)
      } else{
        text(x, y, labels = vec_tag[i], col = "red",
             cex = 0.75, font=2, pos = 2, offset = 0.1)
      }
      i=i+1
    }
  }

  hj_lasso

}
