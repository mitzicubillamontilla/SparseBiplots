#' @title Elastic Net HJ Biplot
#'
#' @description This function is a generalization of the Ridge regularization method and the LASSO penalty. Realizes the representation of the SPARSE HJ Biplot through a combination of LASSO and Ridge, on the data matrix. This means that with this function you can eliminate weak variables completely as with the LASSO regularization or contract them to zero as in Ridge.
#'
#' @usage ElasticNet_HJBiplot(X, lambda = 1e-04, alpha = 1e-04, transform_data = 'scale',
#'     ind_name = FALSE, vec_name = TRUE)
#'
#' @param X array_like; \cr
#'     A data frame with the information to be analyzed
#'
#' @param lambda  float; \cr
#'     Tuning parameter of the LASSO penalty. Higher values lead to sparser components.
#'
#' @param alpha  float; \cr
#'     Tuning parameter of the Ridge shrinkage
#'
#' @param transform_data character; \cr
#'    A value indicating whether the columns of X (variables) should be centered or scaled. Options are: "center" or "scale". For default is "scale".
#'
#' @param ind_name bool; \cr
#'     Logical value, if it is TRUE it prints the name for each row of X. If it is FALSE (default) does not print the names.
#'
#' @param vec_name bool; \cr
#'     Logical value, if it is TRUE (default) it prints the name for each column of X. If it is FALSE does not print the names.
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
#' @author Mitzi Cubilla-Montilla, Carlos Torres, Ana Belen Nieto Librero and Purificacion Galindo Villardon
#'
#' @references
#' \itemize{
#'  \item Galindo, M. P. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio, 10(1), 13-23.
#'  \item Erichson, N. B., Zheng, P., Manohar, K., Brunton, S. L., Kutz, J. N., & Aravkin, A. Y. (2018). Sparse principal component analysis via variable projection. arXiv preprint arXiv:1804.00341.
#'  \item Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.
#' }
#'
#' @examples
#'  data(mtcars)
#'  ElasticNet_HJBiplot(mtcars, 0.2, 0.1, transform_data = 'scale', ind_name=TRUE)
#'
#' @import sparsepca
#'
#' @importFrom graphics abline arrows plot text
#'
#' @export

ElasticNet_HJBiplot = function(X, lambda = 1e-04, alpha = 1e-04, transform_data = 'scale',
                               ind_name = FALSE, vec_name = TRUE) {

  # List of objects that the function returns
  hj_elasticnet = list(loadings = NULL,
                       n_ceros = NULL,
                       eigenvalues = NULL,
                       explvar = NULL,
                       coord_ind = NULL,
                       coord_var = NULL)

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
  V_L=V
  spca=spca(X, k=dim(X)[2], alpha = alpha, beta = lambda, scale = FALSE,
            verbose = FALSE)
  V_L[,c(1:dim(spca$loadings)[2])] = spca$loadings
  rownames(V_L)=colnames(X)

  #Number of null weights by each component
  n_ceros = vector()
  for (i in 1:ncol(V_L)){
    n_ceros = cbind(n_ceros, sum((V_L[,i] == 0) * 1))
  }
  colnames(n_ceros) = PCs[, 1:dim(n_ceros)[2]]

  # Objects returned by the function
  hj_elasticnet$loadings = V_L #CARGAS
  row.names(hj_elasticnet$loadings) = vec_tag #update row
  colnames(hj_elasticnet$loadings) = PCs #update col
  hj_elasticnet$n_ceros = n_ceros #N CEROS
  hj_elasticnet$coord_ind = X%*%hj_elasticnet$loadings #INDIVIDUOS
  hj_elasticnet$coord_var = t(D%*%t(hj_elasticnet$loadings)) #VARIABLES
  colnames(hj_elasticnet$coord_var) = PCs[, 1:dim(hj_elasticnet$coord_var)[2]] #update
  QR = qr(hj_elasticnet$coord_ind) #Descomposicion QR
  R = qr.R( QR )
  hj_elasticnet$eigenvalues = round(abs(diag(R)),digits=4) #AUTOVALORES
  vari=hj_elasticnet$eigenvalues^2
  hj_elasticnet$explvar = round(vari/sum(vari), digits = 4)*100 #VARIANZA EXPLICADA

  # Limits of the plot
  xmin = min(hj_elasticnet$coord_ind[,1], hj_elasticnet$coord_var[,1]) - 0.5
  xmax = max(hj_elasticnet$coord_ind[,1], hj_elasticnet$coord_var[,1]) + 0.5
  ymin = min(hj_elasticnet$coord_ind[,2], hj_elasticnet$coord_var[,2]) - 0.5
  ymax = max(hj_elasticnet$coord_ind[,2], hj_elasticnet$coord_var[,2]) + 0.5

  # x and y labels
  var1 = paste(c("(", hj_elasticnet$explvar[1], "%", ")"), collapse = "")
  axis1 = paste(PCs[1], var1)
  var2 = paste(c("(", hj_elasticnet$explvar[2], "%", ")"), collapse = "")
  axis2 = paste(PCs[2], var2)

  # Plot
  plot(hj_elasticnet$coord_ind, col = "darkgreen",
       xlab = axis1, ylab = axis2, pch = 20,
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  arrows(0, 0, hj_elasticnet$coord_var[,1], hj_elasticnet$coord_var[,2], col="red",length = 0.1, angle = 20)
  abline(h = 0, v = 0, col="gray30", lwd=1, lty = 2)

  # Print sample's tags (if required)
  if (ind_name == TRUE){
    i=1
    for (tag in ind_tag){
      text(hj_elasticnet$coord_ind[i, 1]+0.1, hj_elasticnet$coord_ind[i, 2]+0.1,
           labels = ind_tag[i], col = "darkgreen", cex = 0.75, offset = 0.01, font=2)
      i=i+1
    }
  }

  # Print variables' tags (if required)
  if(vec_name == TRUE){
    i=1
    for(tag in vec_tag){
      x = hj_elasticnet$coord_var[i, 1]
      y = hj_elasticnet$coord_var[i, 2]
      if (x > 0){
        text(x, y, labels = vec_tag[i], col = "red",
             cex = 0.75, pos = 4, offset = 0.1, font=2)
      } else{
        text(x, y, labels = vec_tag[i], col = "red",
             cex = 0.75, pos = 2, offset = 0.1,font=2)
      }
      i=i+1
    }
  }

  hj_elasticnet

}
