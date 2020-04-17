#' @title Plotting Biplot
#'
#' @description \code{Plot_Biplot} initializes a ggplot2-based visualization of the caracteristics presented in the data analized by the Biplot selected.
#'
#' @usage Plot_Biplot(X, ind.name = FALSE, vec.name = TRUE, point.col = "red", arrow.col = "black",
#' axis = c(1,2), angle.vec = TRUE)
#'
#' @param X List containing the output of one of the functions of the package.
#'
#' @param ind.name Logical value, if it is TRUE it prints the name for each row of X. If it is FALSE (default) does not print the names.
#'
#' @param vec.name Logical value, if it is TRUE (default) it prints the name for each column of X. If it is FALSE does not print the names.
#'
#' @param point.col Character indicating the color of the points.
#'
#' @param arrow.col Character indicating the color of the arrows.
#'
#' @param axis Vector with lenght 2 which contains the axis ploted in x and y axis.
#'
#' @param angle.vec Logical value, if it it TRUE (default) it print the vector names with orentation of the angle of the vector. If it is FALSE the angle of all tags is 0.
#'
#' @return Return a \code{\link{ggplot2}} object.
#'
#' @author Mitzi Cubilla-Montilla, Carlos Torres-Cubilla, Ana Belen Nieto Librero and Purificacion Galindo Villardon
#'
#' @seealso \code{\link{HJBiplot}}, \code{\link{Ridge_HJBiplot}}, \code{\link{ElasticNet_HJBiplot}}
#'
#' @examples
#' hj.biplot <- HJBiplot(mtcars, Transform.Data = 'scale')
#' Plot_Biplot(hj.biplot)
#'
#' @import ggplot2
#'
#' @import ggrepel
#'
#' @export

Plot_Biplot <- function(X, groups = NULL, ind.name = FALSE, vec.name = TRUE,
                point.col = "red", arrow.col = "black",
                axis = c(1,2), angle.vec = TRUE){

  #### 1. Params ####

  #### >Groups ####
  if(is.null(groups)){
    groups = "1"
    hide.legend = TRUE
  } else {
    hide.legend = FALSE
  }

  #### >Axis ploted #####
  axis.x <- axis[1]
  axis.y <- axis[2]

  #### >Limits ####
  xmin <- min(X$coord_ind[,axis.x],
              X$coord_var[,axis.x]) - 0.5
  xmax <- max(X$coord_ind[,axis.x],
              X$coord_var[,axis.x]) + 0.5
  ymin <- min(X$coord_ind[,axis.y],
              X$coord_var[,axis.y]) - 0.5
  ymax <- max(X$coord_ind[,axis.y],
              X$coord_var[,axis.y]) + 0.5

  #### >Axis labels ####
  PCs <- names(X$eigenvalues)
  var1 <- paste(c("(", X$explvar[axis.x], "%", ")"), collapse = "")
  var2 <- paste(c("(", X$explvar[axis.y], "%", ")"), collapse = "")
  eje1 <- paste(PCs[axis.x], var1)
  eje2 <- paste(PCs[axis.y], var2)

  #### >Axis names ####
  # It is used to indicate the columns to plot
  x.var <- colnames(X$coord_ind)[axis.x]
  y.var <- colnames(X$coord_ind)[axis.y]

  #### >Angle names ####
  ifelse(
    vec.name == TRUE & angle.vec == TRUE,
    #Angulo para los nombres de las variables
    angle <- atan(X$coord_var[, y.var] / X$coord_var[, x.var]) * 360 / (2 * pi),
    angle <- rep(0, nrow(X$coord_var))
    )

  ##### 2. Plot ####

  biplot <-
    ggplot() +
    #### >Draw axis ####
    geom_hline(
      yintercept = 0
      ) +
    geom_vline(
      xintercept = 0
      ) +
    #### >Axis labels ####
    labs(
      x = eje1,
      y = eje2
      ) +
    #### >Plot points ####
    geom_point(
      aes(x = X$coord_ind[, x.var],
          y = X$coord_ind[, y.var],
          colour = groups)
      ) +
    #### >Plot arrows ####
      geom_segment(
        aes(x = 0,
            y = 0,
            xend = X$coord_var[, x.var],
            yend = X$coord_var[, y.var]
            ),
        arrow = arrow(length = unit(0.5, "cm")),
        colour = arrow.col
        )

  #### >Point names ####
  if (ind.name == TRUE){
    biplot <-
      biplot +
      geom_text_repel(
        aes(
          x = X$coord_ind[, x.var],
          y = X$coord_ind[, y.var],
          label = rownames(X$coord_ind),
          colour = groups)
        )
  }

  #### >Color points ####
  biplot <-
    biplot +
    scale_colour_manual(values = point.col)

  #### >Vector names ####
  biplot <-
    biplot +
    geom_text_repel(
      aes(
        x = X$coord_var[, x.var],
        y = X$coord_var[, y.var],
        label = rownames(X$coord_var),
        angle = angle,
        hjust = ifelse(X$coord_var[, x.var] > 0, 1, 0),
        vjust = ifelse(X$coord_var[, y.var] > 0, 1, 0)
        ),
      col = arrow.col,
      box.padding = 0.1,
      point.padding = 0.1
      )

  #### Hide legend ####
  if(hide.legend){
    biplot <-
      biplot +
      theme(legend.position = "none")
  }

  biplot

  }
