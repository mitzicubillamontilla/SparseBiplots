#' @title Plotting Biplot
#'
#' @description \code{Plot_Biplot} initializes a ggplot2-based visualization of the caracteristics presented in the data analized by the Biplot selected.
#'
#' @usage Plot_Biplot(X, axis = c(1,2),
#'   color = "red", shape = 20, size = 4,
#'   ind.label = FALSE, ind.label.size = 4,
#'   arrow.col = "black", vec.name = TRUE, angle.vec = FALSE)
#'
#' @param X List containing the output of one of the functions of the package.
#'
#' @param axis Vector with lenght 2 which contains the axis ploted in x and y axis.
#'
#' @param color Points colors. It can be a character indicating the color of all the points or a factor to use different colors.
#'
#' @param shape Points shape. It can be a number to indicate the shape of all the points or a factor to indicate different shapes.
#'
#' @param size numeric value indicating the size of points.
#'
#' @param ind.label Logical value, if it is TRUE it prints the name for each row of X. If it is FALSE (default) does not print the names.
#'
#' @param ind.label.size numeric value indicating the size of the labels of points.
#'
#' @param arrow.col Character indicating the color of the arrows.
#'
#' @param vec.name Logical value, if it is TRUE (default) it prints the name for each column of X. If it is FALSE does not print the names.
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
#' Plot_Biplot(hj.biplot, ind.label = TRUE)
#'
#' @import ggplot2
#'
#' @import ggrepel
#'
#' @export

Plot_Biplot <- function(X, axis = c(1,2),
                        color = "red", shape = 19, size = 2,
                        ind.label = FALSE, ind.label.size = 4,
                        arrow.col = "black", vec.name = TRUE, angle.vec = FALSE
                        ){

  #### 1. Params ####

  #### >>Axis ploted #####
  axis.x <- axis[1]
  axis.y <- axis[2]

  #### >>Axis labels ####
  PCs <- names(X$eigenvalues)
  var1 <- paste(c("(", X$explvar[axis.x], "%", ")"), collapse = "")
  var2 <- paste(c("(", X$explvar[axis.y], "%", ")"), collapse = "")
  eje1 <- paste(PCs[axis.x], var1)
  eje2 <- paste(PCs[axis.y], var2)

  #### >>Axis names ####
  # It is used to indicate the columns to plot
  x.var <- colnames(X$coord_ind)[axis.x]
  y.var <- colnames(X$coord_ind)[axis.y]


  #### >>Subsect variables ####
  selection <- (X$coord_var[, x.var] == 0) & (X$coord_var[, y.var] == 0)
  new_coord_var <- subset(X$coord_var, subset = !selection)



  ##### 2. Plot ####

  #### >>Empty plot ####
  biplot <-
    ggplot() +
    #### >>Draw axis ####
    geom_hline(
      yintercept = 0
      ) +
    geom_vline(
      xintercept = 0
      ) +
    #### >>Axis labels ####
    labs(
      x = eje1,
      y = eje2
      )



  #### 3. Plot points ####

  #### >>Shape params ####
  if(length(shape) == 1){
    shape.aes <- factor(1)
  } else {
    shape.aes <- shape
    hide.point.shape <- "legend"}

  #### >>Add points ####
  biplot <-
    biplot +
    geom_point(
      aes(x = X$coord_ind[, x.var],
          y = X$coord_ind[, y.var],
          colour = color
          ,shape = shape.aes
          ),
      size = size
      )

  #### >>Point names ####
  if (ind.label == TRUE){
    biplot <-
      biplot +
      geom_text_repel(
        aes(
          x = X$coord_ind[, x.var],
          y = X$coord_ind[, y.var],
          label = rownames(X$coord_ind),
          colour = color),
        size = ind.label.size
      )
  }

  #### >>Colors ####
  if(length(color) == 1){
    hide.point.color <- FALSE
    biplot <-
      biplot +
      scale_colour_manual(values = color)
  } else {
    if(is.factor(color)) {
      hide.point.color <- "legend"
    } else {
      hide.point.color <- "colorbar"
    }
  }

  #### >>Change shape ####
  if(length(shape) == 1){
    hide.point.shape <- FALSE
    biplot <-
      biplot +
      scale_shape_manual(values = shape)
  }



  ##### 4. Plot arrows #####

  #### >>Add arrows ####
    biplot <- biplot +
    geom_segment(
      aes(x = 0,
          y = 0,
          xend = new_coord_var[, x.var],
          yend = new_coord_var[, y.var]
          ),
      arrow = arrow(length = unit(0.5, "cm")),
      colour = arrow.col
      )

  #### >>Angle names ####
  ifelse(
    vec.name == TRUE & angle.vec == TRUE,
    #Angulo para los nombres de las variables
    angle <- atan(new_coord_var[, y.var] / new_coord_var[, x.var]) * 360 / (2 * pi),
    angle <- rep(0, nrow(new_coord_var))
  )

  #### >>Vector names ####
  if(vec.name == TRUE){
    biplot <-
    biplot +
    geom_text_repel(
      aes(
        x = new_coord_var[, x.var],
        y = new_coord_var[, y.var],
        label = rownames(new_coord_var),
        angle = angle,
        hjust = ifelse(new_coord_var[, x.var] > 0, 1, 0),
        vjust = ifelse(new_coord_var[, y.var] > 0, 1, 0)
        ),
      col = arrow.col,
      box.padding = 0.1,
      point.padding = 0.1
      )
  }



  #### 5. Legend ####
  biplot <- biplot +
    guides(
      colour = hide.point.color,
      shape = hide.point.shape
      )

  biplot

  }
