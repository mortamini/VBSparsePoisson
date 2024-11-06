plot_credible_interval <- function(
  gg_density,  # ggplot object that has geom_density
  bound_left,
  bound_right
) {
  build_object <- ggplot_build(gg_density)
  x_dens <- build_object$data[[1]]$x
  y_dens <- build_object$data[[1]]$y
  
  index_left <- min(which(x_dens >= bound_left))
  index_right <- max(which(x_dens <= bound_right))
  
  gg_density + geom_area(
    data=data.frame(
      x=x_dens[index_left:index_right],
      y=y_dens[index_left:index_right]), 
    aes(x=x,y=y),
    fill="grey",
    alpha=0.6)
}