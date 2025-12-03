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
    alpha=0.6)+
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12))
}