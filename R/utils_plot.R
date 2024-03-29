# Default ggplot theme for all plots of this package
my_ggplot_theme = function(cex = 1) 
{
  theme_light(base_size = 10 * cex) +
    theme(
      legend.position = "bottom",
      legend.key.size = grid::unit(.3 * cex, "cm"),
      panel.background = ggplot2::element_rect(fill = 'white')
    )
}