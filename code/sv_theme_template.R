library(ggplot2)

# Define ggplot theme
theme_sv <- function(){
  theme_bw(base_size=11, base_family="Roboto Light") %+replace%
  theme(
    panel.grid.major =  element_line(
      colour = "grey50",
      size = 0.2,
      linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey97"),
    plot.margin = unit(c(0.2, 1, 0.2, 1), "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10, family = "Roboto Medium", face = "bold"),
    strip.background = element_rect(fill= NULL, colour = "white", linetype = NULL),
    strip.text = element_text(colour = 'grey50', size = 9, vjust = 0.5, family = "Roboto Medium")
  )
}

# Colourblind palette for plots etc.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")