library(ggplot2)
library(extrafont)
library(viridis)

# font_import()

# Define ggplot theme
theme_tn <- function(){
  theme_bw(base_size=11) +
  theme(
    text = element_text(size = 11),
    # panel.grid.major =  element_line(
    #   colour = "grey90",
    #   size = 0.2,
    #   linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.margin = unit(c(0.2, 1, 0.2, 1), "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill= "white", colour = "white", linetype = NULL),
    strip.text = element_text(colour = 'grey10', size = 11, vjust = 0.5, hjust = 0, face = "bold"),
    axis.line = element_line(colour = "black", size = .3),
    axis.text = element_text(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
  legend.position = "bottom"
  )
}