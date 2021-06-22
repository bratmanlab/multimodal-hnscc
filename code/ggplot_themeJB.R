theme_JB <- function(base_size = 9, base_family = "sans", 
                     plot.background = plot.background, ...) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family) + 
      theme(plot.title = element_text(face = "plain",
                                      size = 9, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            panel.border = element_blank(),
            plot.background = element_rect(color = NA),
            axis.title = element_text(size = 9, face = "plain"),
            axis.title.y = element_text(angle=90, vjust = 2,
                                        margin = margin(0,5,0,0)),
            axis.title.x = element_text(vjust = 1,
                                        margin = margin(5,0,0,0)),
            axis.text = element_text(size = 9),
            axis.line = element_line(colour = "black", size = 0.13),
            axis.ticks = element_line(size = 0.13),
            axis.ticks.length = unit(0.15, "cm"),
            panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = margin(0,0,0,0),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face = "plain", size = rel(1)),
            plot.margin = unit(c(5,5,5,5), "mm"),
            strip.background = element_rect(colour = NA, fill = NA),
            strip.text = element_text(face = "plain", vjust = 0.5, hjust = 0.5),
            ...
      ))
}

scale_fill_JB <- function(...){
  library(scales)
  library(ggsci)
  discrete_scale("fill", "Publication", manual_pal(values = c("#ef3b2c","#386cb0","#fdb462","#7fc97f","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_JB <- function(...){
  library(scales)
  library(ggsci)
  discrete_scale("colour", "Publication", manual_pal(values = c("#ef3b2c","#386cb0","#fdb462","#7fc97f","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

palette_JB <- c("#ef3b2c","#386cb0","#fdb462","#7fc97f","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")

theme_JBwBorder <- function(base_size = 9, base_family = "sans", 
                     plot.background = plot.background, ...) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family) + 
      theme(plot.title = element_text(face = "plain",
                                      size = 9, hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(color = "black", size = 0.13),
            panel.border = element_blank(),
            plot.background = element_rect(color = NA),
            axis.title = element_text(size = 9, face = "plain"),
            axis.title.y = element_text(angle=90, vjust = 2,
                                        margin = margin(0,5,0,0)),
            axis.title.x = element_text(vjust = 1,
                                        margin = margin(5,0,0,0)),
            axis.text = element_text(size = 9),
            axis.line = element_line(colour = "black", size = 0.13),
            axis.ticks = element_line(size = 0.13),
            axis.ticks.length = unit(0.15, "cm"),
            panel.grid.major = element_line(colour = NA),
            panel.grid.minor = element_line(colour = NA),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = margin(0,0,0,0),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face = "plain", size = rel(1)),
            plot.margin = unit(c(5,5,5,5), "mm"),
            strip.background = element_rect(colour = NA, fill = NA),
            strip.text = element_text(face = "plain", vjust = 0.5, hjust = 0.5),
            ...
      ))
}
