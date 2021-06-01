# ==========================================================================
# Author:   Aaron
# Aim: some function of my own plots
# Version:  1.0
# Date:     Nov 4, 2019
# ==========================================================================

require(ggplot2)
require(grid)
require(RColorBrewer)

#' GO hist plot
#'
#' @details GO hist plot
#'
#' @param data data.frame of enrichment result from ClusterProfiler enrichXxx family function
#' @param n number of terms to plot, default is 10
#' @param color color of the hist plot
#'
#' @return a ggplot object

myGOhist <- function(data, n = 10, color = NULL) {
    data$Name <- paste0(data$Description, "(", data$Count, ")")
    data$minus.log10p <- -log10(data$p.adjust)

    # check color
    if (is.null(color)) {
        color <- brewer.pal(11, "Spectral")[2]
    }

    # check rows
    if (nrow(data) < n) n <- nrow(data)
    data <- data[1:n, ]
    data$Name <- factor(data$Name, levels = rev(data$Name))

    ggplot(data, aes(x=Name, y = minus.log10p)) +
        geom_bar(width = 0.5, stat = "identity", fill = color, colour = color) +
        coord_flip() +
        labs(x="",y="-log10 FDR")+
        scale_y_continuous(expand = c(0, 0)) +   # remove the gap between geom and axis
        theme(axis.title=element_text(size=8,face="bold"),
              axis.text=element_text(size=8),
              axis.text.x.top = element_text(face = "bold"),
              legend.key.height = unit(8, "points"),
              legend.key.width = unit(10, "points"),
              legend.position = "top",
              legend.title=element_text(size=8,face="bold"),
              legend.text=element_text(face="bold"),
              panel.grid.major  = element_blank(),     # remove the grid
              panel.background = element_blank(),      # remove the grey backgroud
              axis.line=element_line(size = .3, colour="black"),
              strip.text.x=element_text(size=8,face="bold"))
}
