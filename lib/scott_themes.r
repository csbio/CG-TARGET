#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# Some custom themes for ggplot2
library(grid)
library(gridExtra)

scott_theme_1 <- function() {
    theme_classic() +
    theme(axis.title = element_text(size = 16)) +
    theme(axis.title.x = element_text(vjust = -1.5)) +
    theme(axis.title.y = element_text(vjust = 0.2)) +
    theme(plot.title = element_text(size = 20)) +
    theme(axis.text = element_text(size = 14)) +
    theme(plot.margin = unit(c(1,1,1,1), 'cm'))
}

    
    
