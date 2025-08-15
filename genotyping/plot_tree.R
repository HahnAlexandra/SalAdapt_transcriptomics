new_lib <- "~/new_R_library"
.libPaths(c(new_lib, .libPaths()))

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggtree)
    library(treeio)

})

grouping <- read_tsv("./labels.txt") %>% distinct()

tree <- read.nexus(file="./output/tonsa_mb.nex.con.tre")
x <- as_tibble(tree)

dm <- left_join(x, grouping, by="label")



p <- ggtree(tree) + theme_tree()


pout <- p %<+% grouping + 
    geom_tiplab(aes(color=group), size=0.9) +
    theme(legend.position="right")+ 
    #geom_text( show.legend  = F ) +
    geom_tippoint(aes(color=group), size=0.9) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = c("A_hudsonica" = "#1B9E77",
                               "A_lilljeborgi" = "#D95F02",
                               "F" = "#7570B3",
                               "IV" = "#E7298A",
                               "out_group" = "#666666",
                               "S" = "#66A61E",
                               "SB" = "#E6AB02",
                               "X" = "#A6761D"),
                    na.value = "black")

ggsave(pout, file="./output/tree_plot.pdf", h=13, w=10)
