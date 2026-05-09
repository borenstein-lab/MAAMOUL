## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MAAMOUL)
library(igraph)
library(tidyverse)

## -----------------------------------------------------------------------------
data(edges)
head(edges)

## ----fig.height = 2, fig.width = 4--------------------------------------------
g <- graph_from_data_frame(edges, directed = FALSE)

# Degree distribution
deg <- degree(g)
ggplot(data.frame(deg), aes(x = deg)) +
  geom_histogram(color = 'black', fill = 'grey50', bins = 20) +
  ggtitle("Node degree distribution") +
  xlab('Node degree') +
  ylab('Count') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))

# Connected components
comp <- components(g)
sprintf("Number of connected components: %i", comp$no)
sprintf("Largest component size: %i", max(comp$csize))
sprintf("Smallest component size: %i", min(comp$csize))

# Basic stats
sprintf("Number of EC nodes: %i", sum(grepl("^EC", V(g)$name)))
sprintf("Number of metabolite nodes: %i", sum(!grepl("^EC", V(g)$name)))
sprintf("Number of edges: %i", ecount(g))

## -----------------------------------------------------------------------------
data(ec_pvals)
ec_pvals %>% select(feature, pval) %>% head()

## -----------------------------------------------------------------------------
data(mtb_pvals)
mtb_pvals %>% select(feature, pval) %>% head()

## -----------------------------------------------------------------------------
all_net_nodes <- unique(c(edges[[1]], edges[[2]]))

print("Loaded network information and feature p-values.")

sprintf(
  "%i of %i observed EC features (features for which a p-value has been computed) are also in the network.",
  sum(ec_pvals$feature %in% all_net_nodes),
  nrow(ec_pvals)
)

sprintf(
  "%i of %i observed metabolite features (features for which a p-value has been computed) are also in the network.",
  sum(mtb_pvals$feature %in% all_net_nodes),
  nrow(mtb_pvals)
)

