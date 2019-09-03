##########################################################################
# Phylogenetic analyses and visualizations for fish-plastic meta-analysis
##########################################################################

# load packages and data ----
library(tidyverse)
library(gbm)
library(dismo)
library(mgcv)
library(lme4)
library(gamm4)
library(readxl)
library(readr)
library(ggplot2)
library(rotl) #for phylogenetic analyses, get all the species? from Hinchliff et al. 2015 PNAS
library(phytools)

d1 = read_xlsx("Plastics ingestion records fish master_updated.xlsx") %>% 
  mutate(Found = case_when(Habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                           Habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")) %>% 
  filter(!is.na(`Prop w plastic`))

d_sp_summ <- d1 %>% 
  group_by(`Species name`) %>% 
  #filter(str_subset(`Species name`, "sp", negate = TRUE)) %>%  View
  #  filter(N > 500 & Found == "pelagic") %>% 
  summarize(Sp_mean = mean(`Prop w plastic`),
            Sample_size = sum(N),
            num_studies = n_distinct(Source)) %>% 
  arrange(-Sp_mean)

# a=d_sp_summ[-grep(".$", d_sp_summ$`Species name`), ]

d_sp_summ <- subset(d_sp_summ, !grepl(".$", d_sp_summ$`Species name`))

gets the species names
taxa <- d_sp_sum$`Species name`[1:20]  # this is just a subset, will eventually include all species
# 
# resolved_names <- tnrs_match_names(taxa)
# 
# #plots species
# my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
# plot(my_tree, no.margin=TRUE)
# 
# tree<-read.tree(my_tree) # not working, but doesnt seem to matter atm
# 
# plot.tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
# plot.tree <- ladderize(plot.tree, right = TRUE)
# 
# tree.angle <- 270
# tree.start <- 180
# treeheight <- max(nodeHeights(plot.tree))
# 
# # consider adding in colors or bars next
# 
# # This tree plots, but has issues
# plot(plot.tree, type = "fan", open.angle = 360 - tree.angle, rotate = 270, 
#      root.edge = TRUE, 
#      show.tip.label = TRUE, label.offset = .36, cex = 0.6, #tip.color = col1, add back in when we have color working
#      edge.width = 1.5, font = 3, 
#      x.lim = c(-1 * treeheight, 1.2 * treeheight), 
#      y.lim = c(-1 * treeheight, 1.2 * treeheight))
#  
# 
# 
# 
# ## Trying a plot with ggtree
# 
# # jumping through hoops to install ggtree
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # 
# # BiocManager::install("ggtree")
# 
# library(ggtree)
#     
# try <- ggtree(plot.tree) +
#   geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab(size = 2)
# try