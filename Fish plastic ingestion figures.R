########################################################################################
# Preliminary visualizations and analyses for marine fish-plastic ingestion meta-analysis
########################################################################################

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
library(MCMCglmm)

d = read.csv("Plastics ingestion records fish master_updated.csv") %>% 
  mutate(WeightedProp = Prop.w.plastic*N,
         Found = as_factor(case_when(Habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                                  Habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")))

d1 = read_xlsx("Plastics ingestion records fish master_updated.xlsx") %>% 
  mutate(Found = case_when(Habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                           Habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")) %>% 
  filter(!is.na(`Prop w plastic`))

         
# summary tables ----
d_sp_summ <- d1 %>% group_by(`Species name`) %>%  
#  filter(N > 500 & Found == "pelagic") %>% 
  summarize(Sp_mean = mean(`Prop w plastic`),
            Sample_size = sum(N),
            num_studies = n_distinct(Source)) %>% 
  arrange(-Sp_mean)

Fisheries_summ <- d1 %>% 
  filter(`Prop w plastic` > 0) %>% 
  group_by(Commercial) %>% 
  summarize(Sp_mean = mean(NwP/N),
            Sample_size = sum(N),
            num_species = n_distinct(`Species name`))

# fish of concern for humans
concern_fish <- d1 %>% 
  group_by(`Species name`) %>% 
  filter(Commercial %in% c("commercial", "highly commercial") & `Prop w plastic` == 0)

# geographic summary of data
Fish_geo_summ <- d1 %>% 
  group_by(`Oceanographic province (from Longhurst 2007)`) %>% 
  summarize(num_studies = n_distinct(Source),
            num_sp = n_distinct(`Species name`),
            num_ind_studied = sum(N),
            prop_by_region = sum(NwP)/sum(N))


# preliminary plots ----
p <- ggplot(d1, 
           aes(`Trophic level via fishbase`, `Prop w plastic`, size= N, weight = N)) + 
  geom_point(alpha = 0.4) +  # Eventually add in foraging behavior here 
  geom_smooth(col = "blue", method = "loess") +
  xlab("Trophic level") +
  ylab("Proportion of individuals with ingested plastic") +
  annotate("text", x = c(2.75, 3.9), y= -0.05, 
           label = c("planktivorous", "piscivorous")) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold")) 
p + guides(size = FALSE)


p2 <- ggplot(d, aes(Habitat, Prop.w.plastic, size=N)) +
  #geom_jitter() + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2


# proportion ingestion plastic with depth (perhaps remove demersal here)
p3 <- ggplot(filter(d1, Habitat %in% c("demersal", "pelagic-neritic", "pelagic-oceanic", "reef-associated", "benthopelagic")),
            aes(-`Average depth`, `Prop w plastic`, size=N, weight = N)) +  #col=Family
  geom_point(alpha = 0.4, aes(color = Basin)) + 
  geom_smooth(col = "grey20", method = "loess") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Habitat, scales = "free_y", ncol = 1) +
  #ylim(0,1) + 
  #xlim(-500,0) +
  coord_flip() +
  theme_classic()
p3

p3_b <- ggplot(filter(d1, Found != "NA"), 
               aes(-`Average depth`, `Prop w plastic`, size=N, weight = N)) +
  geom_point(alpha = 0.4) + 
  geom_smooth(col = "blue", method = "loess") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Found, scales = "free_y", ncol = 1) +
  #ylim(0,1) + 
  xlim(-1500,0) +
  coord_flip() +
  xlab("Average depth") +
  ylab("Proportion of individuals with ingested plastic") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        strip.text = element_text(size = 10)) 
p3_b + guides(size = FALSE)



p4 <- ggplot(filter(d1, `Oceanographic province (from Longhurst 2007)` != "NA"), 
             aes(`Oceanographic province (from Longhurst 2007)`, `Prop w plastic`)) +
  geom_jitter() + 
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4


p5 <- ggplot(d, aes(Order, Prop.w.plastic)) +
  #geom_jitter() + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p5


p6 <- ggplot(d1, 
             aes(`Vulnerability score (via fishbase from Cheug et al 2005)`, `Prop w plastic`,
                 weight = N)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold"))
p6


p7 <- ggplot(filter(d1, `Prime_forage` != "NA"), 
                aes(`Prime_forage`, `Prop w plastic`)) +
  geom_jitter() + 
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p7


p8 <- ggplot(filter(d1, `Commercial` != "NA"), 
             aes(`Commercial`, `Prop w plastic`)) +
  geom_jitter() + 
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p8




<<<<<<< HEAD
# GLMM ----

# USING THIS ONE WITH SCALED DEPTH, certified by Gemma
glmm_FwP <- glmer(cbind(NwP, N-NwP) ~ scale(`Average depth`)*Found + `Trophic level via fishbase` + 
                    (1|Source) + (1|Order), data = d1, family = binomial)
summary(glmm_FwP)


# MCMC GLMM ----
MCMCglmm_FwP <- MCMCglmm(Prop.w.plastic ~ scale(Average.depth),
                         random = ~ Order,
                         data = d, family = "gaussian",
                         nitt = 1150, thin = 100, burnin = 150,
                         pr = TRUE, # To save the posterior mode of for each level or category in your random effect(s) we use pr = TRUE, which saves them in the $Sol part of the model output.
                         verbose = TRUE)



# trying a gamm ----
gamm_FwP <- gamm(Prop.w.plastic*N ~ s(Trophic.level.via.fishbase,k=5)+s(Average.depth, k=5)+Found, 
                 random=list(Order=~1), data=d)
### $gam to look at gam effects. $lme to look at random effects.
summary(gamm_FwP$gam)
plot(gamm_FwP$gam)


# this seems to be working
gamm_lmer_FwP <- gamm4(cbind(NwP, N-NwP) ~ s(Trophic.level.via.fishbase, k=5) + s(Average.depth, k=5) + Found, 
                        random = ~(1|Order) + (1|Source), data = d, family = binomial)
# summary(gamm_lmer_FwP$gam)
# plot(gamm_lmer_FwP$gam)


gamm_lmer_FwP <- gam(Prop.w.plastic ~ s(Trophic.level.via.fishbase, k=5) + s(Average.depth, k=5), 
                      data = subset(d, Found == "pelagic"), family = gaussian)


# playing with a BRT ----

## I think this is what I want, check with Steph
gbmFwP <- gbm.step(data=d, 
                   gbm.x = c(3,6,16,18,27,35), 
                   gbm.y = 9,   # this is NwP
                   #weights = 8,  # weighted by sample size
                   family = "gaussian", 
                   tree.complexity = 5,
                   learning.rate = 0.001, bag.fraction = 0.5)
summary(gbmFwP)
gbm.plot(gbmFwP)




 # testing phylogenetic analyses ----
library(rotl) #for phylogenetic analyses, get all the species? from Hinchliff et al. 2015 PNAS
library(phytools)


#gets the species names
taxa <- d_sp_sum$`Species name`[1:20]  # this is just a subset, will eventually include all species

resolved_names <- tnrs_match_names(taxa)

#plots species
my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
plot(my_tree, no.margin=TRUE)

tree<-read.tree(my_tree) # not working, but doesnt seem to matter atm

plot.tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
plot.tree <- ladderize(plot.tree, right = TRUE)

tree.angle <- 270
tree.start <- 180
treeheight <- max(nodeHeights(plot.tree))

# consider adding in colors or bars next

# This tree plots, but has issues
plot(plot.tree, type = "fan", open.angle = 360 - tree.angle, rotate = 270, 
     root.edge = TRUE, 
     show.tip.label = TRUE, label.offset = .36, cex = 0.6, #tip.color = col1, add back in when we have color working
     edge.width = 1.5, font = 3, 
     x.lim = c(-1 * treeheight, 1.2 * treeheight), 
     y.lim = c(-1 * treeheight, 1.2 * treeheight))
 



## Trying a plot with ggtree

# jumping through hoops to install ggtree
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
=======
# # GLMM ----
# glmm_FwP <- glmer(cbind(NwP, N-NwP) ~ `Trophic level via fishbase` + `Average depth`*Found + `Oceanographic province (from Longhurst 2007)` +
#                    (1|Source) + (1|Order), data = d1, family = binomial)
# summary(glmm_FwP)
>>>>>>> 017f2f8348b3b4f5fdf48db6b998d9cda0816ac3
# 
# 
# # trying a gamm ----
# gamm_FwP <- gamm(Prop.w.plastic*N ~ s(Trophic.level.via.fishbase,k=5)+s(Average.depth, k=5)+Found, 
#                  random=list(Order=~1), data=d)
# ### $gam to look at gam effects. $lme to look at random effects.
# summary(gamm_FwP$gam)
# plot(gamm_FwP$gam)
# 
# 
# # this seems to be working
# gamm_lmer_FwP <- gamm4(cbind(NwP, N-NwP) ~ s(Trophic.level.via.fishbase, k=5) + s(Average.depth, k=5) + Habitat, 
#                        random = ~(1|Order) + (1|Source), data = d, family = binomial)
# summary(gamm_lmer_FwP$gam)
# plot(gamm_lmer_FwP$gam)
# 
# 
# # playing with a BRT ----
# 
# ## I think this is what I want, check with Steph
# gbmFwP <- gbm.step(data=d, 
#                    gbm.x = c(3,6,16,18,27,34), 
#                    gbm.y = 7,   # this is NwP
#                    weights = 8,  # weighted by sample size
#                    family = "poisson", 
#                    tree.complexity = 5,
#                    learning.rate = 0.001, bag.fraction = 0.5)
# summary(gbmFwP)
# gbm.plot(gbmFwP)
# 
# 
# 
# 
# # testing phylogenetic analyses ----
# library(rotl) #for phylogenetic analyses, get all the species? from Hinchliff et al. 2015 PNAS
# library(phytools)
# 
# 
# #gets the species names
# taxa <- d_sp_sum$`Species name`[1:20]  # this is just a subset, will eventually include all species
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
