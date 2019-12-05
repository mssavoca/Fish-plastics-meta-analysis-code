######################################################## 
# Figures for fish-plastic ingestion meta-analysis paper
########################################################

#load packages, data, and functions ----

library(tidyverse)
library(ggtree)
library(ggstance)
library(ape)
library(phylobase)
library(dismo)
library(readxl)
library(rotl) #for phylogenetic analyses, get all the species? from Hinchliff et al. 2015 PNAS
library(phytools)
library(tidyverse)
library(stringr)
library(mgcv)
library(gamm4)
library(MuMIn)
library(ggrepel)

SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}


# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

d_poll <- as_tibble(read_csv("Spatial Information_microplastics.csv"))

d = read_csv("Plastics ingestion records fish master_final_2019data.csv") %>% 
  janitor::clean_names() %>% 
  rename(NwP = nw_p,
         N = n) %>% 
  mutate(WeightedProp = prop_w_plastic * N,
         Found = as_factor(case_when(habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                                     habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")),
         prime_forage = na_if(prime_forage, "not listed")) %>% 
  rename(ProvCode = "oceanographic_province_from_longhurst_2007") %>% 
  separate(binomial, into = c("genus", "species"), sep = " ", remove = FALSE) %>% 
  left_join(dplyr::select(d_poll, c(ProvCode, adjacency, mean_poll_abund)), by = "ProvCode")

d_full <- d %>%
  filter(includes_microplastic == "Y")


# species summary tables
d_sp_sum <- d_full %>%
  filter(!species %in% c("sp.", "spp.")) %>%
  group_by(binomial, order, commercial) %>%
  summarize(Sp_mean = mean(prop_w_plastic),
            Sample_size = sum(N),
            num_studies = n_distinct(source)) %>% 
  ungroup %>% 
  mutate(commercial = factor(commercial),
         studies_cat = as.double(cut(num_studies, 
                                     c(0, 1, 3, 6),
                                     c(1,2,3))),
         commercial = fct_collapse(commercial,
                                   Commercial = c("commercial", "highly commercial"),
                                   Minor = c("minor commercial", "subsistence"),
                                   None = "none"))%>%
  arrange(-Sp_mean)




# Phylogenetic tree figure for paper, Figure 1 ----

# building the basic tree 

breaks <- c(seq(1,nrow(d_sp_sum),50),nrow(d_sp_sum)+1)  # why are we doing this?

for (i in 1:(length(breaks)-1)){
  taxa <- as.character(d_sp_sum$binomial[breaks[i]:(breaks[i+1]-1)])
  taxa <- taxa[taxa != "" & !is.na(taxa)]
  
  resolved_namest <- tnrs_match_names(taxa)                          # I think this is where all the extra species fall out
  resolved_namest <- resolved_namest[!is.na(resolved_namest$unique_name),]
  if (i==1){
    resolved_namess <- resolved_namest
  } else {
    resolved_namess <- rbind(resolved_namess, resolved_namest)
  }
}
resolved_names <- resolved_namess
resolved_names <- resolved_names[resolved_names$flags!="INCERTAE_SEDIS_INHERITED",]

#plots species
my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id, label_format = "name")

my_tree$tip.label<-gsub("_"," ",my_tree$tip.label) # removes underscore between genus and species names
my_tree$tip.label<-str_extract(my_tree$tip.label, "[A-Z][a-z]+ [a-z]+")

my_tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
my_tree <- ladderize(my_tree, right = TRUE)

#View(my_tree)

#my_tree2 = phylo4d(my_tree, d_sp_sum)


# first plot try, fan layout
p <- ggtree(my_tree, layout="fan", open.angle=0) + 
  #geom_text2(aes(label=label), hjust=-.2, size=4) +
  ggplot2::xlim(-0.6, 1.3) 
p

my_tree$tip.label <- as.factor(my_tree$tip.label)

# Adding data
p %<+% d_sp_sum + 
  aes(color = Sp_mean) +
  geom_tiplab2(aes(label = paste0("italic('", label, "')"),
                   color=Sp_mean, angle = angle), parse = TRUE,
               size = 6, align = FALSE, hjust = -0.05) +
  geom_tippoint(aes(color = Sp_mean, shape = commercial, size = studies_cat)) +
  scale_color_gradientn(colours = c("steelblue4",
                                    "gray40",
                                    "coral", "coral1",
                                    "firebrick2", "firebrick3", "firebrick4"), 
                        name = "Proportion with \ningested plastic") +
  #scale_size(range = c(3, 7)) +
  scale_size_continuous(guide = FALSE, range = c(3, 7)) +
  labs(shape = "Commercial \nstatus") +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(2.5, "cm"),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 26),
        legend.box = "horizontal") +
  guides(shape = guide_legend(override.aes = list(size = 5)))


# save plots
dev.copy2pdf(file="Prelim_phylo_ggtree2.pdf", width=35, height=35)
ggsave("Prelim_phylo_ggtree2.jpg", width = 35, height = 35, units = "in")


# risk/interest plot, either figure 1B or supplementary ----
risk_plot <- d_sp_sum %>% 
  drop_na(commercial) %>% 
  ggplot(aes(log10(Sample_size), Sp_mean)) +
  geom_point(aes(color = Sp_mean, shape = commercial, size = studies_cat), 
             alpha = 0.8) +
  geom_hline(yintercept = 0.25, linetype="dashed", color = "grey50") +
  geom_vline(xintercept = log10(25), linetype="dashed", color = "grey50") +
  xlab("Log[Sample Size]") +
  ylab("Species-specific plastic ingestion incidence (FO)") +
  labs(shape = "Commercial status", shape = "Commercial Status", size =  "Number of studies") +
  scale_color_gradientn(colours = c("steelblue4",
                                    "gray40",
                                    "coral", "coral1",
                                    "firebrick2", "firebrick3", "firebrick4"), 
                        name = "Proportion with \ningested plastic") +
  scale_size_continuous(breaks = seq(from = 1, to = 3, by = 1), 
                        labels = c("Poorly studied (n=1)", "Moderately studied (n=2-3)", "Well studied (n=4-6)"),
                        range = c(1.5, 5)) +
  annotate("text", x = c(0.6, 2.8, 0.6, 2.8),
           y=c(0.9, 0.9, 0.08, 0.08),
           label = c("high incidence, data poor", "high incidence, data rich",
                     "low incidence, data poor", "low incidence, data rich")) +
  theme_classic(base_size = 16)
risk_plot

dev.copy2pdf(file="risk_plot.pdf", width=12, height=7)


# Figure 3, plastic ingestion by depth and habitat ----

Fig_3 <- ggplot(filter(d_full, Found != "NA"), 
               aes(average_depth, prop_w_plastic, size = N, weight = N)) +
  geom_point(alpha = 0.4) + 
  geom_smooth(col = "blue", method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #facet_wrap(~ Found, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_x_reverse() +
  xlim(1500,0) +
  xlab("Average depth") +
  ylab("Proportion of individuals with ingested plastic") +
  theme_classic(base_size = 16)
  # theme(axis.text.x = element_text(size=12),
  #       axis.text.y = element_text(size=12),
  #       axis.title=element_text(size=13, face="bold"),
  #       strip.text = element_text(size = 12)) 
Fig_3 + guides(size = FALSE)

dev.copy2pdf(file="Fig_3.pdf", width=7, height=10)

# Figure 4, species accumulation curves----
## Cumulative unique entries in list of vectors
cum_unique <- function(l) {
  result <- integer(length(l))
  for (i in 1:length(l)) {
    result[i] <- length(unique(unlist(l[1:i])))
  }
  result
}
d_rarefaction_all <-  d_full %>% 
  group_by(publication_year) %>% 
  summarize(species = list(unique(binomial)),
            annual_N = sum(N)) %>% 
  ungroup %>% 
  mutate(cum_species = cum_unique(species),
         cum_n = cumsum(annual_N),
         cpue = cum_species / cum_n)
d_rarefaction_ingest_only <-  d_full %>%  
  filter(prop_w_plastic > 0) %>% 
  group_by(publication_year) %>% 
  summarize(species = list(unique(binomial)),
            annual_N = sum(N)) %>% 
  ungroup %>% 
  mutate(cum_species = cum_unique(species),
         cum_n = cumsum(annual_N),
         cpue = cum_species / cum_n)

rarefaction_plot <- ggplot() +
  geom_line(data = d_rarefaction_all, aes(cum_n, cum_species), color = "dark blue") +
  geom_line(data = d_rarefaction_ingest_only, aes(cum_n, cum_species), color = "dark red") +
  geom_point(data = d_rarefaction_all, aes(cum_n, cum_species)) +
  geom_point(data = d_rarefaction_ingest_only, aes(cum_n, cum_species)) +
  geom_text_repel(data = d_rarefaction_all, 
                  aes(cum_n, cum_species, label = publication_year), 
                  nudge_x = -1500, nudge_y = 10,
                  segment.color = "black") +
  geom_text_repel(data = d_rarefaction_ingest_only, 
                  aes(cum_n, cum_species, label = publication_year), 
                  nudge_x = 1500, nudge_y = -10,
                  segment.color = "black") +
  labs(x = "Cumulative number of individuals sampled",
       y = "Cumuluative number of species sampled") + 
  theme_classic(base_size = 14)
rarefaction_plot






# Figure S1, number of studies over time----
study_hist <- d %>% 
  group_by(publication_year) %>% 
  summarize(n_studies = n_distinct(source)) %>% 
  ggplot(aes(publication_year, n_studies)) + 
  geom_bar(stat = "identity") + 
  geom_smooth(se = FALSE) +
  theme_classic(base_size = 14) +
  xlab("Publication year") + 
  ylab("Number of studies") 
study_hist 


######################################################### 
# Modeling for fish-plastic ingestion meta-analysis paper
#########################################################

# trying out some GAMMs in gamm4, feel free to try with mgcv too

# wrapper fucntion  needed to allow MuMIn to compare models
gamm4 <- function(...) structure(c(gamm4::gamm4(...), list(call = match.call())), class = c("gamm", "list"))  

d_full_wo_gaps <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(average_depth, Found, trophic_level_via_fishbase, prime_forage, NwP, N, 
          adjacency, mean_poll_abund) 



# model testing both ecological and geographic variables
gamm4_FwP_eco_geo <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase, bs="ts") + 
                             s(scale(average_depth), bs="ts") + 
                             s(scale(mean_poll_abund), bs="ts") +
                             Found + prime_forage + adjacency, # mean poll abund (try spline and no spline) instead of this, and adjacency (no spline)
                   random = ~(1|order) + (1|source), 
                   REML = TRUE, # new peice from Elliot
                   data = d_full_wo_gaps, family = binomial)
summary(gamm4_FwP_eco_geo$gam)
summary(gamm4_FwP_eco_geo$dev.expl)
plot(gamm4_FwP_eco_geo$gam)



# model testing ecological variables
gamm4_FwP_eco <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase, bs="ts") + 
                         s(scale(average_depth), bs="ts") + 
                         Found + prime_forage, 
                       random = ~(1|order) + (1|source), 
                       REML = TRUE, # new peice from Elliot
                       data = d_full_wo_gaps, family = binomial)
summary(gamm4_FwP_eco$gam)
summary(gamm4_FwP_eco$dev.expl)
plot(gamm4_FwP_eco$gam)




# model testing geographic variable
gamm4_geo <- gamm4(cbind(NwP, N-NwP) ~ s(scale(mean_poll_abund), bs="ts") +   # oceanographic_province_from_longhurst_2007,  # swap longhurst region with mean poll abund (try spline and no spline) instead of this, and adjacency (no spline)
                   random = ~(1|order) + (1|source), 
                   REML = TRUE,
                   data = d_full_wo_gaps, family = binomial)
summary(gamm4_geo$gam)
plot(gamm4_geo$gam)


gamm4_geo <- gamm4(cbind(NwP, N-NwP) ~ s(scale(mean_poll_abund), bs="ts") + 
                             adjacency, # mean poll abund (try spline and no spline) instead of this, and adjacency (no spline)
                           random = ~(1|order) + (1|source), 
                           REML = TRUE, # new peice from Elliot
                           data = d_full_wo_gaps, family = binomial)
summary(gamm4_geo$gam)
plot(gamm4_geo$gam)


gamm4_ProvCode <- gamm4(cbind(NwP, N-NwP) ~ ProvCode,  # swap longhurst region with mean poll abund (try spline and no spline) instead of this, and adjacency (no spline)
                   random = ~(1|order) + (1|source), 
                   REML = TRUE,
                   data = d_full_wo_gaps, family = binomial)
summary(gamm4_ProvCode$gam)
plot(gamm4_ProvCode$gam)


#testing publication year
gamm4_FwP_pubyear <- gamm4(cbind(NwP, N-NwP) ~ s(publication_year),
                          random = ~(1|order) + (1|source), 
                           data = d_full_wo_gaps, family = binomial)
summary(gamm4_FwP_pubyear$gam)
summary(gamm4_FwP_pubyear$dev.expl)
plot(gamm4_FwP_pubyear$gam)


# null model
gamm4_null <- gamm4(cbind(NwP, N-NwP) ~ 1, 
                   random = ~(1|order) + (1|source), 
                   data = d_full_wo_gaps, family = binomial)
summary(gamm4_null$gam)
plot(gamm4_null$gam)



# Trying a GLMM

glmm_FwP_eco_geo <- glmer(cbind(NwP, N-NwP) ~ trophic_level_via_fishbase + 
                            scale(average_depth) + 
                            scale(mean_poll_abund) +
                            Found + prime_forage + adjacency +
                            (1|order) + (1|source), 
                          na.action = "na.fail",
                          data = d_full_wo_gaps, family = binomial)
summary(glmm_FwP_eco_geo)
r.squaredGLMM(glmm_FwP_eco_geo)



glmm_FwP_eco <- glmer(cbind(NwP, N-NwP) ~ trophic_level_via_fishbase + 
                            scale(average_depth) + 
                            Found + prime_forage +
                            (1|order) + (1|source), 
                          na.action = "na.fail",
                          data = d_full_wo_gaps, family = binomial)
summary(glmm_FwP_eco)
r.squaredGLMM(glmm_FwP_eco)


glmm_FwP_geo <- glmer(cbind(NwP, N-NwP) ~ scale(mean_poll_abund) + adjacency +
                            (1|order) + (1|source), 
                          na.action = "na.fail",
                          data = d_full_wo_gaps, family = binomial)
summary(glmm_FwP_geo)
r.squaredGLMM(glmm_FwP_geo) # NEED TO DO THIS TO TEST FOR FIT



# a GLM with no random effects
glm_FwP_eco_geo <- glm(cbind(NwP, N-NwP) ~ trophic_level_via_fishbase + 
                         scale(average_depth) + 
                         scale(mean_poll_abund) +
                         Found + prime_forage + adjacency,
                       data = d_full_wo_gaps, family = binomial)
summary(glm_FwP_eco_geo)

# multi-model selection using AICc
GLMM_dredge <- dredge(glmm_FwP_eco_geo)

View(GLMM_dredge)
write_csv(GLMM_dredge, "GLMM model selection table.csv")
#subset(GLMM_dredge, delta < 4)
a=model.avg(GLMM_dredge)
summary(a) #The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. An alternative, the ‘full’ average assumes that a variable is included in every model
confint(a) #computes confidence interval








## Trying a BRT
d_full_wo_gaps <- d_full %>% 
  drop_na(order,trophic_level_via_fishbase, habitat, prime_forage, average_depth, oceanographic_province_from_longhurst_2007)

# d_test <- d_full %>%  
#   slice(1:100) %>% 
#   select(prop_w_plastic, order,trophic_level_via_fishbase, habitat, prime_forage, average_depth, oceanographic_province_from_longhurst_2007) %>% 
#   drop_na() %>% 
#   mutate(prop_w_plastic = numeric(prop_w_plastic),
#          trophic_level_via_fishbase = numeric(trophic_level_via_fishbase),
#          average_depth = numeric(average_depth))

gbmFwP <- gbm.step(data=as.data.frame(d_full), 
                   gbm.x = 7, 13,  
                   gbm.y = 10,   # this is Prop w plastic, 34 is Prop w plastic multiplied by the assessment's sample size
                   #weights = 9,  # weighted by sample size
                   family = "gaussian", 
                   tree.complexity = 5,
                   learning.rate = 0.001, bag.fraction = 0.5)
summary(gbmFwP)
gbm.plot(gbmFwP)
