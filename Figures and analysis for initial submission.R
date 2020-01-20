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
library(sjPlot)
library(tibble)
library(tidyr)
library(ggeffects)

SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}


# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}

d_poll <- as_tibble(read_csv("Spatial Information_microplastics.csv"))

d = read_csv("Plastics ingestion records fish master_final.csv") %>% 
  janitor::clean_names() %>% 
  rename(NwP = nw_p,
         N = n) %>% 
  mutate(WeightedProp = prop_w_plastic * N,
         Found = as_factor(case_when(habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                                     habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")),
         prime_forage = na_if(prime_forage, "not listed")) %>% 
  rename(ProvCode = "oceanographic_province_from_longhurst_2007") %>% 
  separate(binomial, into = c("genus", "species"), sep = " ", remove = FALSE) %>% 
  left_join(dplyr::select(d_poll, c(ProvCode, adjacency, mean_poll_abund)), by = "ProvCode") %>% 
  mutate(adjacency = as_factor(case_when(adjacency == 1 ~ "coastal",
                                         adjacency == 0 ~ "oceanic")),
         source = as_factor(source))
           

d_full <- d %>%
  filter(includes_microplastic == "Y")


# summary tables----

# species summary table
d_sp_sum <- d %>%
  filter(!species %in% c("sp.", "spp.")) %>%
  group_by(binomial, order, commercial, iucn_status) %>%
  drop_na(binomial) %>% 
  summarize(Sp_mean = mean(prop_w_plastic, na.rm = TRUE),
            Sample_size = sum(N),
            num_studies = n_distinct(source)) %>% 
  ungroup %>% 
  mutate(commercial = factor(commercial),
         studies_cat = as.double(cut(num_studies, 
                                     c(0, 1, 3, Inf),
                                     c(1,2,3))),
         commercial = fct_collapse(commercial,
                                   Commercial = c("commercial", "highly commercial"),
                                   Minor = c("minor commercial", "subsistence"),
                                   None = "none")) %>%
  #filter(commercial == "Commercial") %>% 
  arrange(-Sp_mean)



d_family_disc <- d %>% 
  filter(family %in%  c("Carangidae", "Mugilidae", "Pleuronectidae", "Soleidae", "Myctophidae")) %>% 
  group_by(family) %>% 
  summarise(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N, na.rm = TRUE),
            num_studies = n_distinct(source)) %>% 
  arrange(desc(FO_plastic))
write_csv(d_family_disc, "Fish families of concern.csv")


# fish of concern for humans
concern_fish <- d %>% 
  group_by(common_name, binomial, family) %>% 
  filter(commercial %in% c("commercial", "highly commercial")) %>%
  summarize(species_avg = sum(NwP)/sum(N),
            mean_num_plast = weighted.mean(mean_num_particles_per_indv, N),
            sample_size = sum(N),
            num_studies = n_distinct(source),
            commercial_status = first(commercial),
            AC_status = first(aquaculture), 
            rec_status = first(recreational)) %>% 
  ungroup %>% 
  #filter(species_avg > 0.25 & sample_size > 25) %>% 
  arrange(-species_avg)
write_csv(concern_fish, "Concerning fish for humans.csv")

# IUCN and vulnerability status of fish
conserve_fish <- d %>% 
  group_by(binomial) %>% 
  #filter(iucn_status == "LC") %>% 
  filter(iucn_status %in% c("NT","VU", "EN", "CR")) %>% 
  summarize(sp_avg = sum(NwP)/sum(N),
            iucn = first(iucn_status),
            sample_size = sum(N)) %>% 
  arrange(-sp_avg)

d_vulnerability <- d_full %>%
  #  filter(family %in% c("Scombridae", "Sphyrnidae", "Carcharhinidae")) %>% 
  drop_na(vulnerability_score_via_fishbase_from_cheug_et_al_2005) %>% 
  mutate(Vulnerability.category = cut(vulnerability_score_via_fishbase_from_cheug_et_al_2005,
                                      breaks=c(-Inf, 20, 40, 60, 80, Inf), 
                                      labels=c("low","low-moderate", "moderate-high", "high-very high", "very high"))) %>% 
  group_by(common_name, binomial, Vulnerability.category, vulnerability_score_via_fishbase_from_cheug_et_al_2005, iucn_status) %>%  
  summarize(FO_plastic = sum(NwP)/sum(N),
            Sample_size = sum(N),
            #    num_sp = n_distinct(binomial), 
            num_studies = n_distinct(source)) %>% 
  ungroup %>% 
  filter(Vulnerability.category %in% c("moderate-high", "high-very high", "very high") & 
           #   iucn_status %in% c("NT","VU", "EN", "CR") & 
           FO_plastic > 0.25 & Sample_size > 25)
write_csv(d_vulnerability, "Vulnerability table.csv")

# geographic summary of data
Fish_geo_summ <- d_full %>% 
  filter(ProvCode %in% c("CHIN", "KURO", "SUND", "INDE")) %>%
  #filter(ProvCode %in% c("NAST E", "NECS", "MEDI")) %>% 
  #filter(ProvCode == "BPRL") %>% 
  group_by(ProvCode) %>% 
  summarize(num_studies = n_distinct(source),
            num_sp = n_distinct(binomial),
            num_w_plast = sum(NwP, na.rm = TRUE),
            num_ind_studied = sum(N, na.rm = TRUE),
            prop_by_region = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            wgt_mean_plast_num = weighted.mean(mean_num_particles_per_indv, N),
            se_plast_num = SE(mean_num_particles_per_indv))
write_csv(Fish_geo_summ, "Fish_plastic geographic summary.csv")




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

# SUPPLEMENTAL FIGURE, ALL SPECIES
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
  scale_shape_discrete(na.translate = F) +
  labs(shape = "Commercial \nstatus") +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(3, "cm"),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 28),
        legend.box = "horizontal") +
  guides(shape = guide_legend(override.aes = list(size = 5)))


# save plots
dev.copy2pdf(file="Fish_plastic_phylo.pdf", width=50, height=50)
ggsave("Prelim_phylo_ggtree2.tiff", width = 42, height = 42, units = "in")
ggsave("Prelim_phylo_ggtree2.eps", width = 45, height = 45, units = "in")
ggsave("Prelim_phylo_ggtree2.jpg", width = 45, height = 45, units = "in")

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
                        labels = c("Poorly studied (n=1)", "Moderately studied (n=2-3)", "Well studied (n=4-9)"),
                        range = c(1.5, 5)) +
  annotate("text", x = c(0.6, 2.8, 0.4, 3.5),
           y=c(0.8, 0.8, 0.08, 0.08),
           label = c("high incidence, data poor", "high incidence, data rich",
                     "low incidence, data poor", "low incidence, data rich")) +
  theme_classic(base_size = 16)
risk_plot

dev.copy2pdf(file="risk_plot.pdf", width=12, height=7)


# Figure 3, plastic ingestion by depth and habitat ----

p <- drop_na(d_full, trophic_level_via_fishbase, prop_w_plastic, N)

p <-   ggplot(filter(d_full, order == "Perciformes"),
    aes(trophic_level_via_fishbase, 
             y = prop_w_plastic, 
             size = N, 
             weight = N)) + 
  geom_point(alpha = 0.4) +  # Eventually add in foraging behavior here 
  geom_smooth(method = "lm", se = FALSE) +
  #geom_smooth(aes(color = source), method = "lm", se = FALSE) +
  #xlim(2, 5) +
  xlab("Trophic level") +
  ylab("Proportion of individuals with ingested plastic") +
  #facet_wrap(~order) +
  # annotate("text", x = c(2.75, 4), y= -0.05,    # Comment out if faceting
  #          label = c("planktivorous", "piscivorous")) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold")) 
p + guides(size = FALSE)

# DO A GLM 


# Figure 1, Family phylogeny----

# to convert phylo object to Newick string
writeTree<-function(tree){
  tree<-reorder.phylo(tree,"cladewise")
  n<-length(tree$tip)
  string<-vector(); string[1]<-"("; j<-2
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,2]<=n){
      string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
      if(!is.null(tree$edge.length)){
        string[j]<-paste(c(":",round(tree$edge.length[i],10)), collapse="")
        j<-j+1
      }
      v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
      while(length(v)>0&&k==v[length(v)]){
        string[j]<-")"; j<-j+1
        w<-which(tree$edge[,2]==tree$edge[k,1])
        if(!is.null(tree$edge.length)){
          string[j]<-paste(c(":",round(tree$edge.length[w],10)), collapse="")
          j<-j+1
        }
        v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
      }
      string[j]<-","; j<-j+1
    } else if(tree$edge[i,2]>=n){
      string[j]<-"("; j<-j+1
    }
  }
  if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)], ";")
  else string<-c(string[1:(length(string)-2)],";")
  string<-paste(string,collapse="")
  return(string)
}




d_family <- d_full %>% 
  group_by(family) %>% 
  summarise(FO_plastic = sum(NwP, na.rm = TRUE)/sum(N, na.rm = TRUE),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N, na.rm = TRUE),
            num_studies = n_distinct(source),
            prop_commercial = sum(commercial %in% c("commercial", "highly commercial"))/n()) %>% 
  drop_na(family) %>% 
  mutate(studies_cat = as.double(cut(num_studies, 
                                     c(0, 1, 3, Inf),
                                     c(1,2,3))),
         commercial_cat = cut(prop_commercial,
                              breaks=c(-Inf, 0.01, 0.25, Inf), 
                              labels=c("None", "Minor", "Commercial"))) %>% 
  arrange(desc(FO_plastic))



breaks <- c(seq(1,nrow(d_family),50),nrow(d_family)+1)  # why are we doing this?

for (i in 1:(length(breaks)-1)){
  taxa <- as.character(d_family$family[breaks[i]:(breaks[i+1]-1)])
  taxa <- taxa[taxa != "" & !is.na(taxa)]
  
  resolved_namest <- tnrs_match_names(taxa, context_name = "Animals")                          # I think this is where all the extra species fall out
  resolved_namest <- resolved_namest[!is.na(resolved_namest$unique_name),]
  if (i==1){
    resolved_namess <- resolved_namest
  } else {
    resolved_namess <- rbind(resolved_namess, resolved_namest)
  }
}
resolved_names <- resolved_namess
resolved_names <- resolved_names[resolved_names$flags!="INCERTAE_SEDIS_INHERITED",]


#toString(fish_familes$unique_name)  #puts commas between all the numbers for below

#dput(fish_familes$unique_name) #puts commas between all the names with quotes for below

#plots FAMILIES
my_tree <- tol_induced_subtree(ott_ids = c(85219, 118776, 816156, 19305, 893503, 960231, 563241, 253371, 486695, 42415, 738324, 561376, 239745, 679814, 749780, 34184, 875698, 
                                           978560, 1003121, 384676, 137651, 739933, 191025, 646019, 804461, 722754, 1089734, 637234, 648754, 37461, 888446, 563518, 856584,
                                           186486, 563230, 479853, 734459, 437596, 1089730, 441564, 292707, 823203, 427507, 42408, 655424, 44805, 577720, 19301, 441571, 232621, 
                                           279762, 1089742, 813991, 1032209, 129124, 804451, 308750, 587923, 583638, 340592, 400235, 749770, 1074732, 701559, 563513, 710014, 
                                           99942, 450143, 609233, 160291, 712841, 739939, 883406, 1089740, 460871, 734202, 818997, 214115, 132684, 72393, 407171, 563254, 74014, 
                                           280947, 99937, 769569, 214115, 175436, 308741, 114163, 363181, 615333, 258647, 1052881, 13838, 555245, 1093612, 930712, 778687, 579429, 
                                           562630, 574822, 280953, 339090, 540474, 859881, 137656, 62639, 811925, 479864, 644001, 95055, 401063, 765787, 715629, 724438, 659851, 
                                           614519, 769567, 892951, 892958, 237354, 250743, 65336, 553102, 1026498, 36225), 
                               label_format = "name")

Gadiformes_tree <- read.tree("subtree-node-mrcaott33200ott659845-Gaidropsaridae--Raniceps.tre")
final_tree <- bind.tree(my_tree + Gadiformes_tree)

Gnathostomata_tree <- read.tree("subtree-node-ott278114-Gnathostomata.tre")



# Adding tree from Rabosky et al. 2018 An inverse latitudinal gradient in speciation rate for marine fishes. Nature
# x <- read.tree("actinopt_12k_treePL.tre.xz")
# New_x <- writeTree(x)
# y <- my_tree
# 
# final_tree <- bind.tree(x, y, position = if (is.null(x$root.edge)) 0 else
#   x$root.edge)
# 
# my_tree <- tol_induced_subtree(ott_ids = c(85219, 118776, 816156, 19305, 893503, 960231, 563241, 253371, 486695, 42415, 738324, 561376, 239745, 679814, 749780, 34184, 875698, 
#                                            978560, 1003121, 384676, 137651, 739933, 191025, 646019, 804461, 722754, 1089734, 637234, 648754, 37461, 888446, 563518, 856584,
#                                            186486, 563230, 479853, 734459, 437596, 1089730, 441564, 292707, 823203, 427507, 42408, 655424, 44805, 577720, 19301, 441571, 232621, 
#                                            279762, 1089742, 813991, 1032209, 129124, 804451, 308750, 587923, 583638, 340592, 400235, 749770, 1074732, 701559, 563513, 710014, 
#                                            99942, 450143, 609233, 160291, 712841, 739939, 883406, 1089740, 460871, 734202, 818997, 214115, 132684, 72393, 407171, 563254, 74014, 
#                                            280947, 99937, 769569, 214115, 175436, 308741, 114163, 363181, 615333, 258647, 1052881, 13838, 555245, 1093612, 930712, 778687, 579429, 
#                                            562630, 574822, 280953, 339090, 540474, 859881, 137656, 62639, 811925, 479864, 644001, 95055, 401063, 765787, 715629, 724438, 659851, 
#                                            614519, 769567, 892951, 892958, 237354, 250743, 65336, 553102, 1026498, 36225), 
#                                file = writeTree(x))
# 
# 

# my_tree$tip.label<-gsub("_"," ",my_tree$tip.label) # removes underscore between genus and species names
# my_tree$tip.label<-str_extract(my_tree$tip.label, "[A-Z][a-z]+ [a-z]+")

my_tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
my_tree <- ladderize(my_tree, right = TRUE)

View(my_tree)

final_clipped_tree <- keep.tip(Gnathostomata_tree,
                               tip = dput(fish_familes$unique_name))


#my_tree2 = phylo4d(my_tree, d_sp_sum)

my_tree$tip.label[my_tree$tip.label == "Belonidae"] <- "Scomberesocidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott801ott480916"] <- "Ariidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott595ott3629079"] <- "Sternoptychidae"
my_tree$tip.label[my_tree$tip.label == "Gonostomatidae_(family_in_Opisthokonta)"] <- "Gonostomatidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott47507ott47509"] <- "Acanthuridae"
my_tree$tip.label[my_tree$tip.label == "mrcaott7529ott77068"] <- "Cottidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott9658ott106653"] <- "Pholidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott15200ott303035"] <- "Nototheniidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott126948ott301378"] <- "Cheilodactylidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott7815ott9211"] <- "Gobiidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott61966ott80330"] <- "Bythitidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott56885ott4134216"] <- "Lotidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott13841ott182974"] <- "Nettastomatidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott65747ott167810"] <- "Ophichthidae"
my_tree$tip.label[my_tree$tip.label == "mrcaott73302ott93287"] <- "Atherinopsidae"
my_tree$tip.label[my_tree$tip.label == "Alepocephaliformes"] <- "Alepocephalidae"
#my_tree$tip.label[my_tree$tip.label == "mrcaott33200ott659845"] <- "Gadidae"

#bind.tip(my_tree, tip.label = "Gadidae", position = "mrcaott33200ott659845")

# first plot try, fan layout
p <- ggtree(my_tree, layout="fan", open.angle=0) + 
  #geom_text2(aes(label=label), hjust=-.2, size=4) +
  ggplot2::xlim(-0.6, 1.3) 
p



# Adding data
shapes <- c("None" = 15, "Minor" = 17, "Commercial" = 16)

p %<+% d_family + 
  aes(color = FO_plastic) +
  geom_tiplab2(aes(label = paste0("italic('", label, "')"),
                   color=FO_plastic), parse = TRUE,
               size = 6, align = FALSE, offset = 0.05) +
  geom_tippoint(aes(color = FO_plastic, shape = commercial_cat, size = studies_cat)) +
  scale_color_gradientn(colours = c("steelblue4",
                                    "gray40",
                                    "coral", "coral1",
                                    "firebrick2", "firebrick3", "firebrick4"), 
                        name = "Proportion with \ningested plastic") +
  #scale_size(range = c(3, 7)) +
  scale_size_continuous(guide = FALSE, range = c(3, 7)) +
  scale_shape_manual(na.translate = F, values = shapes) +
  labs(shape = "Commercial \nstatus") +
  theme(legend.position = c(0.49, 0.49),
        legend.key.size = unit(1.25, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.box = "horizontal") +
  guides(shape = guide_legend(override.aes = list(size = 5)))

dev.copy2pdf(file="Fish_family_plastic_phylo.pdf", width=20, height=20)
ggsave("Fish_family_plastic_phylo.jpg", width = 20, height = 20, units = "in")






m1 <- lm(trophic_level_via_fishbase~prop_w_plastic, data = d_full)
summary(m1)

Fig_3 <- ggplot(filter(d_full, Found != "NA"), 
               aes(average_depth, prop_w_plastic, size = N, weight = N)) +
  geom_point(alpha = 0.4) + 
  geom_smooth(col = "blue", method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Found, scales = "free_y", ncol = 1) +
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


Fig_3_PF <- ggplot(filter(d_full, prime_forage != "NA"), 
                aes(prime_forage, prop_w_plastic, size = N, weight = N)) +
  geom_boxplot() +
  geom_jitter(width = 0.3, alpha = 0.4) + 
  xlab("Primary foraging behaivor") +
  ylab("Proportion of individuals with ingested plastic") +
  theme_classic(base_size = 16)
# theme(axis.text.x = element_text(size=12),
#       axis.text.y = element_text(size=12),
#       axis.title=element_text(size=13, face="bold"),
#       strip.text = element_text(size = 12)) 
Fig_3_PF + guides(size = FALSE)



Fig_3_PA <- ggplot(filter(d_full, primary_aggregation %in% c("Schooling", "Solitary")), 
                   aes(primary_aggregation, prop_w_plastic)) +
  geom_boxplot() +
  geom_jitter(width = 0.3, alpha = 0.4) + 
  xlab("Aggregation behaivor") +
  ylab("Proportion of individuals with ingested plastic") +
  theme_classic(base_size = 16)
# theme(axis.text.x = element_text(size=12),
#       axis.text.y = element_text(size=12),
#       axis.title=element_text(size=13, face="bold"),
#       strip.text = element_text(size = 12)) 
Fig_3_PA + guides(size = FALSE)


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
            annual_N = sum(N, na.rm = TRUE)) %>% 
  ungroup %>% 
  mutate(cum_species = cum_unique(species),
         cum_n = cumsum(annual_N),
         cpue = cum_species / cum_n)
d_rarefaction_ingest_only <-  d_full %>%  
  filter(prop_w_plastic > 0) %>% 
  group_by(publication_year) %>% 
  summarize(species = list(unique(binomial)),
            annual_N = sum(N, na.rm = TRUE)) %>% 
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

dev.copy2pdf(file="rarefaction_plot.pdf", width=8, height=8)




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

dev.copy2pdf(file="studies_by_year.pdf", width=12, height=8)

######################################################### 
# Modeling for fish-plastic ingestion meta-analysis paper
#########################################################

d_full_wo_gaps_TL <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(average_depth, Found, trophic_level_via_fishbase, NwP, N, 
          adjacency, mean_poll_abund, ProvCode) 

#Full GLMM with trophic level----

glmm_FwP_eco_geo_TL <- glmer(cbind(NwP, N-NwP) ~ scale(trophic_level_via_fishbase) +
                            scale(average_depth) + 
                            scale(mean_poll_abund) +
                            Found + adjacency +
                            (1|order) + (1|source), 
                          na.action = "na.fail",
                          data = d_full_wo_gaps_TL, family = binomial)
summary(glmm_FwP_eco_geo_TL)
r.squaredGLMM(glmm_FwP_eco_geo_TL) # use R2c theoretical value, see: https://www.rdocumentation.org/packages/MuMIn/versions/1.43.15/topics/r.squaredGLMM  

# DO PLOT WITH DIFFERENT LINES FOR DIFFERENT ORDER, TL and Plastic ingestion by ORDER
# DO TROPHIC LEVEL AND FORAGING MODE IN SEPARATE MODELS SINCE THEY"RE CONFLATIED
# LOOK AT INTERCEPT ONLY GLMMs, DO RANDOM EFFECTS MAKE SENSE?
# PLOT PRIMARY FORAGING MODE AS SEPARATE MODEL WITH IT AS MAIN OR RANDOM EFFECT 

dat <- ggpredict(glmm_FwP_eco_geo_TL, terms = "trophic_level_via_fishbase")
plot(dat) +
  ggtitle("") +
  xlab("Trophic level") +
  ylab("Plastic ingestion incidence") +
  ylim(0.15,0.6) +
annotate("text", x = c(2.5, 4), y= 0.16, 
         label = c("planktivorous", "piscivorous"))

dev.copy2pdf(file="TL_GLMM response.pdf", width=4.5, height=5)

dat2 <- ggpredict(glmm_FwP_eco_geo_TL, terms = "Found")
plot(dat2) +
  ggtitle("") +
  xlab("Habitat") +
  ylab("Plastic ingestion incidence") +
  ylim(0.1,0.6)

dev.copy2pdf(file="Found_GLMM response.pdf", width=4.5, height=5)

dat2 <- ggpredict(glmm_FwP_eco_geo_TL, terms = "adjacency")
plot(dat2) +
  ggtitle("") +
  xlab("Proximity to continent") +
  ylab("Plastic ingestion incidence") +
  ylim(0,0.6)

dev.copy2pdf(file="Adjacency_GLMM response.pdf", width=4.5, height=5)

a=get_model_data(glmm_FwP_eco_geo_TL, type = "eff")[[2]] 

# plot total model predictors on response varaible 
plot_model(glmm_FwP_eco_geo_TL, sort.est = TRUE) +
  ggtitle("Response: Proportion with ingested plastic") +
  xlab("Predictors") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 1) +
  ylim(0.001,2) +
  scale_x_discrete(labels=c("scale(mean_poll_abund)" = "plastic pollution abundance",
                            "scale(average_depth)" = "average depth found",
                            "scale(trophic_level_via_fishbase)" = "trophic level",
                            "Founddemersal" = "habitat: demersal",
                            "adjacencyoceanic" = "proximity: oceanic"))
dev.copy2pdf(file="Main effects TL_GLMM.pdf", width=7, height=4)

# # plot main effects WOW
# plot_model(glmm_FwP_eco_geo_TL, type = "eff", terms = "adjacency") +
#   ggtitle("") +
#   ylab("Proportion with ingested plastic") +
#   scale_x_discrete(breaks=c("ctrl", "trt1", "trt2"),
#                         labels=c("Control", "Treat 1", "Treat 2"))

# plot random effects WOW
plot_model(glmm_FwP_eco_geo_TL, type = "re",
           grid = FALSE, 
           sort.est = "sort.all")[[2]] +
  ggtitle("") +
  xlab("Order") +
  ylab("Random intercept") +
  #ylim(0.1,10) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = 1)

dev.copy2pdf(file="Random effect of order_GLMM.pdf", width=7, height=11)


 # multi-model selection using AICc
GLMM_dredge <- dredge(glmm_FwP_eco_geo_TL)

View(GLMM_dredge)
write_csv(GLMM_dredge, "GLMM model selection table.csv")
#subset(GLMM_dredge, delta < 4)
a=model.avg(GLMM_dredge)
summary(a) #The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. An alternative, the ‘full’ average assumes that a variable is included in every model
confint(a) #computes confidence interval


#Full GLMM with foraging behavior----

d_full_wo_gaps_PF <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(average_depth, Found, NwP, N, prime_forage,
          adjacency, mean_poll_abund, ProvCode) 

glmm_FwP_eco_geo_PF <- glmer(cbind(NwP, N-NwP) ~ scale(average_depth) + 
                               scale(mean_poll_abund) +
                               prime_forage + Found + adjacency +
                               (1|order) + (1|source), 
                             na.action = "na.fail",
                             data = d_full_wo_gaps_PF, family = binomial)
summary(glmm_FwP_eco_geo_PF)
r.squaredGLMM(glmm_FwP_eco_geo_PF)

# DO PLOT WITH DIFFERENT LINES FOR DIFFERENT ORDER, TL and Plastic ingestion by ORDER
# DO TROPHIC LEVEL AND FORAGING MODE IN SEPARATE MODELS SINCE THEY"RE CONFLATIED
# LOOK AT INTERCEPT ONLY GLMMs, DO RANDOM EFFECTS MAKE SENSE?
# PLOT PRIMARY FORAGING MODE AS SEPARATE MODEL WITH IT AS MAIN OR RANDOM EFFECT 


dat2 <- ggpredict(glmm_FwP_eco_geo_PF, terms = "prime_forage")
plot(dat2) +
  ggtitle("") +
  xlab("Primary foraging behavior") +
  ylab("Proportion with ingested plastic") +
  ylim(0,0.6)




glmm_FwP_eco <- glmer(cbind(NwP, N-NwP) ~ trophic_level_via_fishbase + 
                            scale(average_depth) + 
                            Found + 
                            (1|order) + (1|source), 
                          na.action = "na.fail",
                          data = d_full_wo_gaps_TL, family = binomial)
summary(glmm_FwP_eco)
r.squaredGLMM(glmm_FwP_eco) # use R2c theoretical value, see: https://www.rdocumentation.org/packages/MuMIn/versions/1.43.15/topics/r.squaredGLMM  


d_full_wo_gaps_geo <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(NwP, N, adjacency, mean_poll_abund) 

glmm_FwP_geo <- glmer(cbind(NwP, N-NwP) ~ scale(mean_poll_abund) + adjacency +
                            (1|order) + (1|source), 
                          na.action = "na.fail",
                          data = d_full_wo_gaps_geo, family = binomial)
summary(glmm_FwP_geo)
r.squaredGLMM(glmm_FwP_geo) # use R2c theoretical value, see: https://www.rdocumentation.org/packages/MuMIn/versions/1.43.15/topics/r.squaredGLMM  



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
d_full_droprespNA <- d_full %>% 
  drop_na(prop_w_plastic)

  drop_na(order,trophic_level_via_fishbase, habitat, prime_forage, average_depth, oceanographic_province_from_longhurst_2007)

# d_test <- d_full %>%  
#   slice(1:100) %>% 
#   select(prop_w_plastic, order,trophic_level_via_fishbase, habitat, prime_forage, average_depth, oceanographic_province_from_longhurst_2007) %>% 
#   drop_na() %>% 
#   mutate(prop_w_plastic = numeric(prop_w_plastic),
#          trophic_level_via_fishbase = numeric(trophic_level_via_fishbase),
#          average_depth = numeric(average_depth))

gbmFwP <- gbm.step(data=as.data.frame(d_full_droprespNA), 
                   gbm.x = 36, 37,  
                   gbm.y = 10,   # this is Prop w plastic, 34 is Prop w plastic multiplied by the assessment's sample size
                   #weights = 9,  # weighted by sample size
                   family = "gaussian", 
                   tree.complexity = 5,
                   learning.rate = 0.001, bag.fraction = 0.5)
summary(gbmFwP)
gbm.plot(gbmFwP)
