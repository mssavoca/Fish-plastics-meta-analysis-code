########################################################################################
# Preliminary visualizations and analyses for marine fish-plastic ingestion meta-analysis
########################################################################################

# load packages and data ----
library(dismo)
library(gbm)
library(gamm4)
library(ggrepel)
library(itsadug)
library(lme4)
library(MCMCglmm)
library(mgcv)
library(MuMIn)
library(readxl)
library(readr)
library(tidyverse)

SE = function(x){sd(x)/sqrt(sum(!is.na(x)))}

d = read_csv("Plastics ingestion records fish master_final.csv") %>% 
  janitor::clean_names() %>% 
  rename(NwP = nw_p,
         N = n) %>% 
  mutate(WeightedProp = prop_w_plastic * N,
         Found = as_factor(case_when(habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                                     habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")),
         prime_forage = na_if(prime_forage, "not listed")) %>% 
  separate(binomial, into = c("genus", "species"), sep = " ", remove = FALSE)

d_full <- d %>%
  filter(includes_microplastic == "Y")


# d_old_data = read_csv("Plastics ingestion records fish master_updated.csv") %>% 
#   janitor::clean_names() %>% 
#   rename(NwP = nw_p,
#          N = n) %>% 
#   mutate(WeightedProp = prop_w_plastic * N,
#          Found = as_factor(case_when(habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
#                                      habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")),
#          prime_forage = na_if(prime_forage, "not listed")) %>% 
#   separate(species_name, into = c("genus", "species"), sep = " ", remove = FALSE)
# 
# 
# d_old_data <- d_old_data %>% 
#   drop_na(average_depth, Found, trophic_level_via_fishbase, prime_forage) 
# 
# anti_join(distinct(d_full, source), distinct(d_old_data, source))


# 
# d1 = read_xlsx("Plastics ingestion records fish master_updated.xlsx") %>% 
#   mutate(Found = case_when(Habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
#                            Habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")) %>% 
#   filter(!is.na(`Prop w plastic`))

# average number of plastic particles per individual   
d %>% filter(prop_w_plastic > 0) %>% 
  drop_na(mean_num_particles_per_indv) %>% 
  summarise(avg = mean(mean_num_particles_per_indv), se = SE(mean_num_particles_per_indv)) 

# summary tables ----
# some summary tables
d_sp_sum <- d_full %>%
  filter(!species %in% c("sp.", "spp.")) %>%
  #filter(commercial %in% c("commercial", "highly commercial")) %>% 
  group_by(binomial, order, commercial) %>%
  summarize(Sp_mean = mean(prop_w_plastic),
            Sample_size = sum(N),
            num_studies = n_distinct(source)) %>%
  ungroup %>% 
  arrange(-Sp_mean)


d_phylo_summ <- d_full %>% 
  group_by(order) %>% 
  drop_na(order) %>% 
  summarise(FO_plastic = sum(NwP)/sum(N),
            # sp = first(species),
            # TL = first(trophic_level_via_fishbase),
            wgt_mean_plast_num = mean(mean_num_particles_per_indv),
            se_num_plast = SE(mean_num_particles_per_indv),
            sample_size = sum(N),
            num_studies = n_distinct(source)) %>% 
#  filter(FO_plastic > 0.1) %>% 
  arrange(desc(FO_plastic))


d_family_disc <- d %>% 
  filter(family %in%  c("Carangidae", "Mugilidae", "Pleuronectidae", "Soleidae", "Myctophidae")) %>% 
  group_by(family) %>% 
  summarise(FO_plastic = sum(NwP)/sum(N),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N),
            num_studies = n_distinct(source)) %>% 
  arrange(desc(FO_plastic))
write_csv(d_family_disc, "Fish families of concern.csv")

d_sp_disc <- d_full %>% 
  filter(binomial %in%  c("Thunnus thynnus", "Thunnus alalunga", "Katsuwonus pelamis", "Thunnus obesus", "Thunnus albacares",
                              "Galeocerdo cuvier", "Prionace glauca")) %>% 
  mutate(Vulnerability.category = cut(vulnerability_score_via_fishbase_from_cheug_et_al_2005,
                                      breaks=c(-Inf, 20, 40, 60, 80, Inf), 
                                      labels=c("low","low-moderate", "moderate-high", "high-very high", "very high"))) %>%
  group_by(common_name, binomial) %>% 
  summarise(FO_plastic = sum(NwP)/sum(N),
            mean_num_plast = mean(mean_num_particles_per_indv, na.rm = TRUE),
            sample_size = sum(N),
            num_studies = n_distinct(source),
            commercial_status = first(commercial),
            AC_status = first(aquaculture), 
            rec_status = first(recreational),
            Vulnerability.score = first(vulnerability_score_via_fishbase_from_cheug_et_al_2005),
            IUCN.status = first(iucn_status)) %>% 
  ungroup %>% 
  arrange(desc(FO_plastic))
write_csv(d_sp_disc, "Tuna and shark species of concern.csv")

conserve_overview_fish <- d %>% 
  group_by(iucn_status) %>% 
  summarize(num_sp_per_cat = n_distinct(binomial))

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


# fish of concern for humans
concern_fish <- d %>% 
  group_by(common_name, binomial, family, species) %>% 
  filter(commercial %in% c("commercial", "highly commercial")) %>%
  summarize(species_avg = sum(NwP)/sum(N),
            mean_num_plast = weighted.mean(mean_num_particles_per_indv, N),
            sample_size = sum(N),
            num_studies = n_distinct(source),
            commercial_status = first(commercial),
            AC_status = first(aquaculture), 
            rec_status = first(recreational)) %>% 
  ungroup %>% 
#  filter(species_avg > 0.25 & sample_size > 25) %>% 
  arrange(-species_avg)
write_csv(concern_fish, "Concerning fish for humans.csv")

# geographic summary of data
Fish_geo_summ <- d %>% 
  #filter(oceanographic_province_from_longhurst_2007 %in% c("CHIN", "KURO", "SUND", "INDE")) %>%
  #filter(oceanographic_province_from_longhurst_2007 %in% c("NAST E", "NECS", "MEDI")) %>% 
  #filter(oceanographic_province_from_longhurst_2007 == "BPRL") %>% 
  group_by(oceanographic_province_from_longhurst_2007) %>% 
  summarize(num_studies = n_distinct(source),
            num_sp = n_distinct(binomial),
            num_w_plast = sum(NwP),
            num_ind_studied = sum(N),
            prop_by_region = sum(NwP)/sum(N),
            wgt_mean_plast_num = weighted.mean(mean_num_particles_per_indv, N),
            se_plast_num = SE(mean_num_particles_per_indv))
write_csv(Fish_geo_summ, "Fish_plastic geographic summary.csv")

# preliminary plots ----

# MAX'S RAREFACTION CURVE
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


# Supplemental plot by 5 year bins
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


FO_by_year <- d %>% 
    filter(source != "Lucia et al. 2018") %>% 
    group_by(publication_year) %>% 
    summarize(FO = sum(NwP)/sum(N)) %>% 
    ggplot(aes(publication_year, FO)) + 
    geom_bar(stat = "identity") + 
    geom_smooth(se = FALSE) +
    theme_classic() +
    xlab("Publication year") + 
    ylab("Frequency of occurence")
FO_by_year 

num_plast_part_yr <- d %>% 
  filter(includes_microplastic == "Y") %>% 
  group_by(publication_year) %>% 
  summarize(mean_num_part = mean(mean_num_particles_per_indv, na.rm = TRUE)) %>% 
  ggplot(aes(publication_year, mean_num_part)) + 
  geom_bar(stat = "identity") + 
  geom_smooth(se = FALSE) +
  theme_classic() +
  xlab("Publication year") + 
  ylab("Mean number of particles ingested per year")
num_plast_part_yr

# map plot, THIS WORKS!!
require(raster)
shape <- shapefile("longhurst_v4_2010/Longhurst_world_v4_2010.shp")
plot(shape)

name <- shape$ProvCode
polygon <- shape@polygons


# risk/interest plot, do color by order, shape by Commercial
risk_plot <- d_sp_sum %>% 
#  filter(!species %in% c("sp.", "spp.")) %>%
  group_by(binomial, commercial) %>%
  summarize(Sp_mean = mean(prop_w_plastic),
            Sample_size = sum(N),
            num_studies = n_distinct(source)) %>%
  ungroup %>% 
  mutate(commercial = factor(commercial),
         studies_cat = cut(num_studies, 
                           c(0, 1, 3, 6),
                           c("Poorly studied (n=1)", "Moderately studied (n=2-3)", "Well studied (n=4-6)")),
         studies_cat_num = as.double(cut(num_studies, 
                           c(0, 1, 3, 6),
                           c(1,2,3))),
         commercial = fct_collapse(commercial,
                                   Commercial = c("commercial", "highly commercial"),
                                   Minor = c("minor commercial", "subsistence"),
                                   None = "none")) %>% 
  arrange(num_studies) %>% 
  drop_na(commercial) %>% 
  ggplot(aes(log10(Sample_size), Sp_mean)) +
  geom_point(aes(color = Sp_mean, shape = commercial, size = studies_cat_num), 
             alpha = 0.8) +
  geom_hline(yintercept = 0.25, linetype="dashed", color = "grey50") +
  geom_vline(xintercept = log10(25), linetype="dashed", color = "grey50") +
  #facet_wrap(~ commercial) +
  xlab("Log[Sample Size]") +
  ylab("Species-specific plastic ingestion incidence (FO)") +
  labs(shape = "Commercial status", shape = "Commercial Status", size =  "Number of studies") +
  #scale_color_manual(values = c("red","blue", "black")) +
  scale_color_gradientn(colours = c("steelblue4",
                                    "gray40",
                                    "coral", "coral1",
                                    "firebrick2", "firebrick3", "firebrick4"), 
                        name = "Proportion with \ningested plastic") +
  scale_size_continuous(breaks = seq(from = 1, to = 3, by = 1), 
                        labels = c("Poorly studied (n=1)", "Moderately studied (n=2-3)", "Well studied (n=4-6)"),
                        range = c(1.5, 5)) +
  annotate("text", x = c(0.6, 3.8, 0.6, 3.8),
           y=c(0.9, 0.9, 0.06, 0.055),
           label = c("high incidence, data poor", "high incidence, data rich",
                     "low incidence, data poor", "low incidence, data rich")) +
  theme_classic(base_size = 14)
risk_plot

dev.copy2pdf(file="risk_plot.pdf", width=11, height=7)

binom_plot <- ggplot(filter(d, includes_microplastic != "NA"), 
                     aes(publication_year, includes_microplastic)) +
  geom_point() +
  geom_smooth(method = "glm", method.args = list(family = "binomial"))
binom_plot

p <- ggplot(d_full, 
           aes(trophic_level_via_fishbase, prop_w_plastic, size= N, weight = N)) + 
  geom_point(alpha = 0.4) +  # Eventually add in foraging behavior here 
  geom_smooth(col = "blue", method = "loess") +
  xlim(2, 5) +
  xlab("Trophic level") +
  ylab("Proportion of individuals with ingested plastic") +
  annotate("text", x = c(2.75, 4), y= -0.05, 
           label = c("planktivorous", "piscivorous")) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold")) 
p + guides(size = FALSE)



p3_b <- ggplot(filter(d, Found != "NA"), 
               aes(average_depth, prop_w_plastic, size = N, weight = N)) +
  geom_point(alpha = 0.4) + 
  geom_smooth(col = "blue", method = "loess", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ Found, scales = "free_y", ncol = 1) +
  #ylim(0,1) + 
  coord_flip() +
  scale_x_reverse() +
  xlim(1500,0) +
  xlab("Average depth") +
  ylab("Proportion of individuals with ingested plastic") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold"),
        strip.text = element_text(size = 10)) 
p3_b + guides(size = FALSE)


p_raincloud <- ggplot(d,aes(x=Found,y=prop_w_plastic, fill = Found))+
  geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust = 2)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  
  ylab('Score')+xlab('Group')+coord_flip()+theme_classic()+guides(fill = FALSE)+
  ggtitle('Figure R3: The Basic Raincloud with Colour')
p_raincloud 


p4 <- ggplot(filter(d1, `Oceanographic province (from Longhurst 2007)` != "NA"), 
             aes(`Oceanographic province (from Longhurst 2007)`, `Prop w plastic`)) +
  geom_jitter() + 
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p4


p5 <- ggplot(d, aes(order, mean_num_particles_per_indv)) +
  #geom_jitter() + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p5


p6 <- ggplot(d, 
             aes(Vulnerability.score..via.fishbase.from.Cheug.et.al.2005., Prop.w.plastic,
                 weight = N)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title=element_text(size=13, face="bold"))
p6


p7 <- ggplot(filter(d, prime_forage != "NA"), 
                aes(prime_forage, prop_w_plastic)) +
  geom_jitter(aes(size = N)) + 
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p7


p8 <- ggplot(filter(d1, `Commercial` != "NA"), 
             aes(`Commercial`, `Prop w plastic`)) +
  geom_jitter() + 
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p8


p9 <- ggplot(filter(d, IUCN.status != "NA"), 
             aes(IUCN.status, Prop.w.plastic)) +
  geom_jitter() + 
  geom_boxplot(alpha = 0.1, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p9

<<<<<<< HEAD
# GLMM ----
# USING THIS ONE WITH SCALED DEPTH, certified by Gemma

d_full_wo_gaps <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(average_depth, Found, trophic_level_via_fishbase, prime_forage) 

glmm_FwP <- glmer(cbind(NwP, N-NwP) ~ scale(average_depth)*Found + trophic_level_via_fishbase + prime_forage +
                  (1|source) + (1|order), 
                  data = d_full_wo_gaps, family = binomial, na.action = "na.fail")
summary(glmm_FwP)
confint(ranef(glmm_FwP), method="Wald")
ranef(glmm_FwP)



# multi-model selection using AICc
GLMM_dredge <- dredge(glmm_FwP)

View(GLMM_dredge)
write_csv(GLMM_dredge, "GLMM model selection table.csv")
#subset(GLMM_dredge, delta < 4)
a=model.avg(GLMM_dredge)
summary(a) #The ‘subset’ (or ‘conditional’) average only averages over the models where the parameter appears. An alternative, the ‘full’ average assumes that a variable is included in every model
confint(a) #computes confidence interval



glmm_Commerical <- glmer(cbind(NwP, N-NwP) ~ commercial +
                    (1|source) + (1|order), data = d_full_wo_gaps, family = binomial)
summary(glmm_Commerical)




d_num_plast <- d_full %>%
  select(mean_num_particles_per_indv,
         average_depth, Found, trophic_level_via_fishbase, prime_forage) %>% 
  drop_na(mean_num_particles_per_indv,
          average_depth, Found, trophic_level_via_fishbase, prime_forage) 


%>% 
  filter(mean_num_particles_per_indv > 0)

glmm_num_plast <- glmer(mean_num_particles_per_indv ~ scale(average_depth)*Found + trophic_level_via_fishbase + prime_forage +
                          (1|source) + (1|order), 
                        data = d_num_plast, na.action = "na.fail")
summary(glmm_num_plast)
confint(ranef(glmm_FwP), method="Wald")
ranef(glmm_FwP)



# trying a gamm ----


# wrapper fucntion  needed to allow MuMIn to compare models
gamm4 <- function(...) structure(c(gamm4::gamm4(...), list(call = match.call())), class = c("gamm", "list"))  

d_full_wo_gaps <- d %>%
  filter(includes_microplastic == "Y") %>% 
  drop_na(average_depth, Found, trophic_level_via_fishbase, prime_forage) 

# final models
gamm4_FwP_w_conditional <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + s(scale(average_depth), by = Found) + 
                                   Found + prime_forage, 
                                 random = ~(1|order) + (1|source), 
                                 data = d_full_wo_gaps, family = binomial)

summary(gamm4_FwP_w_conditional$gam)
plot(gamm4_FwP_w_conditional$gam)

# gam_FwP_w_conditional <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + s(scale(average_depth), by = Found) + 
#                                    Found + prime_forage,
#                                  data = d_glmm_full, family = binomial)

gamm4_FwP <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + s(scale(average_depth)) + Found + prime_forage, 
                   random = ~(1|order) + (1|source), 
                  data = d_full_wo_gaps, family = binomial)
summary(gamm4_FwP$gam)
plot(gamm4_FwP$gam)


# gamm4_full <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + s(scale(average_depth)) + Found + prime_forage, 
#                    random = ~(1|order) + (1|source), 
#                    data = d_glmm_full, family = binomial)
# summary(gamm4_FwP$gam)
# plot(gamm4_FwP$gam)

gamm4_FwP_wo_TL <- gamm4(cbind(NwP, N-NwP) ~  + s(scale(average_depth)) + Found + prime_forage, 
                         random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
summary(gamm4_FwP_wo_TL$gam)
plot(gamm4_FwP_wo_TL$gam)

gamm4_FwP_wo_TL_w_interaction <- gamm4(cbind(NwP, N-NwP) ~  + s(scale(average_depth), by = Found) + prime_forage, 
                         random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_wo_depth <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + Found + prime_forage, 
                            random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_wo_depth_w_interaction <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + Found + prime_forage, 
                            random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_wo_Found <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + s(scale(average_depth)) + prime_forage, 
                            random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)

gamm4_FwP_wo_Forage <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + s(scale(average_depth)) + Found, 
                             random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_wo_Forage_w_interaction <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase) + s(scale(average_depth), by = Found) + 
                                             Found, 
                             random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_TLonly <- gamm4(cbind(NwP, N-NwP) ~ s(trophic_level_via_fishbase), 
                          random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_depth_cond_only <- gamm4(cbind(NwP, N-NwP) ~ + s(scale(average_depth), by = Found), 
                                   random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_depth_only <- gamm4(cbind(NwP, N-NwP) ~ + s(scale(average_depth)), 
                              random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_Foundonly <- gamm4(cbind(NwP, N-NwP) ~ Found, 
                          random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_FwP_Forageonly <- gamm4(cbind(NwP, N-NwP) ~ prime_forage, 
                              random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)
gamm4_null <- gamm4(cbind(NwP, N-NwP) ~ 1, 
                    random = ~(1|order) + (1|source), data = d_full_wo_gaps, family = binomial)

# compare models
GAMM_mod.sel <- model.sel(gamm4_FwP_w_conditional, gamm4_FwP,
                       gamm4_FwP_wo_TL, gamm4_FwP_wo_TL_w_interaction,
                       gamm4_FwP_wo_depth, gamm4_FwP_wo_depth_w_interaction,
                       gamm4_FwP_wo_Found, 
                       gamm4_FwP_wo_Forage, gamm4_FwP_wo_Forage_w_interaction,
                       gamm4_FwP_TLonly, gamm4_FwP_depth_only, gamm4_FwP_Foundonly, gamm4_FwP_Forageonly,
                       gamm4_null)
View(GAMM_mod.sel)

GAMM_mod.avg=model.avg(GAMM_mod.sel)
summary(GAMM_mod.avg)


# best model
summary(gamm4_FwP_wo_TL$mer)
summary(gamm4_FwP_wo_depth$gam)
plot(gamm4_FwP_wo_depth$gam)


# mgcv gamms
gamm_FwP_wconditional <- gam(prop_w_plastic ~ s(trophic_level_via_fishbase) + s(average_depth, by = Found) + Found + prime_forage + s(order, bs="re"),
                      data = d, family = gaussian)
gam.check(gamm_FwP_wconditional)

plot_smooth(gamm_FwP_wconditional, view="average_depth", cond=list(Group="Found"))

gamm_num_plast <- gam(mean_num_particles_per_indv ~ s(trophic_level_via_fishbase) + s(average_depth) + Found,
                             data = d_num_plast, family = gaussian)
gam.check(gamm_num_plast)
summary(gamm_num_plast)
plot(gamm_num_plast)

mgcv.models <- model.sel(gamm_FwP_wconditional, gamm_FwP)
foo2=model.avg(mgcv.models)


# gamm_FwP <- gamm(Prop.w.plastic ~ s(trophic_level_via_fishbase,k=10)+s(Average.depth, by = Found, k=10) + Found, data=d)
# ### $gam to look at gam effects. $lme to look at random effects.
# summary(gamm_FwP$gam)
# plot(gamm_FwP$gam)

# some other analyses


# MCMC GLMM ----
MCMCglmm_FwP <- MCMCglmm(mean_num_particles_per_indv ~ scale(average-depth)*Found + trophic_level_via_fishbase + prime_forage,
                         random = ~ order,
                         data = d_full_wo_gaps, family = NULL,
                         nitt = 1150, thin = 100, burnin = 150,
                         pr = TRUE, # To save the posterior mode of for each level or category in your random effect(s) we use pr = TRUE, which saves them in the $Sol part of the model output.
                         verbose = TRUE)

d_full_wo_gaps <- select(d_full_wo_gaps, -family)

MCMCglmm_FwP <- MCMCglmm(prop_w_plastic ~ scale(average_depth) + Found + prime_forage,
                         random = ~(1|order) + (1|source), 
                         data = d_full_wo_gaps, family = NULL,
                         nitt = 1150, thin = 100, burnin = 150,
                         pr = TRUE, # To save the posterior mode of for each level or category in your random effect(s) we use pr = TRUE, which saves them in the $Sol part of the model output.
                         verbose = TRUE)

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
  
# Junk code below here  
#   
#   # proportion ingestion plastic with depth (perhaps remove demersal here)
#   p3 <- ggplot(filter(d1, Habitat %in% c("demersal", "pelagic-neritic", "pelagic-oceanic", "reef-associated", "benthopelagic")),
#                aes(-`Average depth`, `Prop w plastic`, size=N, weight = N)) +  #col=Family
#   geom_point(alpha = 0.4, aes(color = Basin)) + 
#   geom_smooth(col = "grey20", method = "loess") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   facet_wrap(~ Habitat, scales = "free_y", ncol = 1) +
#   #ylim(0,1) + 
#   #xlim(-500,0) +
#   coord_flip() +
#   theme_classic()
# p3  
  
  
