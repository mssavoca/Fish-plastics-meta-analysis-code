####################################
# Phylogenetic tree figure for paper
####################################


# load data and packages ----

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



# Abbreviate a binomial e.g. Balaenoptera musculus -> B. musculus
abbr_binom <- function(binom) {
  paste(str_sub(binom, 1, 1), 
        str_extract(binom, " .*"), 
        sep = ".")
}


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


# species summary tables
d_sp_sum <- d_full %>%
  filter(!species %in% c("sp.", "spp.")) %>%
  group_by(binomial, order) %>%
  summarize(Sp_mean = mean(prop_w_plastic),
            Sample_size = sum(N),
            num_studies = n_distinct(source)) %>%
  arrange(-Sp_mean)



# building the basic tree ----

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

View(my_tree)

#my_tree2 = phylo4d(my_tree, d_sp_sum)


# first plot try, fan layout ----
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
               align = TRUE, size = 6) +
  # geom_tiplab2(aes(color=Sp_mean, angle = angle),
  #              align = TRUE, size = 6) +
  scale_color_gradientn(colours = c("steelblue4","steelblue3",
                                    "coral", "coral1",
                                    "firebrick2", "firebrick3", "firebrick4"), 
                       name = "Proportion with \ningested plastic") +
  theme(legend.position = c(0.5,0.5),
        legend.key.size = unit(1.75, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20))



dev.copy2pdf(file="Prelim_phylo_ggtree.pdf", width=30, height=30)







# other options

# first plot try, vertical layout with data ----
p_stand <- ggtree(my_tree) + 
  #geom_text2(aes(label=label), hjust=-.2, size=4) +
  ggplot2::xlim(0, 1.2)
p_stand

# Adding data
p_stand %<+% d_sp_sum +
  aes(color = Sp_mean) +
  geom_tiplab(aes(color=Sp_mean), align = TRUE, size = 3) +
  scale_color_gradientn(colours = c("#053061" ,"#2166AC", "#4393C3", "#F4A582", "#D6604D", "#B2182B", "#67001F"), 
                         name = "Proportion with \nplastic") +
  geom_hilight(node=1, fill="gold") + 
  theme(legend.position = c(0.2,0.8))


p2 <- facet_plot(p_stand, panel = 'Prop w plastic', data = d_sp_sum, 
                 geom = ggstance::geom_barh, 
                 mapping = aes(x = Sp_mean), 
                 stat='identity') 
p2

geom_facet(

inset(my_tree, d_sp_sum$Sp_mean, width=0.2, height=0.15, hjust=-1)
        
        
         legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys
  
  
  
  geom_tippoint(aes(color=order))



  geom_tiplab(aes(color = Sp_mean))


p %<+% dd + geom_text(aes(color=, label=label), hjust=-0.5)


