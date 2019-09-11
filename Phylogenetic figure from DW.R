###########################################################################
# Phylogenetic tree figure code for fish-plastic meta-analysis from Danuta
###########################################################################


library(dplyr)
library(ggplot2)
library(ggtree)
library(gbm)
library(dismo)
library(mgcv)
library(lme4)
library(gamm4)
library(readxl)

d = read_csv("Plastics ingestion records fish master_final.csv") %>% 
  janitor::clean_names() %>% 
  rename(NwP = nw_p,
         N = n) %>% 
  mutate(WeightedProp = prop_w_plastic * N,
         Found = as_factor(case_when(habitat %in% c("demersal", "reef-associated", "benthopelagic", "bathydemersal") ~ "demersal",
                                     habitat %in% c("pelagic-neritic", "pelagic-oceanic", "mesopelagic", "bathypelagic") ~ "pelagic")),
         prime_forage = na_if(prime_forage, "not listed")) %>% 
  separate(binomial, into = c("genus", "species"), sep = " ", remove = FALSE)

# some summary tables
d_sp_sum <- d %>%
  filter(!species %in% c("sp.", "spp.")) %>%
  group_by(binomial, order) %>%
  summarize(Sp_mean = mean(prop_w_plastic),
            Sample_size = sum(N),
            num_studies = n_distinct(source)) %>%
  arrange(-Sp_mean)


# testing phylogenetic analyses ----
library(rotl) #for phylogenetic analyses, get all the species? from Hinchliff et al. 2015 PNAS
library(phytools)
library(tidyverse)
library(stringr)


#gets the species names
# taxa <- as.character(d_sp_sum$`Species.name`[1:20])  # this is just a subset, will eventually include all species
# taxa <- as.character(d_sp_sum$`Species.name`[1:100])  # this is just a subset, will eventually include all species
# 
# resolved_names <- tnrs_match_names(taxa)
# resolved_names <- resolved_names[!is.na(resolved_names$unique_name),] #removing NAs when no match was found (Genus sp. I think)
breaks <- c(seq(1,nrow(d_sp_sum),50),nrow(d_sp_sum)+1)

##### MAX'S ATTEMPT #####
for (i in 1:(length(breaks)-1)){
  # Subset binomials, convert to character, drop missing
  taxa <- as.character(d_sp_sum$binomial[breaks[i]:(breaks[i+1]-1)])
  taxa <- taxa[taxa!=""]
  # Use rotl to look up ottid's
  resolved_namest <- tnrs_match_names(taxa, "Animals")
  resolved_namest <- resolved_namest[!is.na(resolved_namest$unique_name),]
  # Drop ottid's missing from synthetic tree (e.g. incertae sedis)
  resolved_namest <- resolved_namest[is_in_tree(resolved_namest$ott_id),]
  if (i==1){
    resolved_namess <- resolved_namest
  } else {
    resolved_namess <- full_join(resolved_namess,resolved_namest)
  }
}
my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id, label_format = "name")
#########################

for (i in 1:(length(breaks)-1)){
  taxa <- as.character(d_sp_sum$binomial[breaks[i]:(breaks[i+1]-1)])
  taxa <- taxa[taxa!=""]
  resolved_namest <- tnrs_match_names(taxa)
  resolved_namest <- resolved_namest[!is.na(resolved_namest$unique_name),]
  if (i==1){
    resolved_namess <- resolved_namest
  } else {
    resolved_namess <- full_join(resolved_namess,resolved_namest)
  }
}
resolved_names <- resolved_namess
resolved_names <- resolved_names[resolved_names$flags!="INCERTAE_SEDIS_INHERITED",]

#plots species
my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id, label_format = "name")
# plot(my_tree, no.margin=TRUE)

# tree<-read.tree(my_tree)# not working, but doesnt seem to matter atm

plot.tree <- compute.brlen(my_tree, method = "Grafen", power = 1/2) #add branch lengths to my tree using the Grafen (1989) method
plot.tree <- ladderize(plot.tree, right = TRUE)

# This read/write stage is essential to getting the bars to line up with the tree & labels
write.tree(plot.tree, "Plotting - tree.phy") 
plot.tree <- read.tree("Plotting - tree.phy")

# # Mimicing viridis colour palette for colour blindness
# # https://www.thinkingondata.com/something-about-viridis-library/
# col.PROCELLARIIFORMES <- "#440154FF" # purple 481567FF
# col.SPHENISCIFORMES <- "#453781FF" # purple blue
# col.SULIFORMES <- "#39568CFF" # blue
# col.PELECANIFORMES <- "#2D708EFF" # dark teal
# col.GAVIIFORMES <- "#1F968BFF" # teal
# col.PHAETHONTIFORMES <- "#3CBB75FF" # grean-teal
# col.CHARADRIIFORMES <- "#73D055FF" # green
# col.ANSERIFORMES <- "#DCE319FF" # yellow  copy this to line 557 otherwise one branch is pink clcolr[all]  

# # Mimicing viridis colour palette for colour blindness
# # https://www.thinkingondata.com/something-about-viridis-library/
# col.Clupeiformes <- "#440154FF" # purple 481567FF
# col.Perciformes <- "#453781FF" # purple blue
# col.Pleuronectiformes <- "#39568CFF" # blue
# col.Myctophiformes <- "#2D708EFF" # dark teal
# col.Stomiiformes <- "#1F968BFF" # teal
# col.Aulopiformes <- "#3CBB75FF" # grean-teal
# col.Beloniformes  <- "#73D055FF" # green
# col.Mugiliformes <- "#DCE319FF" # yellow  copy this to line 557 otherwise one branch is pink clcolr[all]  


#Data to plot
plot.data <- d_sp_sum
plot.data <- plot.data[plot.data$binomial!="",]
plot.data$TipLabel <- NA
for(i in 1:nrow(plot.data)){
  oj <- str_locate(plot.data$binomial[i],' ')
  plot.data$TipLabel[i] <- paste(substr(plot.data$binomial[i],1,oj[1]-1),substr(plot.data$binomial[i],oj[1]+1,nchar(as.character(plot.data$binomial[i]))),sep="_") 
}
plot.data$TipLabel <- as.factor(plot.data$TipLabel)
# rownames(plot.data) <- as.character(plot.data$TipLabel)

#Species missing from the tree
plot.data <- droplevels(plot.data)
species.list <- as.character(levels(plot.data$TipLabel))
species.list[!species.list %in% plot.tree$tip.label]

#Remove tips for species that are not in the data set for plotting
plot.tree <- drop.tip(plot.tree, tip = plot.tree$tip.label[!plot.tree$tip.label %in% plot.data$TipLabel])

plot.data.df <- as.data.frame(plot.data)
row.names(plot.data.df) <- make.names(as.character(plot.data.df$TipLabel),unique=TRUE)
plot.data.df <- plot.data.df[plot.tree$tip.label,]

ord <- unique(plot.data.df$order)
library(scales)
q_colors =  length(ord) # for no particular reason
v_colors =  viridisLite::viridis(q_colors,option="E")   #sets colors according to cividis scale



# 
# plot.data$col[plot.data$Order == "ANSERIFORMES"] <- col.ANSERIFORMES
# plot.data$col[plot.data$Order == "CHARADRIIFORMES"] <- col.CHARADRIIFORMES
# plot.data$col[plot.data$Order == "GAVIIFORMES"] <- col.GAVIIFORMES
# plot.data$col[plot.data$Order == "PHAETHONTIFORMES"] <- col.PHAETHONTIFORMES
# plot.data$col[plot.data$Order == "PROCELLARIIFORMES"] <- col.PROCELLARIIFORMES
# plot.data$col[plot.data$Order == "SPHENISCIFORMES"] <- col.SPHENISCIFORMES
# plot.data$col[plot.data$Order == "SULIFORMES"] <- col.SULIFORMES
# plot.data$col[plot.data$Order == "PELECANIFORMES"] <- col.PELECANIFORMES

head(plot.tree$tip.label)
head(plot.data.df)


tree.angle <- 330
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


# #Get the colours
# ANSERIFORMES_edge      <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "ANSERIFORMES")))
# CHARADRIIFORMES_edge   <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "CHARADRIIFORMES")))
# GAVIIFORMES_edge       <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "GAVIIFORMES")))
# PHAETHONTIFORMES_edge  <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "PHAETHONTIFORMES")))
# PROCELLARIIFORMES_edge <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "PROCELLARIIFORMES")))
# SPHENISCIFORMES_edge   <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "SPHENISCIFORMES")))
# SULIFORMES_edge        <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "SULIFORMES")))
# PELECANIFORMES_edge    <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "PELECANIFORMES")))
# all                    <- which.edge(plot.tree, group = row.names(plot.data))
# 
# ANSERIFORMES_edge      <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "ANSERIFORMES")))
# CHARADRIIFORMES_edge   <- which.edge(plot.tree, group = row.names(subset(plot.data, Order != "ANSERIFORMES")))
# PHAETHONTIFORMES_edge  <- which.edge(plot.tree, group = row.names(subset(plot.data, !(Order == "CHARADRIIFORMES" | Order == "ANSERIFORMES"))))
# GAVIIFORMES_edge       <- which.edge(plot.tree, group = row.names(subset(plot.data, !(Order == "ANSERIFORMES" | Order == "CHARADRIIFORMES" |
#                                                                                         Order == "PHAETHONTIFORMES"))))
# PELECANIFORMES_edge    <- which.edge(plot.tree, group = row.names(subset(plot.data, !(Order == "ANSERIFORMES" | Order == "CHARADRIIFORMES" |
#                                                                                         Order == "PHAETHONTIFORMES" | Order == "GAVIIFORMES"))))
# SULIFORMES_edge        <- which.edge(plot.tree, group = row.names(subset(plot.data, !(Order == "ANSERIFORMES" | Order == "CHARADRIIFORMES" |
#                                                                                         Order == "PHAETHONTIFORMES" | Order == "GAVIIFORMES" |
#                                                                                         Order == "PELECANIFORMES"))))
# SPHENISCIFORMES_edge   <- which.edge(plot.tree, group = row.names(subset(plot.data, !(Order == "ANSERIFORMES" | Order == "CHARADRIIFORMES" |
#                                                                                         Order == "PHAETHONTIFORMES" | Order == "GAVIIFORMES" |
#                                                                                         Order == "PELECANIFORMES" | Order == "SULIFORMES"))))
# PROCELLARIIFORMES_edge <- which.edge(plot.tree, group = row.names(subset(plot.data, Order == "PROCELLARIIFORMES")))


plot.data.df$col <- NA
for (i in 1:length(ord)){
  assign(paste("col.",ord[i],sep=""),v_colors[i])
  plot.data.df$col[plot.data.df$order == ord[i]] <- v_colors[i]
  assign(paste(ord[i],"edge",sep="_"),which.edge(plot.tree, group = row.names(subset(plot.data.df, order == ord[i]))))
}
all <- which.edge(plot.tree, group = row.names(plot.data.df))
col1 <- plot.data.df$col

# tree.plot.colours <- data.frame(c(col.ANSERIFORMES, col.CHARADRIIFORMES, col.GAVIIFORMES,
#                                   col.PHAETHONTIFORMES, col.PROCELLARIIFORMES, col.SPHENISCIFORMES, col.SULIFORMES, col.PELECANIFORMES))
# 
# row.names(tree.plot.colours) <- c("ANSERIFORMES", "CHARADRIIFORMES", "GAVIIFORMES", 
#                                   "PHAETHONTIFORMES", "PROCELLARIIFORMES", "SPHENISCIFORMES",
#                                   "SULIFORMES", "PELECANIFORMES")

tree.plot.colours <- data.frame(c(col.Perciformes, col.Tetraodontiformes, col.Lophiiformes, col.Scorpaeniformes,
                                  col.Gasterosteiformes, col.Pleuronectiformes, col.Beloniformes, col.Atheriniformes, col.Mugiliformes,
                                  col.Ophidiiformes, col.Beryciformes, col.Lampriformes, col.Gadiformes, col.Zeiformes, col.Myctophiformes,
                                  col.Aulopiformes, col.Stomiiformes, col.Osmeriformes, col.Salmoniformes, col.Clupeiformes, col.Siluriformes,
                                  col.Anguilliformes, col.Carcharhiniformes, col.Squaliformes, col.Rajiformes, col.Myliobatiformes, col.Torpediniformes))

row.names(tree.plot.colours) <- c("Perciformes", "Tetraodontiformes", "Lophiiformes", "Scorpaeniformes",
                                  "Gasterosteiformes", "Pleuronectiformes", "Beloniformes", "Atheriniformes", "Mugiliformes",
                                  "Ophidiiformes", "Beryciformes", "Lampriformes", "Gadiformes", "Zeiformes", "Myctophiformes",
                                  "Aulopiformes", "Stomiiformes", "Osmeriformes", "Salmoniformes", "Clupeiformes", "Siluriformes",
                                  "Anguilliformes", "Carcharhiniformes", "Squaliformes", "Rajiformes", "Myliobatiformes", "Torpediniformes")


#Colours for the tree
clcolr<-rep(as.character(tree.plot.colours[ord[1],]), dim(plot.tree$edge)[1]) #vector as long as the first edge column and populating it with the colour for insects
clcolr[all]    <- "#DE1738" #as.character(tree.plot.colours["ANSERIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Perciformes_edge]      <-as.character(tree.plot.colours["Perciformes",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Tetraodontiformes_edge]   <-as.character(tree.plot.colours["Tetraodontiformes",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Lophiiformes_edge]  <-as.character(tree.plot.colours["Lophiiformes",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Scorpaeniformes_edge]       <-as.character(tree.plot.colours["Scorpaeniformes",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Gasterosteiformes_edge]    <-as.character(tree.plot.colours["Gasterosteiformes",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Pleuronectiformes_edge]        <-as.character(tree.plot.colours["Pleuronectiformes",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Beloniformes_edge]   <-as.character(tree.plot.colours["Beloniformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Atheriniformes_edge] <-as.character(tree.plot.colours["Atheriniformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Mugiliformes_edge]   <-as.character(tree.plot.colours["Mugiliformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Ophidiiformes_edge] <-as.character(tree.plot.colours["Ophidiiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Beryciformes_edge]   <-as.character(tree.plot.colours["Beryciformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Lampriformes_edge] <-as.character(tree.plot.colours["Lampriformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Gadiformes_edge]   <-as.character(tree.plot.colours["Gadiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Zeiformes_edge] <-as.character(tree.plot.colours["Zeiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Myctophiformes_edge]   <-as.character(tree.plot.colours["Myctophiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Aulopiformes_edge] <-as.character(tree.plot.colours["Aulopiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Stomiiformes_edge]   <-as.character(tree.plot.colours["Stomiiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Osmeriformes_edge] <-as.character(tree.plot.colours["Osmeriformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Salmoniformes_edge]   <-as.character(tree.plot.colours["Salmoniformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Clupeiformes_edge] <-as.character(tree.plot.colours["Clupeiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Siluriformes_edge]   <-as.character(tree.plot.colours["Siluriformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Anguilliformes_edge] <-as.character(tree.plot.colours["Anguilliformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Carcharhiniformes_edge]   <-as.character(tree.plot.colours["Carcharhiniformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Squaliformes_edge] <-as.character(tree.plot.colours["Squaliformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Rajiformes_edge] <-as.character(tree.plot.colours["Rajiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Myliobatiformes_edge] <-as.character(tree.plot.colours["Myliobatiformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
clcolr[Torpediniformes_edge] <-as.character(tree.plot.colours["Torpediniformes",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not

# #Colours for the tree
# clcolr<-rep(as.character(tree.plot.colours["ANSERIFORMES",]), dim(plot.tree$edge)[1]) #vector as long as the first edge column and populating it with the colour for insects
# clcolr[all]    <- "#DCE319FF" #as.character(tree.plot.colours["ANSERIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[ANSERIFORMES_edge]      <-as.character(tree.plot.colours["ANSERIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[CHARADRIIFORMES_edge]   <-as.character(tree.plot.colours["CHARADRIIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[PHAETHONTIFORMES_edge]  <-as.character(tree.plot.colours["PHAETHONTIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[GAVIIFORMES_edge]       <-as.character(tree.plot.colours["GAVIIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[PELECANIFORMES_edge]    <-as.character(tree.plot.colours["PELECANIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[SULIFORMES_edge]        <-as.character(tree.plot.colours["SULIFORMES",])   # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[SPHENISCIFORMES_edge]   <-as.character(tree.plot.colours["SPHENISCIFORMES",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# clcolr[PROCELLARIIFORMES_edge] <-as.character(tree.plot.colours["PROCELLARIIFORMES",]) # this defines the ingroup and then sets up a vector to colour each edge of the phylogeny depending on whether it is in the group or not
# 

par(mar = c(3,3,3,3) + 5, xpd = NA)

# (phylo.plot is from ape package)
plot(plot.tree, type = "fan", open.angle = 360 - tree.angle, rotate = 270, 
     root.edge = TRUE, 
     show.tip.label = T, label.offset = .36, cex = 0.6, tip.color = col1,
     edge.width = 1.5, # font = 3, 
     x.lim = c(-1 * treeheight, 1.2 * treeheight), 
     y.lim = c(-1 * treeheight, 1.2 * treeheight), 
     edge.color = clcolr)
# text(x = -1.55, # par("usr")[1] + 0.99* (par("usr")[2] - par("usr")[1]),
#      y = 1.6,   # par("usr")[3] + 0.99* (par("usr")[4] - par("usr")[3]),
#      labels = "(J)")



### Draw Prediction bars, background and labels

# Draw background
theta.start <- -0.5*(pi/180)*tree.angle/nrow(plot.data.df) + (pi/180)*tree.start
theta.subtend <- abs(2*theta.start)

bar.scale <- 0.3 #Multiplication factor for the bars
bar.offset <- 0.03 #Offset to add white space between the bars and the tree, to separate the colours

x.bg       <- sin(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset)
y.bg       <- cos(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset)
x.bg.outer <- sin(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (bar.scale*treeheight))
y.bg.outer <- cos(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (bar.scale*treeheight))
x.0.25     <- sin(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (0.25*bar.scale*treeheight))
y.0.25     <- cos(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (0.25*bar.scale*treeheight))
x.0.5      <- sin(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (0.5*bar.scale*treeheight))
y.0.5      <- cos(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (0.5*bar.scale*treeheight))
x.0.75     <- sin(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (0.75*bar.scale*treeheight))
y.0.75     <- cos(theta.start+c(-1:(nrow(plot.data.df)-1))*theta.subtend)*(treeheight + bar.offset + (0.75*bar.scale*treeheight))

#par(mfrow = c(1,1), mar = c(4,4,0,0)+0.5)
#plot(y.bg.outer ~ x.bg.outer, type = "n")
polygon(x = c(x.bg, rev(x.bg.outer)),
        y = c(y.bg, rev(y.bg.outer)),
        col = "grey93", border = "grey93") # This could also be "dark grey" 

lines(y.0.25 ~ x.0.25, col = "black", lwd = 1.5) # This could also be "white"
lines(y.0.5 ~ x.0.5, col = "black", lwd = 1.5)
lines(y.0.75 ~ x.0.75, col = "black", lwd = 1.5)

theta.start+c(1:nrow(plot.data))*theta.subtend


# Bars    # NOT WORKING
temp.x <- NA
temp.y <- NA
temp.theta <- NA
for(i in 1:nrow(plot.data.df)){
  #for(i in 1:4){ 
  #outer layer first
  theta1 <- theta.start + (i-1)*theta.subtend
  theta2 <- theta1 - theta.subtend
  
  x1 <- sin(theta1)*(treeheight + bar.offset)
  x2 <- sin(theta2)*(treeheight + bar.offset)
  y1 <- cos(theta1)*(treeheight + bar.offset)
  y2 <- cos(theta2)*(treeheight + bar.offset)
  
  temp.theta[i] <- theta1
  temp.x[i] <- x1
  temp.y[i] <- y1
  
  x3 <- sin(theta1)*(treeheight + bar.offset + (bar.scale*treeheight)*(plot.data.df[i,"Sp_mean"]))
  x4 <- sin(theta2)*(treeheight + bar.offset + (bar.scale*treeheight)*(plot.data.df[i,"Sp_mean"]))
  y3 <- cos(theta1)*(treeheight + bar.offset + (bar.scale*treeheight)*(plot.data.df[i,"Sp_mean"]))
  y4 <- cos(theta2)*(treeheight + bar.offset + (bar.scale*treeheight)*(plot.data.df[i,"Sp_mean"]))
  
  x5 <- sin(theta1)*(treeheight + bar.offset + (bar.scale*treeheight))
  x6 <- sin(theta2)*(treeheight + bar.offset + (bar.scale*treeheight))
  y5 <- cos(theta1)*(treeheight + bar.offset + (bar.scale*treeheight))
  y6 <- cos(theta2)*(treeheight + bar.offset + (bar.scale*treeheight))
  #x3.1 <- sin(theta1)*(treeheight + bar.offset + (bar.scale*treeheight)*(data[i,1]-min(data[,1]))/(max(data[,1])-min(data[,1])) + 
  #                       (bar.scale*treeheight)*(data[i,2]-min(data[,2]))/(max(data[,2])-min(data[,2])))
  #x4.1 <- sin(theta2)*(treeheight + bar.offset + (bar.scale*treeheight)*(data[i,1]-min(data[,1]))/(max(data[,1])-min(data[,1])) + 
  #                       (bar.scale*treeheight)*(data[i,2]-min(data[,2]))/(max(data[,2])-min(data[,2])))
  #y3.1 <- cos(theta1)*(treeheight + bar.offset + (bar.scale*treeheight)*(data[i,1]-min(data[,1]))/(max(data[,1])-min(data[,1])) + 
  #                       (bar.scale*treeheight)*(data[i,2]-min(data[,2]))/(max(data[,2])-min(data[,2])))
  #y4.1 <- cos(theta2)*(treeheight + bar.offset + (bar.scale*treeheight)*(data[i,1]-min(data[,1]))/(max(data[,1])-min(data[,1])) + 
  #                       (bar.scale*treeheight)*(data[i,2]-min(data[,2]))/(max(data[,2])-min(data[,2])))
  
  # polygon(c(x1,x2,x6,x5),
  #         c(y1,y2,y6,y5),
  #         #col = col1[i], border = col1[i],
  #         col = "grey80", border = "grey80")
  
  polygon(c(x1,x2,x4,x3),
          c(y1,y2,y4,y3),
          #col = "darkorange", border = "darkorange",
          col = col1[i], border = col1[i])
  #polygon(c(x3,x4,x4.1,x3.1),c(y3,y4,y4.1,y3.1),
  #        #col = "darkgreen", border = "darkgreen",
  #        col = col2[i], border = col2[i])
}


### Add scales


# Upper
lines(y = rep(-0*bar.offset,2), 
      x = c(-1*(treeheight+bar.offset), 
            -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(y = c(-0*bar.offset,-1*bar.offset), 
      x = c(-1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(y = c(-0*bar.offset,-1*bar.offset), 
      x = c(-1*(treeheight+bar.offset+(0.5*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(0.5*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(y = c(-0*bar.offset,-1*bar.offset), 
      x = c(-1*(treeheight+bar.offset+(0.25*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(0.25*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(y = c(-0*bar.offset,-1*bar.offset), 
      x = c(-1*(treeheight+bar.offset+(0.75*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(0.75*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(y = c(-0*bar.offset,-1*bar.offset), 
      x = c(-1*(treeheight+bar.offset), 
            -1*(treeheight+bar.offset)),
      lwd = 1)

text(y = -1.5*bar.offset,
     x = c(-1*(treeheight+bar.offset), 
           -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))))[1],
     labels = "0",
     #labels = expression(bold("Does not ingest plastic")),
     cex = .8,
     pos = 1)
text(y = -1.5*bar.offset,
     x = c(-1*(treeheight+bar.offset), 
           -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))))[2],
     #labels = expression(bold("Ingests plastic")),
     labels = "1",
     cex = .8,
     pos = 1)

# Lower
lines(x = rep(-1*bar.offset,2), 
      y = c(-1*(treeheight+bar.offset), 
            -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(x = c(-2*bar.offset,-1*bar.offset), 
      y = c(-1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(x = c(-2*bar.offset,-1*bar.offset), 
      y = c(-1*(treeheight+bar.offset+(0.25*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(0.25*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(x = c(-2*bar.offset,-1*bar.offset), 
      y = c(-1*(treeheight+bar.offset+(0.5*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(0.5*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(x = c(-2*bar.offset,-1*bar.offset), 
      y = c(-1*(treeheight+bar.offset+(0.75*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))), 
            -1*(treeheight+bar.offset+(0.75*bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"]))))),
      lwd = 1)
lines(x = c(-2*bar.offset,-1*bar.offset), 
      y = c(-1*(treeheight+bar.offset), 
            -1*(treeheight+bar.offset)),
      lwd = 1)

text(x = -1.5*bar.offset,
     y = c(-1*(treeheight+bar.offset), 
           -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))))[1],
     labels = "0",
     #labels = expression(bold("Does not ingest plastic")),
     cex = .8,
     pos = 2)
# text(x = -1.5*bar.offset,
#      y = c(-1*(treeheight+bar.offset), 
#            -1*(treeheight+bar.offset+(0.5*bar.scale*treeheight)*(max(plot.data[,"mean"]-min(plot.data[,"mean"]))/(max(plot.data[,"mean"])-min(plot.data[,"mean"])))))[2],
#      #labels = expression(bold("Ingests plastic")),
#      labels = "0.5",
#      cex = 1,
#      pos = 2)
text(x = -1.5*bar.offset,
     y = c(-1*(treeheight+bar.offset), 
           -1*(treeheight+bar.offset+(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))))[2],
     #labels = expression(bold("Ingests plastic")),
     labels = "1",
     cex = .8,
     pos = 2)


text(x = -4*bar.offset,
     y = -1*(treeheight+bar.offset+0.15*(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))),
     labels = expression("Mean"),
     cex = .8,
     pos = 2)
text(x = -4*bar.offset,
     y = -1*(treeheight+bar.offset+0.5*(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))),
     labels = expression("Frequency of Occurrence"),
     cex = .8,
     pos = 2)
text(x = -4*bar.offset,
     y = -1*(treeheight+bar.offset+0.85*(bar.scale*treeheight)*(max(plot.data.df[,"Sp_mean"]-min(plot.data.df[,"Sp_mean"]))/(max(plot.data.df[,"Sp_mean"])-min(plot.data.df[,"Sp_mean"])))),
     labels = expression("of plastic ingestion"),
     cex = .8,
     pos = 2)

ggsave(file = "Prelim_phylo_full.png", dpi = 600, width = 8, height = 6, units = "in")



# # Add legend
# legend("bottomleft",
#        legend = row.names(tree.plot.colours),
#        pch=16, pt.cex=1.5, cex = .8, bty='n',
#        col = c(col.Perciformes, col.Tetraodontiformes, col.Lophiiformes, col.Scorpaeniformes,
#                col.Gasterosteiformes, col.Pleuronectiformes, col.Beloniformes, col.Atheriniformes, col.Mugiliformes,
#                col.Ophidiiformes, col.Beryciformes, col.Lampriformes, col.Gadiformes, col.Zeiformes, col.Myctophiformes,
#                col.Aulopiformes, col.Stomiiformes, col.Osmeriformes, col.Salmoniformes, col.Clupeiformes, col.Siluriformes,
#                col.Anguilliformes, col.Carcharhiniformes, col.Squaliformes, col.Rajiformes, col.Myliobatiformes, col.Torpediniformes),
#        title = "Order", title.adj = 0.1)

# if(plot.pdf == TRUE){
#   dev.off()
#   system2(command = "open", args = paste0("C:\\Users\\Danuta\\Dropbox\\matt_phylo\\PHYLOPLAS.Predicted__Figure_",Sys.Date(),".pdf"), wait = FALSE)
# }




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
# # testing one option using ggtree
test.tree <- ggtree(plot.tree, layout = "fan", open.angle = 90) + 
  geom_tiplab(size = 1) + 
  theme_tree2()
