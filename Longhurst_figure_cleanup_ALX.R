## Alex prepping data for longhurst province figure ##
<<<<<<< HEAD
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
=======
>>>>>>> master

data <- read.csv("Plastics ingestion records fish master_UDPATED_AGM-MSS2.csv")

# want to get average proportion of plastic per province
summary(data)
<<<<<<< HEAD
head(data)
# data of interest
data2 <- data[,c("Species.name", "Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic", "NwP", "N", "Source")]
head(data2)
colnames(data2) <- c("Species", "OceanProv", "PropPlastic", "NwP", "N", "Source")
=======

# data of interest
data2 <- data[,c("Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic")]
head(data2)
colnames(data2) <- c("OceanProv", "PropPlastic")
>>>>>>> master
head(data2)
length(table(data2$OceanProv))
data3<- data2[order(data2$OceanProv),]
head(data3)
<<<<<<< HEAD


=======
table(data3$OceanProv)
>>>>>>> master

### set up new dataframe 
prov <- unique(data3$OceanProv)[-1]
prov
<<<<<<< HEAD
length(prov)
aveplast <- rep(NA, length(prov))
numfish <- rep(NA, length(prov))
numstudies <- rep(NA, length(prov))
numspecies <- rep(NA, length(prov))
normalized <- rep(NA, lenght(prov))

## using this to double check the produced averages
sub <- data3[data3$OceanProv=="CHIL",]
sub
unique(sub$Source)
length(unique(sub$Source))
overallprop <- sum(sub$NwP)/sum(sub$N)
overallprop
normalized <- overallprop*sum(sub$N)
normalized
sum(sub$N)
mean(sub$PropPlastic) #3 = .1077803
sum(sub$N)
nrow(sub)
length(unique(sub$Species))



for (i in 1:length(prov)){
  sub <- data3[data3$OceanProv==prov[i],]
  aveplast[i] <- sum(sub$NwP)/sum(sub$N)
  numfish[i] <- sum(sub$N)
  numstudies[i] <- length(unique(sub$Source))
  numspecies[i] <- length(unique(sub$Species))
  normalized[i] <- aveplast[i]*sum(sub$N)
}

aveplast
length(aveplast)
numfish
table(numfish)
numstudies
head(data3)
numspecies


newdat <- data.frame(prov, aveplast, numfish, numstudies, numspecies, normalized)
newdat

## binning data
# ave plastic
newdat$aveplastbin <- rep(NA, length(prov))
newdat$aveplastbin[newdat$aveplast<=.10]= "<.10"
newdat$aveplastbin[newdat$aveplast >.10 & newdat$aveplast<=20]= ".11-.20"
newdat$aveplastbin[newdat$aveplast >.20 & newdat$aveplast<=.30]= ".21-.30"
newdat$aveplastbin[newdat$aveplast >.30 & newdat$aveplast<=.40]= ".31-.40"
newdat$aveplastbin[newdat$aveplast >.40 & newdat$aveplast<=.50]= ".41-.50"
newdat$aveplastbin[newdat$aveplast >.50 & newdat$aveplast<=.60]= ".51-.60"
newdat$aveplastbin[newdat$aveplast >.60 & newdat$aveplast<=.70]= ".61-.70"
newdat$aveplastbin[newdat$aveplast >.70 & newdat$aveplast<=.80]= ".71-.80"
newdat$aveplastbin[newdat$aveplast >.80 & newdat$aveplast<=.90]= ".81-.90"
newdat$aveplastbin[newdat$aveplast >.90 & newdat$aveplast<=.99]= ".90-.99"
newdat$aveplastbin[newdat$aveplast ==1]= "1"

head(newdat$aveplastbin, 50)
sort(newdat$aveplastbin)
# number of studies
summary(newdat$numstudies)
newdat$numstudiesbin <- rep(NA, length(prov))
newdat$numstudiesbin[newdat$numstudies == 1] = "1"
newdat$numstudiesbin[newdat$numstudies >1 & newdat$numstudies<=5]= "2-5"
newdat$numstudiesbin[newdat$numstudies >5 & newdat$numstudies<=10]= "6-10"
newdat$numstudiesbin[newdat$numstudies >=10]= ">10"

head(newdat, 50)

# number of fish/study
summary(newdat$numfish)
newdat$numfishbin <- rep(NA, length(prov))
newdat$numfishbin[newdat$numfish >1 & newdat$numfish<10]= "< 10"
newdat$numfishbin[newdat$numfish >=10 & newdat$numfish<=50]= "10-50"
newdat$numfishbin[newdat$numfish >50 & newdat$numfish<=100]= "51-100"
newdat$numfishbin[newdat$numfish >100 & newdat$numfish<500]= "101-500"
newdat$numfishbin[newdat$numfish >500 & newdat$numfish<=1000]= "501-1000"
newdat$numfishbin[newdat$numfish >1000 & newdat$numfish<=1500]= "1001-1500"
newdat$numfishbin[newdat$numfish >1500]= ">1500"

head(newdat, 50)

# normalized proportion of ingestion bins
summary(newdat$normalized)
newdat$normbin <- rep(NA, length(prov))
newdat$normbin[newdat$normalized >=1 & newdat$normalized<25]= "< 10"
newdat$normbin[newdat$normalized >=10 & newdat$normalized<=50]= "10-50"
newdat$normbin[newdat$normalized >50 & newdat$normalized<=100]= "51-100"
newdat$normbin[newdat$normalized >100 & newdat$normalized<500]= "101-500"
newdat$normbin[newdat$normalized >500 & newdat$normalized<=1000]= "501-1000"
newdat$normbin[newdat$normalized >1000 & newdat$normalized<=1500]= "1001-1500"
newdat$normbin[newdat$normalized >1500]= ">1500"




head(newdat, 50)
newdat$aveplast
sort(newdat$numstudies)
sort(newdat$numfish)
write.csv(newdat, "Longhurst_FishSummaryData_fullbinned.csv")


newdat

### Summary table:
# geographic summary of data
Geo_summ_pt2 <- newdat %>% 
  group_by(`prov`) %>% 
  summarize(num_studies = numstudies,
            num_sp = numspecies,
            num_ind_studied = numfish)

Geo_summ_pt2
nrow(Geo_summ_pt2)
ncol(Geo_summ_pt2)

write.table(Geo_summ_pt2, file = "Geo_summ_pt2.txt", sep = ",", quote = FALSE, row.names = F)
=======
aveplast <- rep(NA, length(prov)) 
na.omit(data3)
table(data3$OceanProv)

## using this to double check the produced averages
sub <- data3[data3$OceanProv=="BRAZ",]
sub
mean(sub$PropPlastic) #3 = .1077803

for (i in 1:length(prov)){
  sub <- data3[data3$OceanProv==prov[i],]
  aveplast[i] <- mean(sub$PropPlastic)
}

aveplast

newdat <- data.frame(prov, aveplast)
write.csv(newdat, file="AveragePlasticCount.csv")

# now to rename according to the QGIS shape file
# BPRL --> BPLR
# HUMB == CHIL on this map
# NAST E --> NASE
>>>>>>> master
