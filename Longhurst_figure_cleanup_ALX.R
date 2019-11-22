## Alex prepping data for longhurst province figure ##
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

data <- read.csv("Plastics ingestion records fish master_final.csv")
nrow(data)
# want to get average proportion of plastic per province
summary(data)
head(data)
# data of interest
data2 <- data[,c("Binomial", "Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic", "NwP", "N", "Source")]
head(data2)
colnames(data2) <- c("Species", "OceanProv", "PropPlastic", "NwP", "N", "Source")
head(data2)
length(table(data2$OceanProv))
data3<- data2[order(data2$OceanProv),]
head(data3)

levels(data3$OceanProv) <- c(levels(data3$OceanProv), "CHIL", "BPLR", "NASE", "NPPF") # need to change these to match the publicly available shape file
data3$OceanProv[data3$OceanProv=="BPRL"] <- "BPLR"
data3$OceanProv[data3$OceanProv=="HUMB"] <- "CHIL"
data3$OceanProv[data3$OceanProv=="NAST E"] <- "NASE"
data3$OceanProv[data3$OceanProv=="NPSE"] <- "NPPF"

table(data3$OceanProv)
data3$OceanProv <- droplevels(data3$OceanProv)
levels(data3$OceanProv)
### set up new dataframe 
prov <- unique(data3$OceanProv)
prov
length(prov)
aveplast <- rep(NA, length(prov))
numfish <- rep(NA, length(prov))
numstudies <- rep(NA, length(prov))
numspecies <- rep(NA, length(prov))
normalized <- rep(NA, length(prov))

## using this to double check the produced averages
sub <- data3[data3$OceanProv=="CHIN",]
sub
unique(sub$Source)
length(unique(sub$Source))
overallprop <- sum(sub$NwP, na.rm=T)/sum(sub$N, na.rm=T)
overallprop
normalized <- overallprop*sum(sub$N, na.rm=T)
normalized
sum(sub$N, na.rm=T)
mean(sub$PropPlastic, na.rm=T) #3 = .1077803
sum(sub$N)
nrow(sub)
length(unique(sub$Species))



for (i in 1:length(prov)){
  sub <- data3[data3$OceanProv==prov[i],]
  aveplast[i] <- sum(sub$NwP, na.rm=T)/sum(sub$N, na.rm=T)
  numfish[i] <- sum(sub$N, na.rm=T)
  numstudies[i] <- length(unique(sub$Source))
  numspecies[i] <- length(unique(sub$Species))
  normalized[i] <- aveplast[i]*sum(sub$N, na.rm=T)
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
head(newdat)



#write.csv(newdat, "Longhurst_FishSummaryData_fullbinned.csv")


############################ Spatial Data Analysis/Mapping #################################
library(sf)
library(rgdal)
library(raster)
library(maptools)
#note that to read in the shape file you need the full address (https://datacarpentry.org/r-raster-vector-geospatial/06-vector-open-shapefile-in-r/)
# ggplot2 will only work with a data.frame object

lh_prov <- st_read("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/longhurst_v4_2010/Longhurst_world_v4_2010.shp")
lh_prov
crs(lh_prov)
lh_prov.2 <- fortify(lh_prov) #fortify makes this into a data frame object

## adding adjacency values: 1 for LH province touching major landmass (not islands), 0 for not touching major landmass

lh_prov.2$adjacency <- ifelse(lh_prov.2$ProvCode=="ALSK"| lh_prov.2$ProvCode=="CCAL"|lh_prov.2$ProvCode=="CAMR"|lh_prov.2$ProvCode=="BERS"|lh_prov.2$ProvCode=="BPLR"|lh_prov.2$ProvCode=="CARB"|lh_prov.2$ProvCode=="CHIL"|lh_prov.2$ProvCode=="FKLD"|lh_prov.2$ProvCode=="BRAZ"|lh_prov.2$ProvCode=="APLR"|lh_prov.2$ProvCode=="EAFR"|lh_prov.2$ProvCode=="BENG"|lh_prov.2$ProvCode=="GUIN"|lh_prov.2$ProvCode=="NWCS"|lh_prov.2$ProvCode== "CNRY"|lh_prov.2$ProvCode=="MEDI"|lh_prov.2$ProvCode== "NECS"|lh_prov.2$ProvCode=="SARC"|lh_prov.2$ProvCode=="REDS"|lh_prov.2$ProvCode== "ARAB"|lh_prov.2$ProvCode== "INDW"|lh_prov.2$ProvCode== "INDE"|lh_prov.2$ProvCode== "AUSW"|lh_prov.2$ProvCode== "BERS"|lh_prov.2$ProvCode=="SUND"|lh_prov.2$ProvCode== "AUSE"|lh_prov.2$ProvCode== "KURO"|lh_prov.2$ProvCode== "CHIN"|lh_prov.2$ProvCode=="TASM"|lh_prov.2$ProvCode=="NEWZ"|lh_prov.2$ProvCode=="SPSG"|lh_prov.2$ProvCode== "ARCH", 1, 0)
sum(lh_prov.2$adjacency)

lh_prov.2$adjacency
head(lh_prov.2)


####### calculating centroids for labels ###
centroids.df <- as.data.frame(st_centroid(lh_prov.2)) 
head(centroids.df)


joined<- right_join(lh_prov.2,newdat, by=c("ProvCode"="prov"))
head(joined)
class(joined)
centroids.df <- as.data.frame(st_centroid(joined))
centroids.df
centroids.coords <- st_coordinates(centroids.df$geometry)
centroids.coords 
nrow(centroids.coords)
### combining all of the polygon data for the provinces for which we have data
full_data <- cbind(joined, centroids.coords)

class(full_data)
full_data #X and Y correspond to the centroid of the polygon
class(full_data)
nrow(full_data)

######### Base map #### 
lh_map <- ggplot(lh_prov) + 
  geom_sf() +
  coord_sf() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     # Remove panel background
                     panel.background = element_blank())
#lh_map
class(lh_map)

## mapping with attributes


## average plastic binned
library(viridis)
aveplastmap<- lh_map+
              geom_sf(data=full_data, aes(fill=aveplastbin)) + scale_fill_brewer(palette="YlOrRd") + labs(fill = "Average Plastic Ingested") #+ geom_text(aes(label = full_data$ProvCode, x = full_data$X, y = full_data$Y))

#aveplastmap

### done here. Need to add labels

## number of studies binned
numstudiesmap <-  lh_map+
              geom_sf(data=joined, aes(fill=aveplastbin)) + scale_fill_brewer(palette="YlGn")
#numstudiesmap

######### average plastic pollution map ################
abund <- read.csv("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/Global pollution_map files/vansebillemodel_abundance.csv")
poll_lat <- read.csv("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/Global pollution_map files/latitudes.csv", header = FALSE)
poll_lon <- read.csv("~/Box Sync/Microplastics/Fish-plastics-meta-analysis-code/Global pollution_map files/longitudes.csv", stringsAsFactors = FALSE)

####### Creating new object for plastic rasters
head(abund)
head(poll_lat)
head(poll_lon)

x_vals <- as.numeric(gsub("X", "", colnames(poll_lon))) #gsub = regular expression that will drop the X from in front of the numeric values in our x data (also, lon reads in as column names)
x_vals[x_vals > 179]  = x_vals[x_vals > 179] - 360 #need to shift the view on the map, so essentially flip the map around the 180. Probably need to get rid of 360 because 0 and 360 are redundant


###### IF YOU WANT TO GET RID OF THE LINE AT 0/360 or 0/0 (if we are looking at 0-180 scale), run this code instead ###
# Creating raster list to generate full object 
#raster_obj <- list(z = as.matrix(abund)[, order(x_vals[-1])], #list with attributes as z, and coordinates as x and y
               #    x = x_vals[-1], #got rid of the initial 0 value, because 0 and 360 (or in this case, if we are going between -180 and 180, 0 and 0) overlap. This gets rid of the weird disconnect in the middle of the map
              #    y = as.numeric(poll_lat$V1))
##############################################################################################################

#to calculate average plastic, will keep the weird line in because we want the plastic values at 0 and 360
raster_obj <- list(z = as.matrix(abund)[, order(x_vals)], #list with attributes as z, and coordinates as x and y
                   x = x_vals, #got rid of the initial 0 value, because 0 and 360 (or in this case, if we are going between -180 and 180, 0 and 0) overlap. This gets rid of the weird disconnect in the middle of the map
                   y = as.numeric(poll_lat$V1))

# Create a new raster() file
poll_raster <- raster(x = raster_obj$z, # Matrix values of plastic
                      
       # Defining endpoints of the raster
       xmn = min(raster_obj$x),
       xmx = max(raster_obj$x),
       ymn = min(raster_obj$y),
       ymx = max(raster_obj$y),
       
       # Setting coordinate reference system to be equivalent to longhurst
       crs = crs(lh_prov))

# Plotting Z on the log10 scale (Van Sebille et al (2015, ERL) paper); see that this works; don't want calculations on this scale
plot(poll_raster)

plot(st_geometry(lh_prov), add = TRUE, fill = NULL)

### Next step: bring out average attribute per polygon
extracted_vals <- extract(poll_raster, lh_prov) #this extracts values in the poll_raster per polygon
str(extracted_vals) #this should give us a list of all of the values in each of 54 different LH objects
mean_poll_abund<- unlist(lapply(extracted_vals, mean, na.rm=T)) #this should give us the mean values for each of the 54 LH objects
lapply(extracted_vals, range, na.rm=T) #note that these units are #/km^2

fullmapdat<- cbind(lh_prov.2, mean_poll_abund)
fullmapdat

fullmap.df <- fortify(fullmapdat)
fullmap.df <- st_drop_geometry(fullmap.df) #need to remove geometry in order to get an actual dataframe (zero spatial component)
head(fullmap.df)
nrow(fullmap.df)

full_sampled_data <- merge(fullmap.df, newdat, by.x="ProvCode", by.y="prov", sort = TRUE)

write.csv(fullmap.df, "Spatial Information_microplastics.csv")
write.csv(full_sampled_data, "Spatial Information_onlysampled.csv")

############# PRETTIER MAP OF POLLUTION (VAN SEBILLE et al. 2015 data) #######################################
## caveat - the plotting function of raster retains aspect ratio (only need to change this if we are plotting)
# need to convert pollution to sf package (multipolygon, or spatial polygons data frame...right now, just a raster)
spdf <- as(log10(poll_raster),'SpatialPolygonsDataFrame')
spdf.real <- as(poll_raster,'SpatialPolygonsDataFrame') # not log-transformed, harder to see
#--- convert to an sf object ---#
sf_poll_data <- st_as_sf(spdf.real)
#--- take a look ---#
plot(sf_poll_data, border = FALSE, xlim = c(-180, 180), ylim = c(-90, 90), reset = FALSE, logz = TRUE) # if you want to get rid of the weird line, scroll up and re-create raster object where indicated
plot(st_geometry(lh_prov), add = TRUE, fill = NULL, xlim = c(-180, 180), ylim = c(-90, 90))

