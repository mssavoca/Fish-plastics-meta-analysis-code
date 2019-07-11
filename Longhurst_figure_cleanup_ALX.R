## Alex prepping data for longhurst province figure ##

data <- read.csv("Plastics ingestion records fish master_UDPATED_AGM-MSS2.csv")

# want to get average proportion of plastic per province
summary(data)

# data of interest
data2 <- data[,c("Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic")]
head(data2)
colnames(data2) <- c("OceanProv", "PropPlastic")
head(data2)
length(table(data2$OceanProv))
data3<- data2[order(data2$OceanProv),]
head(data3)
table(data3$OceanProv)

### set up new dataframe 
prov <- unique(data3$OceanProv)[-1]
prov
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