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
table(data3)

### set up new dataframe 
prov <- unique(data3$OceanProv)[-1]
aveplast <- rep(NA, length(prov)) 
newdat <- data.frame(prov, aveplast)
table(data3$OceanProv)

for (i in 1:nrow(newdat)){
  
  
}
