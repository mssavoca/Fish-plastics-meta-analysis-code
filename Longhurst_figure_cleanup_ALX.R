## Alex prepping data for longhurst province figure ##

data <- read.csv("Plastics ingestion records fish master_UDPATED_AGM-MSS2.csv")

# want to get average proportion of plastic per province
summary(data)
head(data)
# data of interest
data2 <- data[,c("Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic", "NwP", "N", "Source")]
head(data2)
colnames(data2) <- c("OceanProv", "PropPlastic", "NwP", "N", "Source")
head(data2)
length(table(data2$OceanProv))
data3<- data2[order(data2$OceanProv),]
head(data3)



### set up new dataframe 
prov <- unique(data3$OceanProv)[-1]
prov
length(prov)
aveplast <- rep(NA, length(prov))
numfish <- rep(NA, length(prov))
numstudies <- rep(NA, length(prov))

## using this to double check the produced averages
sub <- data3[data3$OceanProv=="AUSW",]
sub
unique(sub$Source)
length(unique(sub$Source))
sum(sub$NwP)/sum(sub$N)
mean(sub$PropPlastic) #3 = .1077803
sum(sub$N)
nrow(sub)


for (i in 1:length(prov)){
  sub <- data3[data3$OceanProv==prov[i],]
  aveplast[i] <- sum(sub$NwP)/sum(sub$N)
  numfish[i] <- sum(sub$N)
  numstudies[i] <- length(unique(sub$Source))
}

aveplast
length(aveplast)
numfish
table(numfish)
numstudies


newdat <- data.frame(prov, aveplast, numfish, numstudies)
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

head(newdat)

table(newdat$prov)

head(newdat)
newdat$aveplast
sort(newdat$numstudies)
sort(newdat$numfish)
write.csv(newdat, "Longhurst_FishSummaryData_fullbinned.csv")



