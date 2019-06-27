## Alex prepping data for longhurst province figure ##

data <- read.csv("Plastics ingestion records fish master_UDPATED_AGM-MSS2.csv")

# want to get average proportion of plastic per province
summary(data)
head(data)
# data of interest
data2 <- data[,c("Oceanographic.province..from.Longhurst.2007.", "Prop.w.plastic", "N")]
head(data2)
colnames(data2) <- c("OceanProv", "PropPlastic", "N")
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
sub <- data3[data3$OceanProv=="NPTG",]
sub
mean(sub$PropPlastic) #3 = .1077803
sum(sub$N)
nrow(sub)


for (i in 1:length(prov)){
  sub <- data3[data3$OceanProv==prov[i],]
  aveplast[i] <- mean(sub$PropPlastic)
  numfish[i] <- sum(sub$N)
  numstudies[i] <- nrow(sub)
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
newdat$aveplastbin[newdat$aveplast <.05] = "<.05"
newdat$aveplastbin[newdat$aveplast >.05 & newdat$aveplast<.10]= ".05-.10)"
newdat$aveplastbin[newdat$aveplast >=.10 & newdat$aveplast<.20]= ".10-.20)"
newdat$aveplastbin[newdat$aveplast >=.20 & newdat$aveplast<.30]= ".20-.30)"
newdat$aveplastbin[newdat$aveplast >=.30 & newdat$aveplast<.40]= ".30-.40)"
newdat$aveplastbin[newdat$aveplast >=.40 & newdat$aveplast<.50]= ".40-.50)"
newdat$aveplastbin[newdat$aveplast >=.50 & newdat$aveplast<.60]= ".50-.60)"
newdat$aveplastbin[newdat$aveplast >=.60 & newdat$aveplast<.70]= ".60-.70)"
newdat$aveplastbin[newdat$aveplast >=.70 & newdat$aveplast<.80]= ".70-.80)"
newdat$aveplastbin[newdat$aveplast >=.80 & newdat$aveplast<.90]= ".80-.90)"
newdat$aveplastbin[newdat$aveplast >=.90 & newdat$aveplast<1]= ".90-.1)"
newdat$aveplastbin[newdat$aveplast ==1]= "1"

head(newdat, 50)

# number of studies
summary(newdat$numstudies)
newdat$numstudiesbin <- rep(NA, length(prov))
newdat$numstudiesbin[newdat$numstudies == 1] = "1"
newdat$numstudiesbin[newdat$numstudies >1 & newdat$numstudies<10]= "1-10)"
newdat$numstudiesbin[newdat$numstudies >=10 & newdat$numstudies<25]= "10-25)"
newdat$numstudiesbin[newdat$numstudies >=25 & newdat$numstudies<50]= "25-50)"
newdat$numstudiesbin[newdat$numstudies >=50 & newdat$numstudies<75]= "50-75)"
newdat$numstudiesbin[newdat$numstudies >=75 & newdat$numstudies<100]= "75-100)"
newdat$numstudiesbin[newdat$numstudies >=100]= "100+"

head(newdat, 50)

# number of fish/study
summary(newdat$numfish)
newdat$numfishbin <- rep(NA, length(prov))
newdat$numfishbin[newdat$numfish >1 & newdat$numfish<10]= "< 10"
newdat$numfishbin[newdat$numfish >=10 & newdat$numfish<50]= "10-50)"
newdat$numfishbin[newdat$numfish >=50 & newdat$numfish<100]= "50-100)"
newdat$numfishbin[newdat$numfish >=100 & newdat$numfish<500]= "100-500)"
newdat$numfishbin[newdat$numfish >=500 & newdat$numfish<1000]= "500-1000)"
newdat$numfishbin[newdat$numfish >=1000 & newdat$numfish<1500]= "1000-1500)"
newdat$numfishbin[newdat$numfish >=1500]= "1500+"

head(newdat, 50)

head(newdat)
write.csv(newdat, "Longhurst_FishSummaryData_fullbinned.csv")



