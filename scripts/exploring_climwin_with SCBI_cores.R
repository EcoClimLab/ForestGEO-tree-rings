
# clear environment ####
rm(list = ls())

# load libraries ####
library(climwin)

# load data ####

# core data
cores <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_cores/cross-dated_cores_CSVformat/all_core_chronologies.csv")

# climate data
climate_variables <- c( "cld", "dtr", "frs", 
                        "pet", "pre", 
                        "tmn", "tmp", "tmx", 
                        "vap", "wet")
for(clim_v in climate_variables) {
  assign(clim_v, read.csv(paste0("https://raw.githubusercontent.com/forestgeo/Climate/master/Gridded_Data_Products/Historical%20Climate%20Data/CRU_v4_01/", clim_v, ".1901.2016-ForestGEO_sites-8-18.csv"))
  )
}

# prepare data ####

## climate data ####

for(clim_v in climate_variables) {
  print(clim_v)
  x <- get(clim_v)
  
  ### subset for SCBI data only
  x <- x[x$sites.sitename %in% "Smithsonian_Conservation_Biology_Institute_(SCBI)", ]

  ### reshape to long format
  x_long <- reshape(x, 
                      times = names(x)[-1], timevar = "Date",
                      varying = list(names(x)[-1]), direction = "long", v.names = clim_v)
  
  
  ### combine all variables in one
  if(clim_v == climate_variables [1]) Clim <- x_long[, c(2,3)]
  else Clim <- merge(Clim, x_long[, c(2,3)], by = "Date", all = T)
  
}
  
### format date to dd/mm/yyyy
Clim$Date <- gsub("X", "", Clim$Date)
Clim$Date <- format(as.Date(Clim$Date , format = "%Y.%m.%d"), "%d/%m/%Y") 


## cores ####

### reshape to long format
Biol <- reshape(cores, idvar = c("tag",  "coreID",   "sp", "status.at.coring"),
        times = names(cores)[-(1:4)], timevar = "Date",
        varying = list(names(cores)[-(1:4)]), direction = "long", v.names = "core_measurement")

### format date to dd/mm/yyyy
Biol$Date <- gsub("X", "", Biol$Date)
Biol$Date <- paste0("15/06/", Biol$Date) # dd/mm/yyyy setting up as june 15, ARBITRATY

head(cores)
head(Biol)
tail(Biol)

cores[cores$tag %in% 182024, "X2017"]

### remove years that are not in climate record (+ first year to be able to look at window before measurement)
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) >= min(as.numeric(substr(Clim$Date, 7, 10)))+1, ]
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) <= max(as.numeric(substr(Clim$Date, 7, 10))), ]

# remove NA
Biol <- Biol[!is.na(Biol$core_measurement), ]


# try to run slidingwin ####


results <- slidingwin( baseline = lm(core_measurement ~ 1 + sp, data = Biol),
                       xvar =list(cld = Clim$cld, 
                                  dtr = Clim$dtr, 
                                  frs = Clim$frs, 
                                  pet = Clim$pet, 
                                  pre = Clim$pre, 
                                  tmn = Clim$tmn, 
                                  tmp = Clim$tmp, 
                                  tmx = Clim$tmx, 
                                  vap = Clim$vap, 
                                  wet = Clim$wet),
                       type = "absolute", 
                       range = c(12, 0),
                       stat = c("mean"),
                       func = c("lin","quad"),
                       refday = c(30, 6),
                       cinterval = "month",
                       cdate = Clim$Date, bdate = Biol$Date) 


randomized1<-randwin(repeats = 10,     
                     baseline = lm(core_measurement ~ 1 + sp, data = Biol),
                     xvar = list(Clim$pre),
                     type = "absolute", 
                     range = c(12, 0),
                     stat = c("mean"),
                     func = c("quad"),
                     refday=c(30, 6),
                     cinterval = "month",
                     cdate = Clim$Date, bdate = Biol$Date,
                     window= "sliding")


output <- results
pvalue(datasetrand = randomized1[[1]], dataset = output[[4]]$Dataset, metric = "AIC", sample.size = nrow(Biol))


plotall(datasetrand = randomized1[[1]],
        dataset = output[[4]]$Dataset, 
        bestmodel = output[[4]]$BestModel,
        bestmodeldata = output[[4]]$BestModelData,
        title=output$combos[4,])

