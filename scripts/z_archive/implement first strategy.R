# implement the plan which is:
# - use climwin to identify the strongest temperature related variable and moisture related variables signals (and time windoe) using mixed model with species and core ID (not nested) as random effects and the individual-levels residuals of core measurements spline.
# - using the best vriable, use a gam on the raw data, with a smoothing on dbh and on, Year at the individual level fo year. On model for each species.
# help about specifying gams https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.models.html


# clear environment ####
rm(list = ls())

# load libraries ####
library(climwin)
library(lme4)
library(mgcv)

# set parameters ####
core_type <- "CSV"
reference_date <- c(30, 8) # refday in slidingwin
window_range <- c(15, 0) #range in slidingwin


baseline = "lmer(core_measurement ~ 1 + dbh + (1 | sp) + (1 | coreID), data = Biol)"
baseline_by_sp <- "lmer(core_measurement ~ 1 + dbh + (1 | coreID), data = Biol[Biol$sp %in% sp,])"


# load data ####

## core data ####


cores <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_cores/cross-dated_cores_CSVformat/all_core_chronologies.csv")

dbh_2008 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem1.csv", stringsAsFactors = FALSE)

dbh_2013 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem2.csv", stringsAsFactors = FALSE)

dbh_2018 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem3.csv", stringsAsFactors = FALSE)


## climate data ####
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

Biol$Year <- as.numeric(format(as.Date(Biol$Date, format = "%d/%m/%Y"), "%Y"))

head(cores)
head(Biol)
tail(Biol)

cores[cores$tag %in% 182024, "X2017"]

### remove NAs
Biol <- Biol[!is.na(Biol$core_measurement), ]


### find out DBH for each year
Biol$dbh = NA

for( t in unique(Biol$tag)) {
  print(t)
  print(which(unique(Biol$tag) == t))
  
  x <- Biol[Biol$tag %in% t , ]
  max(x$Year)
  
  y_dbh <- c("2008" = dbh_2008[dbh_2008$tag %in% t, ]$dbh, 
             "2013" = dbh_2013[dbh_2013$tag %in% t, ]$dbh, 
             "2018" = dbh_2018[dbh_2018$tag %in% t, ]$dbh)
  y_status <- c("2008" = dbh_2008[dbh_2008$tag %in% t, ]$status, 
                "2013" = dbh_2013[dbh_2013$tag %in% t, ]$status, 
                "2018" = dbh_2018[dbh_2018$tag %in% t, ]$status)
  
  if(sum(dbh_2008$tag %in% t) > 1 ) { # if there is more than one stem, keep main stem
    y_StemTag <- c("2008" = dbh_2008[dbh_2008$tag %in% t, ]$StemTag, 
                   "2013" = dbh_2013[dbh_2013$tag %in% t, ]$StemTag, 
                   "2018" = dbh_2018[dbh_2018$tag %in% t, ]$StemTag)
    
    StemTag_main <- dbh_2008[dbh_2008$tag %in% t, ]$StemTag[grepl("M|NULL", dbh_2008[dbh_2008$tag %in% t, ]$codes)]
    
    y_dbh <- y_dbh[y_StemTag == StemTag_main]
    y_status <- y_status[y_StemTag == StemTag_main]
    
    names(y_dbh) <- names(y_status) <- c("2008", "2013", "2018")
    
  }
  
  closest.year.tree.was.cenused.alive <- names(which(y_status == "A"))[which.min(abs(as.numeric(names(which(y_status == "A"))) - max(x$Year)))]
  
  # add last DBH measured to the data
  
  ## if that year exist in the core measurement, just do it
  if(closest.year.tree.was.cenused.alive %in% x$Year)  x[x$Year == closest.year.tree.was.cenused.alive, ]$dbh <- as.numeric(y_dbh[closest.year.tree.was.cenused.alive])
  
  ## if that year does not exist in hte core measurement (cored while alive and in a previous census), calculate average ring increment in the past few years and use that to estimate the dbh on the last year of core
  
  if(!closest.year.tree.was.cenused.alive %in% x$Year) {
    nb_years_to_average_from <- as.numeric(closest.year.tree.was.cenused.alive) - max(x$Year)
    avg_core_increment <- mean(rev(x$core_measurement)[1:nb_years_to_average_from])
    
    x[x$Year ==  max(x$Year), ]$dbh <- as.numeric(y_dbh[closest.year.tree.was.cenused.alive]) - 2*avg_core_increment*nb_years_to_average_from
    
    closest.year.tree.was.cenused.alive <- max(x$Year)
  }  
  
  # populate dbh before measurement
  idx_x_before <- x$Year <= closest.year.tree.was.cenused.alive 
  x[idx_x_before, ]$dbh <-  rev(x[x$Year == closest.year.tree.was.cenused.alive, ]$dbh - c(0, cumsum(  rev(2*x[idx_x_before, ]$core_measurement))))[-1]
  
  # populate dbh after measurement
  idx_x_after <- x$Year > closest.year.tree.was.cenused.alive
  x[idx_x_after, ]$dbh <-  x[x$Year == closest.year.tree.was.cenused.alive, ]$dbh + cumsum(  2*x[idx_x_after, ]$core_measurement)
  
  # double check diff dbh corresponds to core measurements
  if(!all(round(diff(x$dbh),3) - 2* round(x$core_measurement[-1], 3) == 0)) stop("dbh calculated invcorrectly")
  
  # save into Biol
  
  Biol[Biol$tag %in% t , ] <- x
}
  

### remove years that are before climate record (+ first few first months to be able to look at window before measurement)
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) >= min(as.numeric(substr(Clim$Date, 7, 10)))+  window_range[1]/12, ]

### remove years that are after climate record
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) <= max(as.numeric(substr(Clim$Date, 7, 10))), ]


## calculate residuals of spine measurement ~ year for each individual####

Biol$residuals <- NA

for( t in unique(Biol$coreID)) {
  x <- Biol[Biol$coreID %in% t, ]
  
  test <- gam(core_measurement~ s(Year), data = x)
  par(mfrow = c(3,2))
  plot(test)
  gam.check(test,pch=19,cex=.3)
  
  plot(core_measurement~ Year, data = x, main = "Raw data")
  points(test$fitted.values~ x$Year, type = "l")
  
  title(paste(x$sp[1], x$status.at.coring[1], t, sep = " - " ), outer = T, line = -2)
  
  # # save plot
  # dev.print(tiff, paste0('results/explorations/by_tag/', paste(x$sp[1], x$status.at.coring[1], t, sep = "_" ), "_residuals_core_meas_Year_GAM", '.tif'),
  #           width = 8,
  #           height =8,
  #           units = "in",
  #           res = 300)
  
  # save residuals
  x$residuals <- test$residuals 
  
  # save back into Biol
  Biol[Biol$coreID %in% t, ] <- x
  
}

## run slidingwin on all species ####
baseline = "lmer(residuals ~ 1 + dbh + (1 | sp) + (1 | coreID), data = Biol)"

### for temperature variables only ####
results_temp <- slidingwin( baseline = eval(parse(text = baseline)),
                       xvar =list(dtr = Clim$dtr,
                                  pet = Clim$pet, 
                                  tmn = Clim$tmn, 
                                  tmp = Clim$tmp, 
                                  tmx = Clim$tmx
                       ),
                       type = "absolute", 
                       range = window_range,
                       stat = c("mean"),
                       func = c("lin","quad"),
                       refday = reference_date,
                       cinterval = "month",
                       cdate = Clim$Date, bdate = Biol$Date) 

results_temp$combos
best_mod_first_step <- which.min(results_temp$combos$DeltaAICc)
results_temp$combos[best_mod_first_step,]
results_temp[[best_mod_first_step]][[1]]

randomized_temp <- randwin(repeats = 100,     
                       baseline = eval(parse(text = baseline)),
                       xvar = list(Clim[,as.character(results_temp$combos$climate[best_mod_first_step])]),
                       type = results_temp$combos$type[best_mod_first_step], 
                       range = window_range,
                       stat = results_temp$combos$stat[best_mod_first_step],
                       func = results_temp$combos$func[best_mod_first_step],
                       refday = reference_date,
                       cinterval = "month",
                       cdate = Clim$Date, bdate = Biol$Date,
                       window= "sliding")


output <- results_temp
pvalue(datasetrand = randomized_temp[[1]], dataset = output[[best_mod_first_step]]$Dataset, metric = "AIC", sample.size = nrow(Biol))


plotall(datasetrand = randomized_temp[[1]],
        dataset = output[[best_mod_first_step]]$Dataset, 
        bestmodel = output[[best_mod_first_step]]$BestModel,
        bestmodeldata = output[[best_mod_first_step]]$BestModelData,
        title=paste((data.frame(lapply(output$combos[best_mod_first_step,], as.character), stringsAsFactors=FALSE)), collapse = "_"))


### for moisture variables only ####
results_moist <- slidingwin( baseline = eval(parse(text = baseline)),
                            xvar =list(cld = Clim$cld, 
                                       pre = Clim$pre, 
                                       wet = Clim$wet
                            ),
                            type = "absolute", 
                            range = window_range,
                            stat = c("mean"),
                            func = c("lin","quad"),
                            refday = reference_date,
                            cinterval = "month",
                            cdate = Clim$Date, bdate = Biol$Date) 

results_moist$combos
best_mod_first_step <- which.min(results_moist$combos$DeltaAICc)
results_moist$combos[best_mod_first_step,]
results_moist[[best_mod_first_step]][[1]]

randomized_moist <- randwin(repeats = 100,     
                           baseline = eval(parse(text = baseline)),
                           xvar = list(Clim[,as.character(results_moist$combos$climate[best_mod_first_step])]),
                           type = results_moist$combos$type[best_mod_first_step], 
                           range = window_range,
                           stat = results_moist$combos$stat[best_mod_first_step],
                           func = results_moist$combos$func[best_mod_first_step],
                           refday = reference_date,
                           cinterval = "month",
                           cdate = Clim$Date, bdate = Biol$Date,
                           window= "sliding")


output <- results_moist
pvalue(datasetrand = randomized_moist[[1]], dataset = output[[best_mod_first_step]]$Dataset, metric = "AIC", sample.size = nrow(Biol))


plotall(datasetrand = randomized_moist[[1]],
        dataset = output[[best_mod_first_step]]$Dataset, 
        bestmodel = output[[best_mod_first_step]]$BestModel,
        bestmodeldata = output[[best_mod_first_step]]$BestModelData,
        title=paste((data.frame(lapply(output$combos[best_mod_first_step,], as.character), stringsAsFactors=FALSE)), collapse = "_"))



## save outputs ####
# dev.print(tiff, paste0('results/ALL_species_mixed_model_', paste((data.frame(lapply(output$combos[best_mod_first_step,], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.tif'),
#           width = 10,
#           height =8,
#           units = "in",
#           res = 300)


# save(results, best_mod_first_step, file = paste0('results/ALL_species_mixed_model_', ifelse(core_type == "CSV", "core_meas", "res"), ".RData"))
# 
# write.csv(results$combos[order(results$combos$DeltaAICc),], 
#           file = paste0('results/ALL_species_mixed_model_', ifelse(core_type == "CSV", "core_meas_", "res_"), "combos.csv"))


# run slidingwin one species at a time ####
best_mod_first_step_by_sp <- NULL

for ( sp in species) {
  print(sp)
  if(core_type == "CSV") sp <- tolower(sp)
  
  results <- slidingwin( baseline = eval(parse(text = baseline_by_sp)),
                         xvar =list(cld = Clim$cld, 
                                    # dtr = Clim$dtr, 
                                    # frs = Clim$frs,# removed wet because model failed to converge 
                                    pet = Clim$pet, 
                                    pre = Clim$pre, 
                                    tmn = Clim$tmn, 
                                    tmp = Clim$tmp, 
                                    tmx = Clim$tmx, 
                                    # vap = Clim$vap, 
                                    wet = Clim$wet
                         ),
                         type = "absolute", 
                         range = window_range,
                         stat = c("mean"),
                         func = c("lin","quad"),
                         refday = reference_date,
                         cinterval = "month",
                         cdate = Clim$Date, bdate = Biol$Date[Biol$sp %in% sp]) 
  
  results$combos
  best_mod_first_step <- which.min(results$combos$DeltaAICc)
  results$combos[best_mod_first_step,]
  results[[best_mod_first_step]][[1]]
  
  best_mod_first_step_by_sp <- rbind(best_mod_first_step_by_sp, data.frame(sp = sp, results$combos[best_mod_first_step,]))
  
  randomized1 <- randwin(repeats = 100,     
                         baseline = eval(parse(text = baseline_by_sp)),
                         xvar = list(Clim[,as.character(results$combos$climate[best_mod_first_step])]),
                         type = results$combos$type[best_mod_first_step], 
                         range = window_range,
                         stat = results$combos$stat[best_mod_first_step],
                         func = results$combos$func[best_mod_first_step],
                         refday = reference_date,
                         cinterval = "month",
                         cdate = Clim$Date, bdate = Biol$Date[Biol$sp %in% sp],
                         window= "sliding")
  
  
  output <- results
  pvalue(datasetrand = randomized1[[1]], dataset = output[[best_mod_first_step]]$Dataset, metric = "AIC", sample.size = nrow(Biol))
  
  
  plotall(datasetrand = randomized1[[1]],
          dataset = output[[best_mod_first_step]]$Dataset, 
          bestmodel = output[[best_mod_first_step]]$BestModel,
          bestmodeldata = output[[best_mod_first_step]]$BestModelData,
          title = paste(sp, paste((data.frame(lapply(output$combos[best_mod_first_step,], as.character), stringsAsFactors=FALSE)), collapse = "_"), sep = " - ")
  )
  
  
  ## save outputs ####
  dev.print(tiff, paste0('results/', sp, "_", paste((data.frame(lapply(output$combos[best_mod_first_step,], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.tif'),
            width = 10,
            height =8,
            units = "in",
            res = 300)
  
  
  save(results, best_mod_first_step, file = paste0('results/', sp, "_", ifelse(core_type == "CSV", "core_meas", "res"), ".RData"))
  
  write.csv(results$combos[order(results$combos$DeltaAICc),], 
            file = paste0('results/', sp, "_", ifelse(core_type == "CSV", "core_meas_", "res_"), "combos.csv"))
  
  
}



