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
library(splines)
library(grid)

# set parameters ####
core_type <- "CSV"
reference_date <- c(30, 8) # refday in slidingwin
window_range <- c(15, 0) #range in slidingwin


# create function ####
calculate_bark_thickness_ln <- function(dbh, sp){
  # this is based off of allometric equations developped by Ian mackgregor (see https://github.com/SCBI-ForestGEO/McGregor_climate-sensitivity-variation/blob/master/manuscript/tables_figures/tableS1_bark_regression.csv) from data published here: https://datadryad.org/stash/dataset/doi:10.5061/dryad.6nc8c
 
  sp <- as.character(sp)
  dbh <- as.numeric(dbh)
  if(!sp %in% c("caco", "cagl", 
                "caovl", "cato", "fram", "juni", "litu", "qual", "qupr", "quru", 
                "quve", "frni", "fagr", "pist")) stop('sp has to be one of c("caco", "cagl", "caovl", "cato", "fram", "juni", "litu", "qual", "qupr", "quru", "quve", "frni", "fagr", "pist")')
  
 warning("- Allometries of Fraxinus americanus were used for Fraxinus nigra.\n  - Bark thickness of Fagus grandifolia is considered to be 0.\n  - Allometries for Pinus strobus are not developped yet (update this code when they are) so Allometry of Quercus rubra is used for now...")
  
  if(!is.na(dbh) & dbh < 0) stop(paste("dbh = ", dbh, "- dbh has to be a positive value"))
  
  switch (sp,
          "caco" = -1.56+0.416*log(dbh),
          "cagl" = -0.393+0.268*log(dbh),
          "caovl" = -2.18+0.651*log(dbh),
          "cato" = -0.477+0.301*log(dbh),
          "fram" = 0.418+0.268*log(dbh),
          "juni" = 0.346+0.279*log(dbh),
          "litu" = -1.14+0.463*log(dbh),
          "qual" = -2.09+0.637*log(dbh),
          "qupr" = -1.31+0.528*log(dbh),
          "quru" = -0.593+0.292*log(dbh),
          "quve" = 0.245+0.219*log(dbh),
          "frni" = 0.418+0.268*log(dbh), # used fram equation
          "fagr" = -Inf, # considered non existent (exp(-Inf) = 0)
          "pist" = -0.593+0.292*log(dbh), # use quru equation for now
          NA
          
  )
  }

# load data ####
## core data ####

cores <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_cores/cross-dated_cores_CSVformat/all_core_chronologies.csv")

dbh_2008 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem1.csv", stringsAsFactors = FALSE)

dbh_2013 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem2.csv", stringsAsFactors = FALSE)

dbh_2018 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem3.csv", stringsAsFactors = FALSE)

## bark data ####

# bark <- read.csv("data/traits/SCBI_bark_depth.csv")
bark <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/McGregor_climate-sensitivity-variation/master/data/traits/SCBI_bark_depth.csv?token=AEWDCIOKCNXPQAECCRNPP7S6MEGGQ")

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


## bark data (do what Ian did in https://github.com/SCBI-ForestGEO/McGregor_climate-sensitivity-variation ####

### keep only species we care about
bark <- droplevels(bark[bark$species %in% unique(cores$sp), ])

### Calculate diameter_nobark for 2008 = DBH.mm.2008-2*bark.depth.mm
bark$diam_nobark_2008.mm <- bark$DBH.mm.2008 - 2*bark$bark.depth.mm 

### predict bark thickness with both dbh and dbh_no_bark and see if that makes a big difference
bark$predict_barkthick_from_no_bark <- apply(bark, 1, function(x) exp(calculate_bark_thickness_ln(dbh = x[["diam_nobark_2008.mm"]], sp = x[["species"]])))
bark$predict_barkthick_from_dbh <- apply(bark, 1, function(x) exp(calculate_bark_thickness_ln(dbh = x[["DBH.mm.2008"]], sp = x[["species"]])))
bark$predicit_dbh_from_no_bark <- bark$diam_nobark_2008.mm + 2*bark$predict_barkthick_from_no_bark 
### Take exponent of bark.depth.mm and make sure predicted values look good.
plot(bark$predict_barkthick_from_no_bark, bark$predict_barkthick_from_dbh)
abline(0,1)
summary(bark$predict_barkthick_from_no_bark - bark$predict_barkthick_from_dbh)
plot(c(bark$predict_barkthick_from_no_bark - bark$predict_barkthick_from_dbh) ~ bark$DBH.mm.2008)
abline(lm(c(bark$predict_barkthick_from_no_bark - bark$predict_barkthick_from_dbh) ~ bark$DBH.mm.2008))
summary(lm(c(bark$predict_barkthick_from_no_bark - bark$predict_barkthick_from_dbh) ~ bark$DBH.mm.2008)) # always underestimtes the bark thickness but no signicicant trend with dbh...
summary(lm(2*c(bark$predict_barkthick_from_no_bark - bark$predict_barkthick_from_dbh) ~ bark$DBH.mm.2008)) # always underestimtes the bark thickness but no signicicant trend with dbh...


plot(c(bark$predicit_dbh_from_no_bark ) ~ bark$DBH.mm.2008)
summary(lm((c(bark$predicit_dbh_from_no_bark ) ~ bark$DBH.mm.2008)))

abline(lm(c(bark$predicit_dbh_from_no_bark ) ~ bark$DBH.mm.2008))
abline(0,1)
table(bark[(bark$predict_barkthick.mm/bark$bark.depth.mm) >1.2,]$species) # evenly overestimated or underestimated
table(bark[(bark$predict_barkthick.mm/bark$bark.depth.mm) <0.8,]$species)# evenly overestimated or underestimated

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

### remove NAs
Biol <- Biol[!is.na(Biol$core_measurement), ]


### cut back to 2000 when tree was cored dead 

for( t in unique(Biol[Biol$status.at.coring == "dead", ]$tag)) {
  x <- Biol[Biol$tag %in% t,]
  
  rows_to_remove <- rownames(x)[x$Year > 2000]
  Biol <- Biol[!rownames(Biol) %in% rows_to_remove, ]
  
}

### find out cores with outliers
tag_with_outliers <- unique(Biol[Biol$core_measurement>10, ]$tag)

tag_with_outliers_within_first_X_years <- NULL
first_years_oulier_limit = 15

par(mfrow = c(4,5), mar = c(2,1,3,1), oma = c(2,3,0,0))
for( t in tag_with_outliers) {
  x <- Biol[Biol$tag %in% t,]
  
  # is the outlier within the first X years?
  within_first_ten_years <- any(which(x$core_measurement > 10) < first_years_oulier_limit)
  
  if(within_first_ten_years) tag_with_outliers_within_first_X_years <- c(tag_with_outliers_within_first_X_years, t)
  
  plot(core_measurement~Year, data = x, type = "l", main = paste(t, x$status.at.coring[1]), xlab = "", col = ifelse(max(x$core_measurement) > 15, "red", "black"))
  abline(h = 10, lty = 2)
  abline(v = sort(x$Year)[first_years_oulier_limit], lty = 2)
  if(within_first_ten_years) {
    arrows(x0 = par("usr")[c(1,2)],
           y0 = par("usr")[c(3,3)],
           x1 = par("usr")[c(2,1)],
           y1 = par("usr")[c(4,4)],code = 0,
           lwd = 2,
              )
  }
  mtext(side = 2, "Core measurement", outer = T, line = 1.5)
  mtext(side = 1, "Year", outer = T)
}

### remove the first X year of tag_with_outliers_within_first_15_years
for( t in tag_with_outliers_within_first_X_years) {
  x <- Biol[Biol$tag %in% t,]
  
  rows_to_remove <- rownames(x)[x$Year <= min (x$Year) + first_years_oulier_limit]
  Biol <- Biol[!rownames(Biol) %in% rows_to_remove, ]
  
}

## find out DBH for each year ####
Biol$dbh = NA
Biol$dbh_no_bark <- NA
Biol$bark_thickness <- NA

tags_with_dbh_issues <- NULL
for( t in unique(Biol$tag)) {
  # print(t)
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
  

  # find out closest year DBH measured
  closest.year.tree.was.cenused.alive <- names(which(y_status == "A"))[which.min(abs(as.numeric(names(which(y_status == "A"))) - max(x$Year)))]
  
  ## if that year exist in the core measurement, just use it as "anchor dbh"
  if(closest.year.tree.was.cenused.alive %in% x$Year)  dbh <- as.numeric(y_dbh[closest.year.tree.was.cenused.alive])
  
  ## if that year does not exist in the core measurement (cored while alive and in a previous census), calculate average ring increment in the past few years and use that to estimate the dbh on the last year of core
  
  if(!closest.year.tree.was.cenused.alive %in% x$Year) { 
    nb_years_to_average_from <- as.numeric(closest.year.tree.was.cenused.alive) - max(x$Year)
    avg_core_increment <- mean(rev(x$core_measurement)[1:nb_years_to_average_from])
    
    
    dbh <- as.numeric(y_dbh[closest.year.tree.was.cenused.alive]) - 2*avg_core_increment*nb_years_to_average_from
    # closest.year.tree.was.cenused.alive <- max(x$Year)
   
    # x[x$Year == closest.year.tree.was.cenused.alive, ]$dbh <- as.numeric(y_dbh[closest.year.tree.was.cenused.alive]) - 2*avg_core_increment*nb_years_to_average_from
    

  }  
  
  # calculate bark thickness and dbh_no_bark
  bark_thickness_mm <- exp(calculate_bark_thickness_ln(dbh, sp = unique(x$sp)))
  dbh_no_bark <- dbh - 2* bark_thickness_mm

  # find out what core measurement were before and after dbh measurements
  idx_x_before <- x$Year <= closest.year.tree.was.cenused.alive 
  idx_x_after <- x$Year > closest.year.tree.was.cenused.alive
  
  # find out if dbh_no_bark a t0 is <0
  min_dbh_no_bark <- dbh_no_bark  - 2*sum(x[idx_x_before, ]$core_measurement)
  
  # if dbh_no_bark <0, redistribute the negatives proportionaly to all core_measurements
  if(min_dbh_no_bark < 0) {
    tags_with_dbh_issues <- c(tags_with_dbh_issues, t)
    
    total_amount_to_correct <- min_dbh_no_bark / 2
    proportion_each_tree_ring <- x$core_measurement / sum(x$core_measurement)
    
    new_core_measurement <- x$core_measurement + (total_amount_to_correct * proportion_each_tree_ring)
    
  } else {
    new_core_measurement <- x$core_measurement
  }
  
  # populate dbh_no_bark before measurement
  x[idx_x_before, ]$dbh_no_bark <-  rev(dbh_no_bark - c(0, cumsum(  rev(2* new_core_measurement[idx_x_before]))))[-1]
  
  
  # populate dbh_no_bark after measurement
  x[idx_x_after, ]$dbh_no_bark <-  dbh_no_bark + cumsum(  2*new_core_measurement[idx_x_after])
  
  # double check diff dbh corresponds to core measurements
  if(!all(round(diff(x$dbh_no_bark),3) - round(2*new_core_measurement[-1], 3) == 0)) stop("dbh calculated invcorrectly")
  
  # caclcuate bark thickness with dbh_no_bark (this will underestimate it but for all years so we are fine I think)
  x$bark_thickness <- exp(calculate_bark_thickness_ln(dbh = x$dbh_no_bark, sp = x$sp[1]))
  
  x$dbh <- x$dbh_no_bark + 2 * x$bark_thickness
  
  # save into Biol
  
  Biol[Biol$tag %in% t , ] <- x
}
  

par(mfrow = c(4,4), mar = c(3,1,3,3), oma = c(1,3,0,1))
for( t in tags_with_dbh_issues) {
  x <- Biol[Biol$tag %in% t,]
  
  plot(core_measurement~Year, data = x, type = "l", main = paste(t, x$status.at.coring[1]), xlab = "", col = ifelse(max(x$core_measurement) > 15, "red", "black"))
  par(new = T)
  plot(dbh~Year, data = x, yaxt = "n", col = "red")
  axis(4, col = "red")
  abline(h = 0, lty = 2, col = "red")
  abline(v = sort(x$Year)[first_years_oulier_limit], lty = 2)
  mtext(side = 2, "Core measurement", outer = T, line = 1.5)
  mtext(side = 4, "dbh", outer = T, line = 0, col = "red")
  mtext(side = 1, "Year", outer = T)
}

## calculate Biomass using equations developped by Erika for now ####
x <- Biol

x$agb <- NA

x$agb <- ifelse(x$sp == "caco", 10^(-1.326 + 2.762 * log10(x$dbh * 0.1)) *1.005, x$agb)# new equation added by Erika 1/13/2020                                          
x$agb <- ifelse(x$sp == "cagl", 10^(-1.326 + 2.762 * log10(x$dbh * 0.1)) *1.005, x$agb)# new equation added by Erika 1/13/2020  
x$agb <- ifelse(x$sp == "caovl", 10^(-1.326 + 2.762 * log10(x$dbh * 0.1)) *1.005, x$agb)# new equation added by Erika 1/13/2020                                              
x$agb <- ifelse(x$sp == "cato", 10^(-1.326 + 2.762 * log10(x$dbh * 0.1)) *1.005, x$agb)# new equation added by Erika 1/13/2020                                              
x$agb <- ifelse(x$sp == "fagr", 10^(2.1112 + 2.462 * log10(x$dbh * 0.1)) / 1000, x$agb)# new equation added by Erika 1/13/2020
x$agb <- ifelse(x$sp == "fram", (2.3626 * (x$dbh * 0.03937)^2.4798) * 0.45359, x$agb)
x$agb <- ifelse(x$sp == "frni", 0.1634 * (x$dbh * 0.1)^2.348, x$agb)
x$agb <- ifelse(x$sp == "juni", exp(-2.5095 + 2.5437 * log(x$dbh * 0.1)), x$agb)
x$agb <- ifelse(x$sp == "litu", (10^(-1.236 + 2.635 * (log10(x$dbh * 0.1)))) * 1.008, x$agb) # new equation given by Erika on Tue 4/2/2019 11:57
x$agb <- ifelse(x$sp == "pist", (exp(5.2831 + 2.0369 * log(x$dbh * 0.1))) / 1000, x$agb)
x$agb <- ifelse(x$sp == "qual", (1.5647 * (x$dbh * 0.03937)^2.6887) * 0.45359, x$agb)
x$agb <- ifelse(x$sp == "qupr", (1.5509 * (x$dbh * 0.03937)^2.7276) * 0.45359, x$agb)
x$agb <- ifelse(x$sp == "quru", (2.4601 * (x$dbh * 0.03937)^2.4572) * 0.45359, x$agb)
x$agb <- ifelse(x$sp == "quve" & (x$dbh * 0.1) < 30, exp(-0.34052 + 2.65803 * log(x$dbh * 0.03937)), x$agb)
x$agb <- ifelse(x$sp == "quve" & (x$dbh * 0.1) >= 30, (10^(1.00005 + 2.10621 * (log10(x$dbh * 0.03937)))) * 0.45359, x$agb)


#Convert from kg to Mg
x$agb <- x$agb / 1000 

# Convert to C
x$agb <-  x$agb * .47
  
Biol <- x

## calculate agb increment for each individual ####
Biol$agb_inc <- NA
for( t in unique(Biol$coreID)) {
  x <- Biol[Biol$coreID %in% t, ]
  x$agb_inc <- c(NA, diff(x$agb))
  
  Biol[Biol$coreID %in% t, ] <- x
}

## remove years that are before climate record (+ first few first months to be able to look at window before measurement) ####
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) >= min(as.numeric(substr(Clim$Date, 7, 10)))+  window_range[1]/12, ]

## remove years that are after climate record ####
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) <= max(as.numeric(substr(Clim$Date, 7, 10))), ]


## calculate residuals of spine measurement ~ year for each individual####

for(what in c("log_core_measurement", "log_agb_inc")) {
  
Biol$residuals <- NA

for( t in unique(Biol$coreID)) {
  x <- Biol[Biol$coreID %in% t, ]

  x$Y <- x[, switch(what, log_core_measurement = "core_measurement", log_agb_inc = "agb_inc")]
  x <- x[!is.na(x$Y),] #remove NA (only first year of measurement for agb_inc)
  
  test <- gam(Y~ s(Year), data = x)
  par(mfrow = c(3,2))
  plot(test)
  gam.check(test,pch=19,cex=.3)
  
  plot(Y~ Year, data = x, main = "Raw data")
  points(test$fitted.values ~ x$Year, type = "l")
  
  title(paste(x$sp[1], x$status.at.coring[1], t, sep = " - " ), outer = T, line = -2)
  
  # save plot
  if(rbinom(1, 1, 0.1)==1) {
    dev.print(tiff, paste0('results/explorations/residuals_by_tag/', paste(x$sp[1], x$status.at.coring[1], t, sep = "_" ), "_", gsub("log_", "", what), "_Year_GAM", '.tif'),
              width = 8,
              height =8,
              units = "in",
              res = 300)
    
  }

  # save residuals
  x$residuals <- test$residuals 
  
  # save back into Biol
  if(what %in% "log_agb_inc") {
    Biol[Biol$coreID %in% t, ]$residuals <- c(NA, x$residuals) 
  } else {
    Biol[Biol$coreID %in% t, ] <- x
  }
 
  
}

## run slidingwin on residuals to find best time window and lin or quad for each variable ####
baseline = "lmer(residuals ~ 1 + (1 | sp) + (1 | coreID), data = Biol)"

results <- slidingwin( baseline = eval(parse(text = baseline)),
                       xvar =list(dtr = Clim$dtr,
                                  pet = Clim$pet, 
                                  tmn = Clim$tmn, 
                                  tmp = Clim$tmp, 
                                  tmx = Clim$tmx,
                                  cld = Clim$cld, 
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

### find best statistique and window for each variable
results$combos

best_results_combos <- do.call(rbind, by(results$combos, results$combos$climate, function(x) data.frame(model_ID = as.numeric(rownames(x)), x, stringsAsFactors = F)[which.min(x$DeltaAICc),]))

best_results_combos <- best_results_combos[order(best_results_combos$DeltaAICc),]

### plot the results and save the signal into Biol ####
for(i in best_results_combos$model_ID) {
  
  # plot the results
  plotall(dataset = results[[i]]$Dataset, 
          bestmodel = results[[i]]$BestModel,
          bestmodeldata = results[[i]]$BestModelData,
          title=paste((data.frame(lapply(results$combos[i,], as.character), stringsAsFactors=FALSE)), collapse = "_"))
  # save the plot
  dev.print(tiff, paste0('results/ALL_species_mixed_model_on_residuals/ALL_species_mixed_model_on_', gsub("log_", "", what), "_", paste((data.frame(lapply(results$combos[i,], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.tif'),
            width = 10,
            height =8,
            units = "in",
            res = 300)
  
  # save the climate signal in Biol
  if(any(grepl("I\\(climate\\^2\\)", names( results[[i]]$BestModelData)))) {
    columns_to_add <- results[[i]]$BestModelData[, c("climate", "I(climate^2)")]
    names(columns_to_add) <- paste0(results$combos[i,]$climate, c("", "^2"))
  } else {
    columns_to_add <- results[[i]]$BestModelData[, c("climate")]
    names(columns_to_add) <- results$combos[i,]$climate
  }
  
  Biol <- cbind(Biol, columns_to_add)
}

names(Biol)

# look at collinearity between climate variables ####
X <- Biol[, c("pre", "wet", "pet", "cld", "dtr", "tmx", "tmp", "tmn")]
X <- X[!duplicated(X),]

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(X, upper.panel = panel.cor)

library(usdm)
vif(X)
(vif_res <- vifstep(X, th = 3))
variables_to_keep <- as.character(vif_res@results$Variables)
pairs(X[, variables_to_keep], upper.panel = panel.cor)

# After discussion, we like to get rid of collinearity by picking what variables made more sense biologically to us so here is the set we move on with:
variables_to_keep <- c("pre", "wet", "cld", "tmx", "tmn")

vifstep(Biol[, variables_to_keep], th = 3) #--> all good

# now do a species by species gam using log of raw measuremets, spline on dbh and year ####
library(MuMIn)


  # create tge gam formula
  
   full_model_formula <- switch(what, "log_core_measurement" =  paste("log_core_measurement ~ s(dbh, k = 3) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")),
                               log_agb_inc = paste("log_agb_inc ~ s(dbh, k = 3) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")))
  
  
  ## identify what variables we should keep for each species, looking at the sum of AIC weights ####
  start_time <- Sys.time()
  
  sum_of_weights_for_each_term_by_sp <- NULL
  for(sp in unique(Biol$sp)) {
    print(sp)
    x <- Biol[Biol$sp %in% sp,]
    x$tag <- factor(x$tag)
    x$log_core_measurement <- log(x$core_measurement+0.1)
    x$log_agb_inc <-  log(x$agb_inc + 0.1)
    x <- x[, c("dbh", "Year", "tag", what, variables_to_keep)]
    
    x <- x[!is.na(x[, what]), ]
    fm1 <- gam(eval(parse(text = full_model_formula)), data = x, na.action = "na.fail") # using mixed model makes all variables important which I find suspicous. I feel that s(Year, by = tag) is enough to add some randomness by tag... + plus I am not even sure that actuallydoes what we want.
    dd <- dredge(fm1)
    dd$cw <- cumsum(dd$weight)
    
    sum_of_weights_for_each_term <- dd[, grepl(paste(c(variables_to_keep, "dbh", "Year"), collapse = "|"), names(dd))]
    sum_of_weights_for_each_term <- apply(sum_of_weights_for_each_term, 2, function(x) sum(dd$weight[!is.na(x)]))
    sum_of_weights_for_each_term
    sum_of_weights_for_each_term_by_sp <- rbind(sum_of_weights_for_each_term_by_sp,
                                                c(sum_of_weights_for_each_term))
    
    
    
    # get the results of the model that includes the variables that have sum of weight > 0.9
    best_model <- gam(eval(parse(text = paste( what,  "~", paste(names(sum_of_weights_for_each_term)[sum_of_weights_for_each_term > 0.9], collapse = " + ")))), data = x, na.action = "na.fail")
    
    
    # save results for individual species
    assign(paste0(sp, "_dd"), dd)
    assign(paste0(sp, "_best_model"), best_model)
    
    # remove x
    rm(x)
  }
  end_time <- Sys.time()
  
  
  (ellapsed_time <- difftime(end_time, start_time))
  
  rownames(sum_of_weights_for_each_term_by_sp) <- unique(Biol$sp)
  
  sum_of_weights_for_each_term_by_sp
  library(lattice)
  
  levelplot(t(sum_of_weights_for_each_term_by_sp), 
            scales=list(x=list(rot=45)), 
            xlab = "parameter", 
            ylab = "species",
            legend = list(top = list(fun = grid::textGrob("Sum of Weights", y=0, x=1.09))))
  
  # save the plot
  dev.print(tiff, paste0('results/Species_by_species_GAMS_on_raw_data/Sum_of_AICweights_', what, '.tif'),
            width = 10,
            height =8,
            units = "in",
            res = 300)
  
  
  
  
  # double check the best models we have by species is correct
  
  best_models_by_species <- apply(sum_of_weights_for_each_term_by_sp, 1, function(x)paste(names(x)[x > 0.9], collapse = " + "))
  
  formula(caco_best_model) == as.formula(paste(what, "~", best_models_by_species[["caco"]]))
  
  
  ## second, plot response curves for each species and variables ####
  
  for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
    
    print(sp)
    
    X <- Biol[Biol$sp %in% sp, ]
    best_model <- get(paste0(sp, "_best_model"))
    
    variables_to_look_at <- names((best_model$var.summary))[!names((best_model$var.summary)) %in% c("Year", "tag")]
    
    n_row <- ifelse(length(variables_to_look_at) <= 2, 1, 2)
    n_col <- length(variables_to_look_at) %/% 2 +  length(variables_to_look_at)%% 2
    
    if(length(variables_to_look_at) > 0) {
      par(mfrow=c(ifelse(n_row == 0, 1, n_row), n_col))
      
      
      for(v in variables_to_look_at) {
        
        X$varying_v <- X[, v]
        
        varying_x <- data.frame(floor(min(X$varying_v)): ceiling(max(X$varying_v))) 
        colnames(varying_x) <- v
        
        X$Y <- X[, switch(what, log_core_measurement = "core_measurement", log_agb_inc = "agb_inc")]
        plot(Y+0.1 ~ varying_v, data = X, log = "y", 
             pch = 16,
             # bg = rgb(0,0,0,0.2),
             col = rainbow(length(unique(X$tag)), s = 0.8, alpha = 0.2)[c(1:length(unique(X$tag)))[match(X$tag, unique(X$tag))]],
             main = paste0(sp[1], " - ", v, ifelse(v %in% best_results_combos$climate, paste0("\nfrom ",
                                                                                              paste(month.abb[reference_date[2] - as.numeric(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")])], collapse = " to ")), "")),
             xlab = v,
             ylab = switch(what, log_core_measurement = "core measurement (mm)", log_agb_inc = "AGB increment (Mg C)"),
             border = "grey")
        
        if(length(variables_to_look_at) > 1) {
          constant_variables <- variables_to_look_at[!variables_to_look_at %in% v]
      
          newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(X$", constant_variables, ")", collapse = ", "), ",  Year = median(X$Year), tag = factor(X$tag[1]))"))), varying_x)
          
        } else {
          newd <- data.frame(varying_x,  Year = median(X$Year), tag = factor(X$tag[1]))
        }
        
        ## prediction
        pt <- predict.gam(best_model, newd, type = "response", exclude =grep("Year", sapply(best_model$smooth, "[[", "label"), value = T), se.fit = T)
        
        ## add preditive line
        lines(exp(pt$fit) ~ varying_x[[1]], lwd = 2)
        lines(exp(pt$fit - 1.96 * pt$se.fit) ~ varying_x[[1]], lwd = 1, lty = 2)
        lines(exp(pt$fit + 1.96 * pt$se.fit) ~ varying_x[[1]], lwd = 1, lty = 2)
        
      }
      
    }
    
    
    
    # save plot
    dev.print(tiff, paste0('results/Species_by_species_GAMS_on_raw_data/GAM_results_raw_', sp, "_", what, ".tif"),
              width = 8,
              height =8,
              units = "in",
              res = 300)
    
  }
  
  
  ## third,  plot response curves for each variable, with one curve per species ####
  #Create a custom color scale
  
  for(v in c("dbh", variables_to_keep)) {
    print(v)
    
    ## predictions
    pt <- NULL
    for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
      best_model <- get(paste0(sp, "_best_model"))
      
      
      varying_x <- data.frame(varying_x = seq(min(Biol[Biol$sp %in% sp, v]), max(Biol[Biol$sp %in% sp, v]), length.out = 100)) ; colnames(varying_x) <- v
      constant_variables <- c("dbh", variables_to_keep)[!c("dbh", variables_to_keep) %in% v]
      
      newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(Biol$", constant_variables, ")", collapse = ", "), ",  Year = median(Biol$Year), tag = factor(Biol[Biol$sp %in% sp,]$tag[1]))"))), varying_x)
      
      if(v %in% names(best_model$var.summary)) {
        pt <- rbind(pt, data.frame(newd, variable = v, species = sp, varying_x = newd[, v], predict.gam(best_model, newd, type = "response", exclude =grep("Year", sapply(best_model$smooth, "[[", "label"), value = T), se.fit = T)))
        
      }
      
    } # ignore errors
    pt$species <- factor(pt$species, levels = rownames(sum_of_weights_for_each_term_by_sp))
    pt$expfit <- exp(pt$fit)
    pt$lwr <- exp(pt$fit - 1.96 * pt$se.fit)
    pt$upr <- exp(pt$fit + 1.96 * pt$se.fit)
    
    p <- ggplot(data = pt, aes(x = varying_x, y = expfit))
    if(v != "dbh") p <- p + geom_rect(xmin = mean(Biol[, v]) - sd(Biol[, v]), ymin = min(pt$lwr), xmax = mean(Biol[, v]) + sd(Biol[, v]), ymax = max(pt$upr), fill = "grey" , alpha=0.01) + geom_vline(xintercept = mean(Biol[, v]), col = "grey")
    
    p <- p + geom_line(aes(group = species, col = species)) +
      # scale_x_continuous(trans= ifelse(v %in% "dbh", 'log','identity')) +
      labs(title = paste0(v, ifelse(v %in% best_results_combos$climate, paste0("\nfrom ",
                                                                               paste(month.abb[reference_date[2] - as.numeric(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")])], collapse = " to ")), "")),
           x = v,
           y = "") + #"core measurements") +
      geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = species), alpha=0.25) + 
      scale_colour_hue(drop = F) + scale_fill_hue(drop = F) + 
      theme_classic()
    
    assign(paste0("p_", v), p +
             theme(legend.position="none"))
  }
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  
  
  grid.arrange(do.call(arrangeGrob, c(lapply(paste0("p_",  c("dbh", variables_to_keep)), function(x)  get(x)), ncol = 3)),
               g_legend(p),
               nrow = 1,
               widths = c(10, 1))
  
  grid.text(switch (what, log_core_measurement = "core measurement (mm)",
                    log_agb_inc = "AGB increment (Mg C)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90)
  
  
  
  # save plot
  dev.print(tiff, paste0('results/Species_by_species_GAMS_on_raw_data/ALL_variables_', what, '.tif'),
            width = 8,
            height =8,
            units = "in",
            res = 300)
  
}



# save environment ####
save.image(file = "Analysis_workspace.RData")
 