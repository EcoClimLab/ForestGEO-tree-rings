# ---run climwin and GLS for all sites ---####

# clear environment ####
rm(list = ls())

# load libraries ####
library(climwin)
library(lme4)
library(mgcv)
library(splines)
library(gridExtra)
library(snow)
library(MuMIn)
library(tidyverse)

# prepare parameters ####
## paths to data ####

path_to_climate_data <- "https://raw.githubusercontent.com/forestgeo/Climate/master/Climate_Data/CRU/CRU_v4_04/" # "https://raw.githubusercontent.com/forestgeo/Climate/master/Gridded_Data_Products/Historical%20Climate%20Data/CRU_v4_01/" # 

path_to_climate_data_NM <- "C:/Users/HerrmannV/Dropbox (Smithsonian)/GitHub/EcoClimLab/ForestGEO_dendro/data/cores/NM/CRU_climate/" # *TO BE EDITED* because the dendro repo is private... I ave to give absolute path (tokens are changing all the time in Github with private repos.... so it would be a pain to have to change them everytime...)

path_to_BCI_pre <- "https://raw.githubusercontent.com/forestgeo/Climate/master/Climate_Data/Met_Stations/BCI/El_Claro_precip_starting_1929/pre_BCI.csv"
path_to_BCI_wet <- "https://raw.githubusercontent.com/forestgeo/Climate/master/Climate_Data/Met_Stations/BCI/El_Claro_precip_starting_1929/wet_BCI.csv"

path_to_CO2 <- "https://raw.githubusercontent.com/forestgeo/Climate/master/Other_environmental_data/CO2_data/CO2_MOANA_NOAA_combined.csv"

climate_variables <- c( "pre", "wet",
                        "tmp", "tmn", "tmx", "pet"#,
                        # "dtr", "cld") # "frs", "vap"
                              )


sites.sitenames <- c(BCI = "Barro_Colorado_Island", 
                     CedarBreaks = "Utah_Forest_Dynamics_Plot",
                     HarvardForest = "Harvard_Forest",
                     LillyDickey = "Lilly_Dickey_Woods",
                     SCBI = "Smithsonian_Conservation_Biology_Institute",
                     ScottyCreek = "Scotty_Creek",
                     Zofin = "Zofin",
                     HKK = "Huai_Kha_Khaeng",
                     NewMexico = "New_Mexico")

## analysis parameters ####
reference_dates <- list(BCI = c(30, 12),
                        CedarBreaks = c(30, 8),
                        HarvardForest = c(30, 8),
                        HKK = c(30, 12),
                        LillyDickey = c(30, 8),
                        NewMexico = c(30, 8),
                        SCBI = c(30, 8),
                        ScottyCreek = c(30, 8),
                        Zofin= c(30, 8)) # refday in slidingwin
window_ranges <- list(BCI = c(15, 0),
                      CedarBreaks = c(15, 0),
                      HarvardForest = c(15, 0),
                      HKK = c(15, 0),
                      LillyDickey = c(15, 0),
                      NewMexico = c(15, 0),
                      SCBI = c(15, 0),
                      ScottyCreek = c(15, 0),
                      Zofin = c(15, 0))# range in slidingwin
# reference_date <- c(30, 8) # refday in slidingwin
# window_range <- c(15, 0) #range in slidingwin

clim_var_group <- list(water = c("pre", "wet"),
                       temperature = c("tmp", "tmn", "tmx", "pet")#,
                       # c("dtr", "cld", "pet")
) # see issue 14, PET is in both the TMP and DTR groups. If it comes out as the best in both groups (should always be for the same time frame), then there are only 2 candidate variables for the GLS

clim_gap_threshold <- 5 # 5%

## variable units ####
variables_units <- c(pre = "(mm mo-1)",
                     wet = "(days mo-1)", 
                     tmp = "(C)", 
                     tmn = "(C)", 
                     tmx = "(C)", 
                     pet = "(mm day-1)", 
                     dtr = "(C)", 
                     cld = "%",
                     dbh = "(cm)",
                     CO2 = "(ppm)")
# load data ####
## climate data ####

for(clim_v in climate_variables) {
  assign(clim_v, 
         rbind(
           read.csv(paste0(path_to_climate_data, clim_v,  ".1901.2019-ForestGEO_sites-6-03.csv")), #forestGEO sites
           read.csv(paste0(path_to_climate_data_NM, clim_v, ".1901.2019-NM_site-7-10.csv")) # NM site
           )
  )
  
}

BCI_pre <- read.csv(path_to_BCI_pre, stringsAsFactors = F)
BCI_wet <- read.csv(path_to_BCI_wet, stringsAsFactors = F)

clim_gaps <- read.csv("https://raw.githubusercontent.com/forestgeo/Climate/master/Climate_Data/CRU/scripts/CRU_gaps_analysis/all_sites.reps.csv")
clim_gaps <- clim_gaps[clim_gaps$start_climvar.class %in% climate_variables, ]

## CO2 data ####
CO2 <- read.csv(path_to_CO2)

## core data ####
all_Biol <- read.csv("https://raw.githubusercontent.com/EcoClimLab/ForestGEO_dendro/master/data_processed/all_site_cores.csv?token=AEWDCIMHCL3MDWFTL6JDCEK7I6YUC")

all_Biol <- split(all_Biol, all_Biol$site)

sites <- names(all_Biol)

# prepare data ####

## climate data ####
for(clim_v in climate_variables) {
  print(clim_v)
  x <- get(clim_v)
  
  ### subset for the sites we care about
  x <- x[x$sites.sitename %in% sites.sitenames, ]
  
  ### reshape to long format
  x_long <- reshape(x, 
                    times = names(x)[-1], timevar = "Date",
                    varying = list(names(x)[-1]), direction = "long", v.names = clim_v)
  
  ### format date
  x_long$Date <- gsub("X", "", x_long$Date)
  x_long$Date <- as.Date(x_long$Date , format = "%Y.%m.%d")
  
  ### combine all variables in one
  if(clim_v == climate_variables [1]) all_Clim <- x_long[, c(1:3)]
  else all_Clim <- merge(all_Clim, x_long[, c(1:3)], by = c("sites.sitename", "Date"), all = T)
  
}

### replace BCI pre and wet by local data
idx_BCI <- all_Clim$sites.sitename %in% "Barro_Colorado_Island"

all_Clim$pre[idx_BCI] <- BCI_pre$climvar.val[match(format(all_Clim$Date[idx_BCI], "%Y-%m"), format(as.Date(BCI_pre$Date), "%Y-%m"))]
all_Clim$wet[idx_BCI] <- BCI_wet$climvar.val[match(format(all_Clim$Date[idx_BCI], "%Y-%m"), format(as.Date(BCI_wet$Date), "%Y-%m"))]

### remove BCI pre and wet from clim gap as it is not relevant anymore
clim_gaps <- clim_gaps[!(clim_gaps$start_sites.sitename %in% "Barro_Colorado_Island" & clim_gaps$start_climvar.class %in% c("pre", "wet")), ]

# keep only complete rows (this will remove BCI dat for years where we don't have pre data)

all_Clim <- all_Clim[complete.cases(all_Clim), ]


### format date to dd/mm/yyyy
all_Clim$Date <- format(all_Clim$Date, "%d/%m/%Y") 

### add site column
all_Clim$site <- names(sites.sitenames)[match(all_Clim$sites.sitename, sites.sitenames) ]


### save prepared clim data 
write.csv(all_Clim, "processed_data/Climate_data_all_sites.csv", row.names = F)


## cores ####


for (site in sites) {
  Biol <- all_Biol[[site]]
  Clim <- droplevels(all_Clim[all_Clim$sites.sitename %in% sites.sitenames[site], ])
  
  ### format date to dd/mm/yyyy ####
  Biol$Date <- paste0("15/06/", Biol$Year) # dd/mm/yyyy setting up as june 15, ARBITRATY
  
  ## remove years that are before climate record (+ first few first months to be able to look at window before measurement) ####
  Biol <- Biol[Biol$Year >= min(as.numeric(substr(Clim$Date, 7, 10)))+  window_ranges[[site]][1]/12, ]
  
  ## remove years that are after climate record ####
  Biol <- Biol[Biol$Year <= max(as.numeric(substr(Clim$Date, 7, 10))), ]
  
  
  
  ## remove any cores that have less than 30 years ####
  Biol <- Biol[Biol$coreID %in% names(which(table(Biol$coreID)>= 30)), ]
  
  
  
  ## consider coreID as factor ####
  Biol$coreID <- factor(as.character(Biol$coreID))
  ## save back into Biol ####
  
  all_Biol[[site]] <- Biol
  
  
}

# summarise variable-month combinations gaps ####
adjusted_clim_gaps <- list()
variables_to_drop <- list() # these will automatically be dropped because no matter what time window is best, all momths have > 5% records missing
for(site in sites) {
  
  Biol <- all_Biol[[site]]
  Clim <- droplevels(all_Clim[all_Clim$sites.sitename %in% sites.sitenames[site], ])
  clim_gap <- clim_gaps[clim_gaps$start_sites.sitename %in% sites.sitenames[site], ]
  
  if(nrow(clim_gap) > 0) {
    Clim$Year <- as.numeric(sapply(strsplit(Clim$Date, "/"), "[[", 3))
  # if(min(Biol$Year) != 1903) stop("I did not code for this eventuallity")
  start_year <- min(min(Biol$Year), min(Clim$Year)) # taking the min because technically Clim year starts earlier to be able to look before first measurement
  end_year <- max(Biol$Year) # taking Biol because already trimed to max date.
  n_years <- end_year - start_year
  max_number_months_missing <- floor(n_years * clim_gap_threshold/100)

  # adjust the number of years missing when start and end Dates are before or after (resp) the start and end year of the analysis
  adjusted_clim_gap <- clim_gap
  idx_starts_before <- clim_gap$start_Date < start_year
  idx_ends_after <- clim_gap$end_Date > end_year
  
  adjusted_clim_gap[idx_starts_before, ]$nyears <- adjusted_clim_gap[idx_starts_before,  ]$nyears - (start_year - clim_gap$start_Date[idx_starts_before])
  
  adjusted_clim_gap[idx_ends_after, ]$nyears <- adjusted_clim_gap[idx_ends_after,  ]$nyears -(clim_gap$end_Date[idx_ends_after] - end_year) 
  
  # replace negative values by 0
  adjusted_clim_gap$nyears[adjusted_clim_gap$nyears <0] <- 0 #adjusted_clim_gap <- adjusted_clim_gap[adjusted_clim_gap$nyears > 0, ]
  
  # make month and variable factors so that they are filled with 0 missing years
  adjusted_clim_gap$month <- factor(adjusted_clim_gap$month, levels = 1:12)
  adjusted_clim_gap$start_climvar.class <- factor(adjusted_clim_gap$start_climvar.class)
  
  # sum number of nyears per variable-month combo
  adjusted_clim_gap <- aggregate(formula = nyears ~ start_climvar.class + month,
                                 data = adjusted_clim_gap, 
                                 FUN = sum,
                                 drop = F)
  
  adjusted_clim_gap$nyears[is.na(adjusted_clim_gap$nyears)] <- 0
  
  # get to % missing data
  adjusted_clim_gap$percent_misssing <- adjusted_clim_gap$nyears *100/ n_years
  
  # get if month passes threshold
  adjusted_clim_gap$remove <- adjusted_clim_gap$nyears > max_number_months_missing
  
  # change factors back to characters
  adjusted_clim_gap$month <- as.character(adjusted_clim_gap$month)
  adjusted_clim_gap$start_climvar.class <- as.character(adjusted_clim_gap$start_climvar.class)
  
  # save what variables shold be removed and new adjusted gaps data
  variables_to_drop[[site]] <- names(which(tapply(adjusted_clim_gap$remove, adjusted_clim_gap$start_climvar.class, sum)==12))
  
  adjusted_clim_gaps[[site]] <- adjusted_clim_gap
  
  } 
}


# add CO2 data to all_Biol ####
all_Biol <- lapply(all_Biol, function(Biol) {
  Biol$CO2 <- CO2$CO2_ppm[match(Biol$Year, CO2$year)] 
  return(Biol)
  })

## Run the Analysis ####
variables_dropped <- list() # this is to store the variables that were not conserdered for best model because the average % gap of the window >=5%
summary_data <- NULL # this will hold n_tree, n_cores and range of years
climate_interactions <- NULL # this is just to show case
best_model_summaries <- NULL
# ylim_p <- list() # this is to later readjust ylim in the GLS results plots, to standardize ylim across sites

## save every object names up until now to erase other stuff before runing each new site
data_to_keep <- c(ls(), "data_to_keep")


for(site in sites) {
  
  
  rm(list = ls()[!ls() %in% data_to_keep])
  cat(paste("running analysis for", site, "...\n" ))
  
  
  reference_date <- reference_dates[[site]]
  window_range <-  window_ranges[[site]]
  
  file.remove(list.files("results/explorations/residuals_by_tag/", pattern = site, full.names = T))
  
  variables_dropped[[site]] <- list()
  
  for(what in switch(as.character(any(!is.na(all_Biol[[site]]$dbh))), "TRUE" = c("log_core_measurement", "log_core_measurement_dbh", "log_agb_inc_dbh", "log_BAI_dbh"), "FALSE" = "log_core_measurement")) {
    
    Biol <- all_Biol[[site]]
    Clim <- droplevels(all_Clim[all_Clim$sites.sitename %in% sites.sitenames[site], ])

    ## create folders if don't exist ####
    dir.create(paste0("results/", what, "/", site), showWarnings = F, recursive = T)
    dir.create(paste0("results/with_CO2/", what, "/", site), showWarnings = F, recursive = T)
    
    ## remove all files so that we start clear 
    file.remove(list.files(paste0("results/", what, "/", site), full.names = T))
    file.remove(list.files(paste0("results/with_CO2", what, "/", site), full.names = T))
    
    ## if we use dbh, remove any missing dbh ####
    if(grepl("dbh", what)) Biol <- droplevels(Biol[!is.na(Biol$dbh), ])
    
    ## if we use agb_inc, remove any missing agb_inc ####
    if(grepl("agb_inc", what)) Biol <- droplevels(Biol[!is.na(Biol$agb_inc), ])
    
    ## if we use BAI, remove any missing agb_inc ####
    if(grepl("BAI", what)) Biol <- droplevels(Biol[!is.na(Biol$BAI), ])

    
    ## remove any species that has less than 7 different treeID (even if we use coreID as random effect) ####
    Biol <- droplevels(Biol[Biol$species_code %in% names(which(tapply(Biol$treeID, Biol$species_code, function(x) length(unique(x))) >= 7)), ])
    
    ## calculate residuals of spline measurement ~ year for each individual####
    Biol$residuals <- NA
    
    for( t in unique(Biol$coreID)) {
      x <- Biol[Biol$coreID %in% t, ]
      
      x$Y <- x[, switch(gsub("_dbh", "", what), 
                        log_core_measurement = "core_measurement", 
                        log_agb_inc = "agb_inc",
                        log_BAI = "BAI")]
      
      first_year_removed <- any(is.na(x$Y))
      x <- x[!is.na(x$Y),] #remove NA (only first year of measurement for agb_inc)
      
      test <- gam(Y~ s(Year), data = x)
      
      
      if(rbinom(1, 1, 0.05)==1) {
        png(paste0('results/explorations/residuals_by_tag/', paste(x$species_code[1], x$tree_status[1], t, sep = "_" ), "_", gsub("log_", "", what), "_Year_GAM", "_", site, '.png'),
            width = 8,
            height =8,
            units = "in",
            res = 300)


        par(mfrow = c(3,2))
        plot(test)
        gam.check(test,pch=19,cex=.3)

        plot(Y~ Year, data = x, main = "Raw data")
        points(test$fitted.values ~ x$Year, type = "l")

        title(paste(x$species_code[1], x$tree_status[1], t, sep = " - " ), outer = T, line = -2)

        # save plot
        dev.off()
      }
      
      # save residuals
      x$residuals <- test$residuals 
      
      # remove Y
      x$Y <- NULL
      
      # save back into Biol
      if(gsub("_dbh", "", what) %in% c("log_agb_inc", "log_BAI") & first_year_removed) {
        Biol[Biol$coreID %in% t, ]$residuals <- c(NA, x$residuals) 
      } else {
        Biol[Biol$coreID %in% t, ] <- x
      }
      
      
    }
    
    
    ## run slidingwin on residuals to find best time window with quad for each variable ####
    if(nlevels(Biol$species_code) > 3) {
      baseline = "lmer(residuals ~ 1 + (1 | species_code) + (1 | coreID), data = Biol[!is.na(Biol$residuals), ])"
      
    } else 
    {
      baseline = "lmer(residuals ~ 1 + (1 | coreID), data = Biol[!is.na(Biol$residuals), ])"
    }
    
    
    results <- slidingwin( baseline = eval(parse(text = baseline)),
                           xvar =list(#dtr = Clim$dtr,
                                      pet = Clim$pet, 
                                      tmn = Clim$tmn, 
                                      tmp = Clim$tmp, 
                                      tmx = Clim$tmx,
                                      #cld = Clim$cld, 
                                      pre = Clim$pre, 
                                      wet = Clim$wet
                           ),
                           type = "absolute", 
                           range = window_range,
                           stat = c("mean"),
                           func = "quad", # c("lin","quad")
                           refday = reference_date,
                           cinterval = "month",
                           cdate = Clim$Date, bdate = Biol[!is.na(Biol$residuals), ]$Date,
                           cmissing = "method1") 
    
    # Check if best windows for each variable meet the maximum gap requirement (average of month-variale combination <= 5% of the time period) ####
    adjusted_clim_gap <- adjusted_clim_gaps[[site]]
    variables_dropped_site <- NULL
    
    for(i in 1:nrow(results$combos)) {
      
      months_to_avg <- results$combos$WindowOpen[i]:results$combos$WindowClose[i]
      months_to_avg <- reference_date[2] - months_to_avg
      # months_to_avg <- months_to_avg[months_to_avg!= 0 ]
      months_to_avg <- ifelse(months_to_avg<=0, rev(1:12)[abs(months_to_avg)+1], months_to_avg)
      
      idx_v <- adjusted_clim_gap$start_climvar.class %in% results$combos$climate[i]
      
      if(any(idx_v) & mean(adjusted_clim_gap[idx_v, ][match(months_to_avg, adjusted_clim_gap$month[idx_v]),]$percent_misssing) > clim_gap_threshold ) {
        variables_dropped_site <- c(variables_dropped_site, as.character(results$combos$climate[i]))
      } 
    }
    variables_dropped[[site]][[what]] <- variables_dropped_site
    
    ### find best window for each variable group ignoring the variables that need to be dropped because of gap filling issues (best way to not get confused later with the ordering of te slindingwin order) ####
    clim_var_group_site <- lapply(clim_var_group, function(x) x[!x %in% variables_dropped[[site]][[what]]])
  
    best_results_combos <- do.call(rbind, lapply(clim_var_group_site, function(X) {
      x <- results$combos[results$combos$climate %in% X,]
      data.frame(model_ID = as.numeric(rownames(x)), x, stringsAsFactors = F)[which.min(x$DeltaAICc),]
    }))
    ### plot climwin results ####
    for(i in best_results_combos$model_ID) {
      png(paste0('results/', what, '/', site, '/climwin_', paste((data.frame(lapply(results$combos[i, c(2, 5, 7, 8)], as.character), stringsAsFactors=FALSE)), collapse = "_"), "_", site, '.png'),
          width = 10,
          height =8,
          units = "in",
          res = 300)
      # plot the results
      plotall(dataset = results[[i]]$Dataset, 
              bestmodel = results[[i]]$BestModel,
              bestmodeldata = results[[i]]$BestModelData,
              title=paste((data.frame(lapply(results$combos[i,], as.character), stringsAsFactors=FALSE)), collapse = "_"), arrow = T)
      
      # dev.off()
      dev.off()
      
    }
    
    # # if PET does not come out as top variable in both T and CLD group, drop it (commeted out) ####
    # if(any(duplicated(best_results_combos))) { best_results_combos <- best_results_combos[!duplicated(best_results_combos), ]# remove one pet as it came out best in both T and cld group
    # } else {
    #   best_results_combos <- best_results_combos[!best_results_combos$climate %in% "pet", ] # remove pet as it was bitten in either T or CLD group
    # }
    # best_results_combos <- best_results_combos[order(best_results_combos$DeltaAICc),]
    # 
    ### Save the signal into Biol ####
    for(i in best_results_combos$model_ID) {
      print(paste("adding climate data to Biol for model", i))
     
      
      # if(any(grepl("I\\(climate\\^2\\)", names( results[[i]]$BestModelData)))) {
      #   columns_to_add <- results[[i]]$BestModelData[, c("climate", "I(climate^2)")]
      #   names(columns_to_add) <- paste0(results$combos[i,]$climate, c("", "^2"))
      # } else {
        columns_to_add <- results[[i]]$BestModelData[c("climate")]
        names(columns_to_add) <- results$combos[i,]$climate
      # } if(any(grepl("I\\(climate\\^2\\)", names( results[[i]]$BestModelData)))) 
      
      to_add_to_Biol <-  cbind(Biol[!is.na(Biol$residuals), ], columns_to_add)
      Biol[, names(columns_to_add)] <- to_add_to_Biol[match(rownames(Biol), rownames(to_add_to_Biol)), names(columns_to_add)]
    }
    
    ## Output Biol and best_results_combos to use in different analysis ####
    dir.create(paste0("processed_data/core_data_with_best_climate_signal/", what, "/"), recursive = T, showWarnings = F)
    write.csv(Biol, file = paste0("processed_data/core_data_with_best_climate_signal/", what, "/", site, ".csv"), row.names = F)
  
    
   
    # now do a species by species gls using log of raw measurements, spline on dbh and year (WITH AND WITHOUT CO2) ####
    
    for(with_CO2 in c(FALSE, TRUE)) {
      ## look at collinearity between climate variables ( and CO2n and dbh when relevant) and remove any variable with vif > 10 ####
    if(grepl("dbh", what) & with_CO2) X <- Biol[, c(as.character(best_results_combos$climate), "CO2", "dbh")]
    if(grepl("dbh", what) & !with_CO2) X <- Biol[, c(as.character(best_results_combos$climate), "dbh")]
    if(!grepl("dbh", what) & with_CO2) X <- Biol[, c(as.character(best_results_combos$climate), "CO2")]
    if(!grepl("dbh", what) & !with_CO2) X <- Biol[, c(as.character(best_results_combos$climate))]
    
    X <- X[!duplicated(X),]
    
    
    usdm::vif(X)
    (vif_res <-  usdm::vifstep(X, th = 3))
    sink(paste0("results/", ifelse(with_CO2, "with_CO2/", ""), what, "/", site, "/VIF.txt"))
    print(vif_res)
    sink()
    variables_to_keep <- as.character(vif_res@results$Variables)
    
    
    ## create the gls formula ####
    
    full_model_formula <- switch(gsub("_dbh", "", what), 
                                 "log_core_measurement" = paste("log_core_measurement ~", paste(variables_to_keep, "+", paste0("I(", variables_to_keep, "^2)"), collapse = " + ")),
                                 "log_agb_inc" = paste("log_agb_inc ~", paste(variables_to_keep, "+", paste0("I(", variables_to_keep, "^2)"), collapse = " + ")),
                                 "log_BAI" = paste("log_BAI ~", paste(variables_to_keep, "+", paste0("I(", variables_to_keep, "^2)"), collapse = " + ") ))
    
    # prepare subset for dredge
    sexpr <-parse(text = paste("!(", paste0("(`I(", variables_to_keep, "^2)` & !", variables_to_keep, ")", collapse = "|"), ")"))
    
    # remove second order term for CO2
    full_model_formula <- gsub("I\\(CO2\\^2\\) \\+ ", "", full_model_formula)
    sexpr <- parse(text =gsub("\\| \\(`I\\(CO2\\^2\\)` & \\!CO2\\)", "", sexpr))
   
    
    ## identify what variables we should keep for each species, looking at the sum of AIC weights ####
    
    sum_of_weights_for_each_term_by_sp <- NULL
    for(sp in unique(Biol$species_code)) {
      start_time <- Sys.time()
      x <- Biol[Biol$species_code %in% sp,]
      
      cat("Running GLS and dredging for species", sp , "and its", length(unique(x$coreID)), "trees...\n")
      
      # log transform the response variables
      x$log_core_measurement <- log(x$core_measurement) 
      x$log_agb_inc <-  log(x$agb_inc)
      x$log_BAI <-  log(x$BAI)
      
      # replace the undefined log values by the log of half the smallest value
    x$log_core_measurement[x$core_measurement %in% 0] <- log(min( x$core_measurement[x$core_measurement >0])/2)
    x$log_agb_inc[x$agb_inc %in% 0] <- log(min( x$agb_inc[x$agb_inc >0])/2)
    x$log_BAI[x$BAI %in% 0] <- log(min( x$BAI[x$BAI >0])/2)
      # x <- x[, c("dbh", "Year", "tag", what, variables_to_keep)]
      
      x <- x[, c("Year", "treeID", "coreID", "dbh", "CO2", gsub("_dbh", "", what), variables_to_keep)]
      
      x <- droplevels(x[!is.na(x[, gsub("_dbh", "", what)]), ])
      
      fm1 <- lme((eval(parse(text = full_model_formula))), random = ~1|coreID, correlation = corCAR1(form=~Year|coreID), data = x, na.action = "na.fail", method = "ML") 
      
      clust <- makeCluster(getOption("cl.cores", 2), type = "SOCK")
      clusterEvalQ(clust, library(splines))
      clusterEvalQ(clust, library(nlme))
      clusterExport(clust, list  = c("x", "fm1", "sexpr"))
      dd <- pdredge(fm1, trace = 2, cluster = clust, subset = sexpr)
      stopCluster(clust)
      dd$cw <- cumsum(dd$weight)
      
       # save dd
      write.csv(dd, file = paste0("results/", ifelse(with_CO2, "with_CO2/", ""), what, "/", site, "/", sp, "_model_comparisons.csv"), row.names = F)
      
      # get sum of weights
      # sum_of_weights_for_each_term <- dd[, grepl(paste(c(variables_to_keep, "dbh", "Year"), collapse = "|"), names(dd))]
      sum_of_weights_for_each_term <- dd[, grepl(paste(c(variables_to_keep,  "dbh", "Year", "coreID"), collapse = "|"), names(dd))]
      sum_of_weights_for_each_term <- apply(sum_of_weights_for_each_term, 2, function(x) sum(dd$weight[!is.na(x)]))
      sum_of_weights_for_each_term
      sum_of_weights_for_each_term_by_sp <- rbind(sum_of_weights_for_each_term_by_sp,
                                                  c(sum_of_weights_for_each_term))
      
      
      
      # update the results of the best model to REML = T
      best_model <- update(get.models(dd, 1)[[1]], method = "REML")
      # get the results of the model that includes the variables that have sum of weight > 0.9
      # if(sum(sum_of_weights_for_each_term[!grepl("(coreID|Year)", names(sum_of_weights_for_each_term))] > 0.9) == 0 ) {
      #   best_model <- update(fm1, fixe = ~1, method = "REML") #intercept only model if none of the variables have 0.9 sum of weitgh
      # } else {
      #   best_model <- update(fm1, fixe = paste( "~", paste(names(sum_of_weights_for_each_term)[sum_of_weights_for_each_term > 0.9], collapse = " + ")), method = "REML")
      # }
      
      
      
      # save results for individual species
      assign(paste0(sp, "_dd"), dd)
      assign(paste0(sp, "_best_model"), best_model)
      
      # save into best_model_summaries
      best_model_summaries <- rbind(
        best_model_summaries,
        data.frame(with_CO2,
                   what,
                   site,
                   sp,
                   parameter = rownames(summary(best_model)$tTable),
                   summary(best_model)$tTable)
      )
      
      # keeping track of timing
      end_time <- Sys.time()
      (ellapsed_time <- difftime(end_time, start_time, units = "mins"))
  
      
      # try interaction between DBH and climate variable ####
      if(grepl("dbh", what) & !with_CO2 & "dbh" %in% variables_to_keep) {
        
        for(g in names(clim_var_group)) {
          v_int <-   names(which(sapply(clim_var_group[[g]], grepl, as.character(best_model$call[2]))))
        
          if(length(v_int) > 0 ) {
            bm_clim_int <- update(best_model, fixe = as.formula(paste0("~ . - I(dbh^2) + dbh:", v_int)))
            climate_interactions <- 
              rbind(climate_interactions,                                         data.frame(with_CO2,
                                                                                             what,
                                                                                             site,
                                                                                             species = sp,
                                                                                             climate_group = g,
                                                                                             model_formula = as.character(bm_clim_int$call[2]),
                                                                                             variable = v_int,
                                                                                             interaction_coef = fixef(bm_clim_int)[grepl(":", names(fixef(bm_clim_int)))],
                                                                                             p_value = summary(bm_clim_int)$tTable[grepl(":", rownames(summary(bm_clim_int)$tTable)), "p-value"])
              )
          }
          
          
          } #  for(g in names(clim_var_group)
       
        
        
        
      }
      
      # remove x
      rm(x)
    }
    
    
    rownames(sum_of_weights_for_each_term_by_sp) <- unique(Biol$species_code)
    
    sum_of_weights_for_each_term_by_sp
    
    # save the plot
    png(paste0('results/', ifelse(with_CO2, "with_CO2/", ""), what, "/", site, "/GLS_Sum_of_AICweights_", site, '.png'),
        width = 10,
        height =8,
        units = "in",
        res = 300)
    print(lattice::levelplot(t(sum_of_weights_for_each_term_by_sp), 
                       scales=list(x=list(rot=45)), 
                       xlab = "parameter", 
                       ylab = "species",
                       legend = list(top = list(fun = grid::textGrob("Sum of Weights", y=0, x=1.09)))))
    
    Sys.sleep(time = 1) # to give time for plot to show up
    #dev.off()
    dev.off()
    
    
    
    ## Plot response curves for each variable, with one curve per species ####
    # first remove any object starting by p_
    rm(list = ls()[grepl("^p_", ls())])
    #Create a custom color scale
    ylim_p <- list()
    for(v in c( variables_to_keep)) { 
      print(v)
      
      ## predictions
      pt <- NULL
      for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
        best_model <- get(paste0(sp, "_best_model"))
        x <- Biol[Biol$species_code %in% sp,]
        
        varying_x <- data.frame(varying_x = seq(min(x[, v], na.rm = T), max(x[, v], na.rm = T), length.out = 100)) ; colnames(varying_x) <- v
        # constant_variables <- c("dbh", variables_to_keep)[!c("dbh", variables_to_keep) %in% v]
        constant_variables <- c(variables_to_keep)[!c(variables_to_keep) %in% v]
        
        newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(x$", constant_variables, ", na.rm = T)", collapse = ", "), ")"))), varying_x)  # newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(Biol$", constant_variables, ", na.rm = T)", collapse = ", "), ",  Year = median(Biol$Year, na.rm = T), coreID = factor(Biol[Biol$species_code %in% sp,]$coreID[1]))"))), varying_x)
        
        # if(v %in% names(best_model$var.summary)) {
        # pt <- rbind(pt, data.frame(newd, variable = v, species = sp, varying_x = newd[, v], do.call(rbind, lapply(1:nrow(newd), function(i) data.frame(predict(best_model, newd[i,], type = "link", level = 0, se.fit = T)))), draw = any(grepl(v,names(best_model$coefficients$fixed)))))
       
        
        pt <- rbind(pt, 
                    data.frame(newd, 
                               variable = v, 
                               species = sp, 
                               varying_x = newd[, v], data.frame(predict(best_model, newd, type = "link", level = 0, se.fit = T)), 
                               draw = any(grepl(v,names(best_model$coefficients$fixed))),
                               sig = ifelse(summary(best_model)$tTable[rev(grep(v, rownames(summary(best_model)$tTable)))[1], "p-value"] < .05, "solid", "dashed")))
        
        # }
        
      } # ignore warnings
      
      # if(!is.null(pt)) {
      pt$species <- factor(pt$species, levels = rownames(sum_of_weights_for_each_term_by_sp))
      pt$species <- factor(paste0(pt$species, " (",tapply(Biol$coreID,  Biol$species_code, function(x) length(unique(x)))[as.character(pt$species)], ")"))
      pt$expfit <- exp(pt$fit)
      pt$lwr <- exp(pt$fit - 1.96 * pt$se.fit)
      pt$upr <- exp(pt$fit + 1.96 * pt$se.fit)
     
      
      p <- ggplot(data = pt[pt$draw,], aes(x = varying_x, y = expfit, group = species))
      if(v != "dbh") p <- p + geom_rect(xmin = mean(Biol[, v], na.rm = T) - sd(Biol[, v], na.rm = T), ymin = min(pt$lwr), xmax = mean(Biol[, v], na.rm = T) + sd(Biol[, v], na.rm = T), ymax = max(pt$upr), fill = "grey" , alpha=0.01) + geom_vline(xintercept = mean(Biol[, v], na.rm = T), col = "grey")
      
      time_window <- reference_date[2] - as.numeric(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")])
      time_window_prev <- time_window <= 0
      time_window <- ifelse(time_window_prev, rev(1:12)[abs(time_window) +1], time_window)
      
      time_window_text <- paste(paste0(ifelse(time_window_prev, "p.", "c."),month.abb[time_window]), collapse = "-") #  paste0("\nfrom ", paste(paste0(ifelse(time_window_prev, "prev. ", "curr. "), month.abb[time_window]), collapse = "\nto "))
      
      p <- p + geom_line(aes(col = species), linetype = as.character(p$data$sig)) +
        # scale_x_continuous(trans= ifelse(v %in% "dbh", 'log','identity')) +
        labs(#title = paste0(v, ifelse(v %in% best_results_combos$climate, time_window_text, "")),
             x = paste(v, ifelse(v %in% best_results_combos$climate, time_window_text, ""), variables_units[[v]]),
             y = "") + #"core measurements") +
        geom_ribbon(aes(ymin=lwr, ymax=upr, bg = species), alpha=0.25) +
        scale_colour_hue(drop = F) + scale_fill_hue(drop = F) + 
        theme_classic()
      
      if(any(pt$draw)) {
        assign(paste0("p_", v), p +
               theme(legend.position="none")) 
      }
     
      # save the ylimits ####
      ylim_p[[v]] <- range(pt$lwr, pt$upr)
      # if(!with_CO2) ylim_p[[what]][[v]][[site]] <- range(pt$lwr, pt$upr)
    }
    
    g_legend <- function()
      {
      a.gplot <- ggplot(data = pt, aes(x = varying_x, y = expfit)) + geom_line(aes(group = species, col = species)) +  geom_ribbon(aes(ymin=lwr, ymax=upr, bg = species), alpha=0.25)
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
      }
    
    # existing_plots <- paste0("p_",  c("dbh", variables_to_keep))
    existing_plots <- ls()[grepl("^p_", ls())]
    
    grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  get(x) + ylim(range(ylim_p))), ncol = ifelse(length(existing_plots) %in% 4, 2, length(existing_plots)))),
                 g_legend(),
                 nrow = 1,
                 widths = c(10, 1))
    
    grid::grid.text(switch (gsub("_dbh", "", what), log_core_measurement = "core measurement (mm)",
                            log_agb_inc = "AGB increment (Mg C)", log_BAI = "BAI (cm2)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90)
    
    
    
    # save plot
    png(paste0('results/', ifelse(with_CO2, "with_CO2/", ""), what, "/", site, "/GLS_ALL_variables_responses_", site, '.png'),
        width = 10,
        height =8,
        units = "in",
        res = 300)
    
    grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  get(x) + ylim(range(ylim_p))), ncol = ifelse(length(existing_plots) %in% 4, 2, length(existing_plots)))),
                 g_legend(),
                 nrow = 1,
                 widths = c(10, 1))
    
    grid::grid.text(switch (gsub("_dbh", "", what), log_core_measurement = "core measurement (mm)",
                            log_agb_inc = "AGB increment (Mg C)",
                            log_BAI = "BAI (cm2)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90)
    dev.off()
     
    # save plots at this point to later fetch them  ####
    save(list = grep("^p_|pt|clim_var_group$|ylim_p", ls(), value = T), file = paste0('results/', ifelse(with_CO2, "with_CO2/", ""), what, "/", site, "/env.RData"))
    
   
    } # for(with_CO2 in c(FALSE, TRUE)) 
     
   # summarize data used in this analysis ####
    summary_data <- rbind(summary_data, data.frame(site = site, what = what,  
                          Biol %>% group_by(species_code) %>%
                            summarize(
                                      n_trees = n_distinct(treeID),
                                      n_cores = n_distinct(coreID),
                                      start_year = min(Year),
                                      end_year = max(Year)                                    )))
    
    
  } # for (what in ...)
  
  
  # save environment ####
  save.image(file = paste0("results/", site, "_all_env.RData"))
} # for sites in ..

# save summary_data ####
write.csv(summary_data, "results/summary_cores_analyzed.csv", row.names = F)

# save best_models summaries ####
write.csv(best_model_summaries, file = "results/best_model_summaries.csv", row.names = F)
# save and summarize climate_interactions ####
write.csv(climate_interactions, file = "results/climate_interactions_coeficients.csv", row.names = F)

climate_interactions_summary <- aggregate(p_value ~ what + site + climate_group, data = climate_interactions, FUN = function(x) round(sum(x<0.05)*100/length(x), 2))
names(climate_interactions_summary) <- gsub("p_value", "freq_sig", names(climate_interactions_summary))

write.csv(climate_interactions_summary, file = "results/climate_interactions_summary.csv", row.names = F)

# save variables_dropped ####
variables_dropped

sink("results/variables_dropped_at_each_site_due_to_gap.txt")
variables_dropped
sink()

# save(ylim_p, file = paste0("results/ylims_for_GLS_plots.RData"))


# see VIF situations for all sites ####
VIF_files <- list.files("results", pattern = "VIF.txt", recursive = T, full.names = T)
VIF_files <- VIF_files[!grepl("CO2", VIF_files)]

# check if any variables were kicked because of linearity issues
for(f in VIF_files) {
  x <- readLines(f)[1]
  if(substr(x, 1, 2) != "No") print(f)
  # cat(f, "\n")
  # cat(x, "\n\n")
} # if nothng shows up, no variables were kicked because of collinarity issues

# look at VIF (if any are close to 3)
for(f in VIF_files) {
  x <- readLines(f)
  x <- x[7:length(x)]
  # if(substr(x, 1, 2) != "No") print(f)
  cat("\n\n", f, "\n")
   print(x)
} # if nothng shows up, no variables were kicked because of collinarity issues


# summary plot for climate  ####

A <- pivot_longer(all_Clim, climate_variables, "climate_var")
A$Year <- as.numeric(gsub("\\d\\d/\\d\\d/", "", A$Date))

A <- aggregate(value ~ sites.sitename + Year + climate_var, data = A, FUN =mean)
A <- rbind(A, data.frame(sites.sitename = "All",
                         Year = CO2$year,
                         climate_var = "CO2",
                         value = CO2$CO2_ppm))
variables_units_full <- variables_units
variables_units_full[] <- paste(names(variables_units), variables_units, sep = "\n")

png("results/Climate_variables_yearly_mean.png", width = 8, height = 10, res = 300, units = "in")
ggplot(A, aes(x = Year, y = value, color = sites.sitename)) +
  geom_line() + 
  facet_wrap(vars(climate_var), ncol =1, scales = "free_y", strip.position = "left", labeller = as_labeller(variables_units_full)) +
  ylab(NULL) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.placement = "outside")
 
dev.off()

