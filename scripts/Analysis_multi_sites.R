# ---runing analysis without dbh ---####

# clear environment ####
rm(list = ls())

# load libraries ####
library(climwin)
library(lme4)
library(mgcv)
library(splines)
library(gridExtra)
library(grids)
library(MuMIn)
library(allodb) # remotes::install_github("forestgeo/allodb")

# prepare parameters ####
## paths to data ####
path_to_core_measurement_files <- "C:/users/herrmannv/Dropbox (Smithsonian)/GitHub/EcoClimLab/ForestGEO_dendro/data_processed/"

sites <- c( "BCI", 
            "HarvardForest", "LillyDickey", "SCBI", 
            "ScottyCreek") # "CedarBreaks", 

path_to_climate_data <-"https://raw.githubusercontent.com/forestgeo/Climate/master/Gridded_Data_Products/Historical%20Climate%20Data/CRU_v4_01/"
climate_variables <- c( "pre", "wet",
                        "tmp", "tmn", "tmx", "pet",
                        "dtr", "cld") # "frs", "vap"

sites.sitenames <- c(BCI = "Barro_Colorado_Island_(BCI)", 
                        CedarBreaks = "",
                        HarvardForest = "Harvard_Forest",
                        LillyDickey = "Lilly_Dickey_Woods",
                        SCBI = "Smithsonian_Conservation_Biology_Institute_(SCBI)",
                        ScottyCreek = "Scotty_Creek")[sites]

## analysis parameters ####

reference_date <- c(30, 8) # refday in slidingwin
window_range <- c(15, 0) #range in slidingwin

first_years_oulier_limit = 15

clim_var_group <- list(c("pre", "wet"),
                       c("tmp", "tmn", "tmx", "pet"),
                       c("dtr", "cld", "pet")
                       ) # see issue 14, PET is in both the TMP and DTR groups. If it comes out as the best in both groups (should always be for the same time frame), then there are only 2 candidate variables for the GAM.

# load data ####
## climate data ####

for(clim_v in climate_variables) {
  assign(clim_v, read.csv(paste0(path_to_climate_data, clim_v, ".1901.2016-ForestGEO_sites-8-18.csv"))
  )
}

## core data ####
all_Biol <- list()
for (site in sites) {
  all_Biol[[site]] <- read.csv(paste0(path_to_core_measurement_files, "cores_", site, ".csv"))
}


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
  
  
  ### combine all variables in one
  if(clim_v == climate_variables [1]) all_Clim <- x_long[, c(1:3)]
  else all_Clim <- merge(all_Clim, x_long[, c(1:3)], by = c("sites.sitename", "Date"), all = T)
  
}

### format date to dd/mm/yyyy
all_Clim$Date <- gsub("X", "", all_Clim$Date)
all_Clim$Date <- format(as.Date(all_Clim$Date , format = "%Y.%m.%d"), "%d/%m/%Y") 

### save prepared clim data ####
write.csv(all_Clim, "processed_data/Climate_data_all_sites.csv", row.names = F)


## cores ####
for (site in sites) {
  Biol <- all_Biol[[site]]
  
  ### format date to dd/mm/yyyy
  Biol$Date <- paste0("15/06/", Biol$Year) # dd/mm/yyyy setting up as june 15, ARBITRATY
  
  ### cut back 10 years tree was cored dead STILL NEED TO CODE ####
  
  
  ### find out cores with outliers ####
  coreID_with_outliers <- unique(Biol[Biol$core_measurement>10, ]$coreID)
  
  
  coreID_with_outliers_within_first_X_years <- NULL
  
  par(mfrow = c(4,5), mar = c(2,1,3,1), oma = c(2,3,0,0))
  for( t in coreID_with_outliers) {
    x <- Biol[Biol$coreID %in% t,]
    
    # is the outlier within the first X years?
    within_first_ten_years <- any(which(x$core_measurement > 10) < first_years_oulier_limit)
    
    if(within_first_ten_years) coreID_with_outliers_within_first_X_years <- c(coreID_with_outliers_within_first_X_years, t)
    
    plot(core_measurement~Year, data = x, type = "l", main = paste(site, t, x$tree_status[1]), xlab = "", col = ifelse(max(x$core_measurement) > 15, "red", "black"))
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
  
  for( t in coreID_with_outliers_within_first_X_years) {
    x <- Biol[Biol$coreID %in% t,]
    
    rows_to_remove <- rownames(x)[x$Year <= min(x$Year) + first_years_oulier_limit]
    Biol <- droplevels(Biol[!rownames(Biol) %in% rows_to_remove, ])
    
  }

  ## remove years that are before climate record (+ first few first months to be able to look at window before measurement) ####
  Biol <- Biol[Biol$Year >= min(as.numeric(substr(all_Clim$Date, 7, 10)))+  window_range[1]/12, ]
  
  ## remove years that are after climate record ####
  Biol <- Biol[Biol$Year <= max(as.numeric(substr(all_Clim$Date, 7, 10))), ]
  
  
  
  ## remove any cores that ave less than 30 years ####
  Biol <- Biol[Biol$coreID %in% names(which(table(Biol$coreID)>= 30)), ]
  ## save back into Biol ####
  
  all_Biol[[site]] <- Biol
  
  
}

## Run the Analysis ####
for(site in sites) {
  
  cat(paste("running analysis for", site, "...\n" ))
  Biol <- all_Biol[[site]]
  Clim <- droplevels(all_Clim[all_Clim$sites.sitename %in% sites.sitenames[site], ])
  
  for(what in c("log_core_measurement")) { # , "log_agb_inc"
    
    ## calculate residuals of spine measurement ~ year for each individual####
    Biol$residuals <- NA
    
    for( t in unique(Biol$coreID)) {
      x <- Biol[Biol$coreID %in% t, ]
      
      x$Y <- x[, switch(what, log_core_measurement = "core_measurement", log_agb_inc = "agb_inc")]
      
      first_year_removed <- any(is.na(x$Y))
      x <- x[!is.na(x$Y),] #remove NA (only first year of measurement for agb_inc)
      
      test <- gam(Y~ s(Year), data = x)
      
      
      if(rbinom(1, 1, 0.1)==1) {
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
      if(what %in% "log_agb_inc" & first_year_removed) {
        Biol[Biol$coreID %in% t, ]$residuals <- c(NA, x$residuals) 
      } else {
        Biol[Biol$coreID %in% t, ] <- x
      }
      
      
    }
    
    
    ## run slidingwin on residuals to find best time window and lin or quad for each variable ####
    if(nlevels(Biol$species_code) > 1) {
      baseline = "lmer(residuals ~ 1 + (1 | species_code) + (1 | coreID), data = Biol[!is.na(Biol$residuals), ])"
    
      } else 
        {
      baseline = "lmer(residuals ~ 1 + (1 | coreID), data = Biol[!is.na(Biol$residuals), ])"
    }
   
    
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
                           func = "quad", # c("lin","quad")
                           refday = reference_date,
                           cinterval = "month",
                           cdate = Clim$Date, bdate = Biol[!is.na(Biol$residuals), ]$Date) 
    
    ### find best window for each variable group
    results$combos
    best_results_combos <- do.call(rbind, lapply(clim_var_group, function(X) {
      x <- results$combos[results$combos$climate %in% X,]
      data.frame(model_ID = as.numeric(rownames(x)), x, stringsAsFactors = F)[which.min(x$DeltaAICc),]
    }))
    best_results_combos <- best_results_combos[!duplicated(best_results_combos), ]# remove one pet in case it is best in both tmp and drt groups.
    # best_results_combos <- do.call(rbind, by(results$combos, results$combos$climate, function(x) data.frame(model_ID = as.numeric(rownames(x)), x, stringsAsFactors = F)[which.min(x$DeltaAICc),]))
    
    best_results_combos <- best_results_combos[order(best_results_combos$DeltaAICc),]
    
    ### plot the results and save the signal into Biol ####
    for(i in best_results_combos$model_ID) {
      print(paste("adding climate data to Biol for model", i))
      #clear plotting device 
      dev.off()
      # plot the results
      plotall(dataset = results[[i]]$Dataset, 
              bestmodel = results[[i]]$BestModel,
              bestmodeldata = results[[i]]$BestModelData,
              title=paste((data.frame(lapply(results$combos[i,], as.character), stringsAsFactors=FALSE)), collapse = "_"))
      
      # save the plot
      dev.print(png, paste0('results/ALL_species_mixed_model_on_residuals/ALL_species_mixed_model_on_', gsub("log_", "", what), "_", paste((data.frame(lapply(results$combos[i,], as.character), stringsAsFactors=FALSE)), collapse = "_"), "_", site, '.png'),
                width = 10,
                height =8,
                units = "in",
                res = 300)
      
      # save the climate signal in Biol
      if(any(grepl("I\\(climate\\^2\\)", names( results[[i]]$BestModelData)))) {
        columns_to_add <- results[[i]]$BestModelData[, c("climate", "I(climate^2)")]
        names(columns_to_add) <- paste0(results$combos[i,]$climate, c("", "^2"))
      } else {
        columns_to_add <- results[[i]]$BestModelData[c("climate")]
        names(columns_to_add) <- results$combos[i,]$climate
      }
      
      to_add_to_Biol <-  cbind(Biol[!is.na(Biol$residuals), ], columns_to_add)
      Biol[, names(columns_to_add)] <- to_add_to_Biol[match(rownames(Biol), rownames(to_add_to_Biol)), names(columns_to_add)]
    }
    
    ## Output Biol and best_results_combos to use in different analysis ####
    write.csv(Biol, file = paste0("processed_data/core_data_with_best_climate_signal", ifelse(what %in% "log_agb_inc", "_AGB", ""), "_", site, ".csv"), row.names = F)
    # saveRDS(best_results_combos, file =  paste0("processed_data/best_results_combos", ifelse(what %in% "log_agb_inc", "_AGB", ""), "_", site, ".rds"))
    
    ## look at collinearity between climate variables and remove any variable with vif > 10 ####
    X <- Biol[, as.character((best_results_combos$climate))]
    X <- X[!duplicated(X),]
    
   
    usdm::vif(X)
    (vif_res <-  usdm::vifstep(X, th = 10))
    variables_to_keep <- as.character(vif_res@results$Variables)

    
    ## now do a species by species gam using log of raw measuremets, spline on dbh and year ####
    
    # create the gam formula
    
    full_model_formula <- switch(what, "log_core_measurement" =  paste("log_core_measurement ~  s(Year, bs ='re', by = treeID) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")),
                                 log_agb_inc = paste("log_agb_inc ~ s(dbh, k = 3) + s(Year, bs ='re', by = treeID) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")))
    
    
    # full_model_formula <- switch(what, "log_core_measurement" =  paste("log_core_measurement ~ s(dbh, k = 3) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")),
    #                              log_agb_inc = paste("log_agb_inc ~ s(dbh, k = 3) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")))
    
    
    ## identify what variables we should keep for each species, looking at the sum of AIC weights ####
    start_time <- Sys.time()
    
    sum_of_weights_for_each_term_by_sp <- NULL
    for(sp in unique(Biol$species_code)) {
      print(sp)
      x <- Biol[Biol$species_code %in% sp,]
      # x$tag <- factor(x$tag)
      x$log_core_measurement <- log(x$core_measurement+0.1)
      # x$log_agb_inc <-  log(x$agb_inc + 0.1)
      # x <- x[, c("dbh", "Year", "tag", what, variables_to_keep)]
      
      x <- x[, c("Year", "treeID", what, variables_to_keep)]
      
      x <- x[!is.na(x[, what]), ]
      fm1 <- gam(eval(parse(text = full_model_formula)), data = x, na.action = "na.fail") # using mixed model makes all variables important which I find suspicous. I feel that s(Year, by = tag) is enough to add some randomness by tag... + plus I am not even sure that actuallydoes what we want.
      dd <- dredge(fm1)
      dd$cw <- cumsum(dd$weight)
      
      # sum_of_weights_for_each_term <- dd[, grepl(paste(c(variables_to_keep, "dbh", "Year"), collapse = "|"), names(dd))]
      sum_of_weights_for_each_term <- dd[, grepl(paste(c(variables_to_keep,  "Year"), collapse = "|"), names(dd))]
      sum_of_weights_for_each_term <- apply(sum_of_weights_for_each_term, 2, function(x) sum(dd$weight[!is.na(x)]))
      sum_of_weights_for_each_term
      sum_of_weights_for_each_term_by_sp <- rbind(sum_of_weights_for_each_term_by_sp,
                                                  c(sum_of_weights_for_each_term))
      
      
      
      # get the results of the model that includes the variables that have sum of weight > 0.9
      
      if(sum(sum_of_weights_for_each_term > 0.9) == 0 ) {
        best_model <- gam(update(as.formula(full_model_formula), .~1), data = x, na.action = "na.fail") #intercept only model if none of the variables have 0.9 sum of weitgh
      } else {
        best_model <- gam(eval(parse(text = paste( what,  "~", paste(names(sum_of_weights_for_each_term)[sum_of_weights_for_each_term > 0.9], collapse = " + ")))), data = x, na.action = "na.fail")
      }
      
      
      
      # save results for individual species
      assign(paste0(sp, "_dd"), dd)
      assign(paste0(sp, "_best_model"), best_model)
      
      # remove x
      rm(x)
    }
    end_time <- Sys.time()
    
    
    (ellapsed_time <- difftime(end_time, start_time))
    
    rownames(sum_of_weights_for_each_term_by_sp) <- unique(Biol$species_code)
    
    sum_of_weights_for_each_term_by_sp
    library(lattice)
    
    levelplot(t(sum_of_weights_for_each_term_by_sp), 
              scales=list(x=list(rot=45)), 
              xlab = "parameter", 
              ylab = "species",
              legend = list(top = list(fun = grid::textGrob("Sum of Weights", y=0, x=1.09))))
    
    # save the plot
    dev.print(png, paste0('results/Species_by_species_GAMS_on_raw_data/Sum_of_AICweights_', what, "_", site, '.png'),
              width = 10,
              height =8,
              units = "in",
              res = 300)
    
    
    
    
   
    
    ## Plot response curves for each variable, with one curve per species ####
    # first remove any object starting by p_
    rm(list = ls()[grepl("^p_", ls())])
    #Create a custom color scale
    
    for(v in c( variables_to_keep)) { # c("dbh", variables_to_keep)
      print(v)
      
      ## predictions
      pt <- NULL
      for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
        best_model <- get(paste0(sp, "_best_model"))
        
        
        varying_x <- data.frame(varying_x = seq(min(Biol[Biol$species_code %in% sp, v], na.rm = T), max(Biol[Biol$species_code %in% sp, v], na.rm = T), length.out = 100)) ; colnames(varying_x) <- v
        # constant_variables <- c("dbh", variables_to_keep)[!c("dbh", variables_to_keep) %in% v]
        constant_variables <- c(variables_to_keep)[!c(variables_to_keep) %in% v]
        
        newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(Biol$", constant_variables, ", na.rm = T)", collapse = ", "), ",  Year = median(Biol$Year, na.rm = T), treeID = factor(Biol[Biol$species_code %in% sp,]$treeID[1]))"))), varying_x)
        
        # if(v %in% names(best_model$var.summary)) {
          pt <- rbind(pt, data.frame(newd, variable = v, species = sp, varying_x = newd[, v], predict.gam(best_model, newd, type = "response", exclude =grep("Year", sapply(best_model$smooth, "[[", "label"), value = T), se.fit = T), draw = v %in% names(best_model$var.summary)))
          
        # }
        
      } # ignore errors
      
      # if(!is.null(pt)) {
        pt$species <- factor(pt$species, levels = rownames(sum_of_weights_for_each_term_by_sp))
        pt$expfit <- exp(pt$fit)
        pt$lwr <- exp(pt$fit - 1.96 * pt$se.fit)
        pt$upr <- exp(pt$fit + 1.96 * pt$se.fit)
        
        
        p <- ggplot(data = pt[pt$draw,], aes(x = varying_x, y = expfit))
        if(v != "dbh") p <- p + geom_rect(xmin = mean(Biol[, v], na.rm = T) - sd(Biol[, v], na.rm = T), ymin = min(pt$lwr), xmax = mean(Biol[, v], na.rm = T) + sd(Biol[, v], na.rm = T), ymax = max(pt$upr), fill = "grey" , alpha=0.01) + geom_vline(xintercept = mean(Biol[, v], na.rm = T), col = "grey")
        
        time_window <- reference_date[2] - as.numeric(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")])
        time_window_prev <- time_window < 0
        time_window <- ifelse(time_window_prev, rev(1:12)[abs(time_window)], time_window)
        
        time_window_text <- paste0("\nfrom ", paste(paste0(ifelse(time_window_prev, "prev. ", "curr. "), month.abb[time_window]), collapse = "\nto "))
        
        p <- p + geom_line(aes(group = species, col = species)) +
          # scale_x_continuous(trans= ifelse(v %in% "dbh", 'log','identity')) +
          labs(title = paste0(v, ifelse(v %in% best_results_combos$climate, time_window_text, "")),
               x = v,
               y = "") + #"core measurements") +
          geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = species), alpha=0.25) + 
          scale_colour_hue(drop = F) + scale_fill_hue(drop = F) + 
          theme_classic()
        
        assign(paste0("p_", v), p +
                 theme(legend.position="none"))  
      #}
      
    }
    
    g_legend<-function(){
      a.gplot <- ggplot(data = pt, aes(x = varying_x, y = expfit)) + geom_line(aes(group = species, col = species)) +  geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = species), alpha=0.25)
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}
    
    
    # existing_plots <- paste0("p_",  c("dbh", variables_to_keep))
    existing_plots <- paste0("p_",  c(variables_to_keep))
    
    grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  get(x)), ncol = 3)),
                 g_legend(),
                 nrow = 1,
                 widths = c(10, 1))
    
    grid::grid.text(switch (what, log_core_measurement = "core measurement (mm)",
                            log_agb_inc = "AGB increment (Mg C)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90)
    
    
    
    # save plot
    dev.print(png, paste0('results/Species_by_species_GAMS_on_raw_data/ALL_variables_', what, "_", site, '.png'),
              width = 8,
              height =8,
              units = "in",
              res = 300)
    
  }
  
  
  # save environment ####
  save.image(file = paste0("results/", site, "_all_env.RData"))
}



