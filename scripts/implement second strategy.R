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

### plot the results and save the signal into Biol
for(i in best_results_combos$model_ID) {
  
  # plot the results
  plotall(dataset = results[[i]]$Dataset, 
          bestmodel = results[[i]]$BestModel,
          bestmodeldata = results[[i]]$BestModelData,
          title=paste((data.frame(lapply(results$combos[i,], as.character), stringsAsFactors=FALSE)), collapse = "_"))
  # save the plot
  dev.print(tiff, paste0('results/ALL_species_mixed_model_on_', paste((data.frame(lapply(results$combos[i,], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.tif'),
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

# look at collinearity betwwen climate variables
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

# now do a species by species gam with log normal familly, spline on dbh and year ####
library(MuMIn)

start_time <- Sys.time()
full_model_formula <- paste("log_core_measurement ~ s(dbh, k = 4) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + "))

# full_model_formula <- paste("log_core_measurement ~ ns(log(dbh), 2) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + "))

sum_of_weights_for_each_term_by_sp <- NULL
for(sp in unique(Biol$sp)) {
  print(sp)
  x <- Biol[Biol$sp %in% sp,]
  x$tag <- factor(x$tag)
  x$log_core_measurement <- log(x$core_measurement+0.1)
  x <- x[, c("dbh", "Year", "tag", "log_core_measurement", variables_to_keep)]
  fm1 <- gam(eval(parse(text = full_model_formula)), data = x, na.action = "na.fail") # using mixed model makes all variables important which I find suspicous. I feel that s(Year, by = tag) is enough to add some randomness by tag... + plus I am not even sure that actuallydoes what we want.
  dd <- dredge(fm1)
  dd$cw <- cumsum(dd$weight)
  
  sum_of_weights_for_each_term <- dd[, grepl(paste(c(variables_to_keep, "dbh", "Year"), collapse = "|"), names(dd))]
  sum_of_weights_for_each_term <- apply(sum_of_weights_for_each_term, 2, function(x) sum(dd$weight[!is.na(x)]))
  sum_of_weights_for_each_term
  sum_of_weights_for_each_term_by_sp <- rbind(sum_of_weights_for_each_term_by_sp,
                                              c(sum_of_weights_for_each_term))
  
  
 
  # get the results of the model that includes the variables that have sum of weight > 0.9
  best_model <- gam(eval(parse(text = paste( "log_core_measurement ~", paste(names(sum_of_weights_for_each_term)[sum_of_weights_for_each_term > 0.9], collapse = " + ")))), data = x, na.action = "na.fail")
  
  
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


# douvle check the best models we have by species is correct

best_models_by_species <- apply(sum_of_weights_for_each_term_by_sp, 1, function(x)paste(names(x)[x > 0.9], collapse = " + "))

formula(caco_best_model) == as.formula(paste("log_core_measurement ~", best_models_by_species[["caco"]]))


for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
  print(sp)
  X <- Biol[Biol$sp %in% sp, ]
  best_model <- get(paste0(sp, "_best_model"))

  variables_to_look_at <- names((best_model$var.summary))[!names((best_model$var.summary)) %in% c("Year", "tag")]
  
  n_row <- ifelse(length(variables_to_look_at) <= 2, 1, 2)
  n_col <- length(variables_to_look_at) %/% 2 +  length(variables_to_look_at)%% 2
  
  
  par(mfrow=c(ifelse(n_row == 0, 1, n_row), n_col))
  
  for(v in variables_to_look_at) {
    
    X$varying_v <- X[, v]
    varying_x <- data.frame(floor(min(X$varying_v)): ceiling(max(X$varying_v))) ; colnames(varying_x) <- v
    
    plot(core_measurement+0.1 ~ varying_v, data = X, log = "y", 
         pch = 16,
         # bg = rgb(0,0,0,0.2),
         col = rainbow(length(unique(X$tag)), s = 0.8, alpha = 0.2)[c(1:length(unique(X$tag)))[match(X$tag, unique(X$tag))]],
         main = paste0(sp[1], " - ", v, ifelse(v %in% best_results_combos$climate, paste0("\nfrom ",
                       paste(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")], collapse = " to ")), "")),
         xlab = v)
    
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
  
  
  
  # save plot
  dev.print(tiff, paste0('results/explorations/GAM_results_raw_', sp, ".tif"),
            width = 8,
            height =8,
            units = "in",
            res = 300)
  
}
















# x <- Biol[Biol$sp %in% sp, ]
# 
# plot(core_measurement ~ Year, data = x)
# 
# length(unique(Biol$tag[Biol$sp %in% "caco"]))
# gam.check(caco_best_model)
# summary(caco_best_model)
# caco_best_model$pterms
# caco_best_model$pred.formula
# family(caco_best_model)
# plot(caco_best_model)
# 
# vis.gam(caco_best_model, view = c("dbh", "pre"),  ticktype = "detailed")
# vis.gam(caco_best_model, view = c("dbh", "tmx"),  ticktype = "detailed")
# vis.gam(caco_best_model, view = c("Year", "dbh"),  ticktype = "detailed")
# vis.gam(caco_best_model, view = c("tmx", "pre")) #, cond=list(tag = "20605", dbh = 300))
# vis.gam(caco_best_model, view = c("tmx", "Year"))
# vis.gam(caco_best_model, view = c("tmx", "tag"))
# 
# caco_best_model$fitted.values <- exp(caco_best_model$fitted.values)
# vis.gam(caco_best_model, view = c("Year", "pre"), cond=list(tag = "20605", dbh = 300), ticktype = "detailed")
# 
# 
# caco_best_model$smooth[[1]]
# 
# caovl_best_model
# # visualize best models
# 
# t = 162252
# x <- Biol[Biol$tag %in% 20605, ]
# plot(log(x$core_measurement) ~ x$Year)
#
#
# all_models <- get.models(dd, subset = TRUE)
# summary(all_models[[1]])
# plot(all_models[[1]], pages = 1)
# plot(all_models[[1]])
# gam.check(all_models[[1]])
#
#
# 
# models_to_try <- c(null = "gam(log(core_measurement+0.1)~ 1 + s(dbh, k = 4) + s(Year, by = tag),
#             data= x)",
#                    pre = "gam(log(core_measurement+0.1)~ ns(pre, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    wet = "gam(log(core_measurement+0.1)~ ns(wet, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    cld = "gam(log(core_measurement+0.1)~ ns(cld, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    dtr = "gam(log(core_measurement+0.1)~ ns(dtr, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    tmp = "gam(log(core_measurement+0.1)~ ns(tmp, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    tmn = "gam(log(core_measurement+0.1)~ ns(tmn, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    pre_wet = "gam(log(core_measurement+0.1)~ ns(pre, 2) + ns(wet, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    pre_cld = "gam(log(core_measurement+0.1)~ ns(pre, 2) + ns(cld, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    pre_dtr = "gam(log(core_measurement+0.1)~ ns(pre, 2) + ns(dtr, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    pre_tmp = "gam(log(core_measurement+0.1)~ ns(pre, 2) + ns(tmp, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
#                    pre_tmn = "gam(log(core_measurement+0.1)~ ns(pre, 2) + ns(tmn, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)",
# )
# 
# 
# test_null <- gam(log(core_measurement+0.1)~ 1 + s(dbh, k = 4) + s(Year, by = tag),
#             data= x)
# 
# test_pre <- gam(log(core_measurement+0.1)~ ns(pre, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                data= x)
# 
# test_wet <- gam(log(core_measurement+0.1)~ ns(wet, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                 data= x)
# 
# test_pre_wet <- gam(log(core_measurement+0.1)~ ns(pre, 2) + ns(wet, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                 data= x)
# 
# test_pre_int <- gam(log(core_measurement+0.1)~ s(dbh, k = 4, by = (pre + I(pre^2))) + s(Year, by = tag),
#                 data= x)
# 
# test_wet_pre_int <- gam(log(core_measurement+0.1)~ poly(wet, 2) + s(dbh, k = 4, by = (pre + I(pre^2))) + s(Year, by = tag),
#                     data= x)
# 
# test_wet_int_pre_int <- gam(log(core_measurement+0.1)~s(dbh, k = 4, by = (pre + I(pre^2)) + wet + I(wet^2)) + s(Year, by = tag),
#                         data= x)
# 
# test_prewet <- gam(log(core_measurement+0.1)~ ns(pre, 2) * ns(wet, 2) + s(dbh, k = 4) + s(Year, by = tag),
#                     data= x)
# 
# plot(test_null, pages = 1)
# plot(test_pre, pages = 1)
# plot(test_wet, pages = 1)
# plot(test_pre_wet, pages = 1)
# plot(test_pre_int, pages = 1)
# plot(test_wet_pre_int)
# plot(test_wet_int_pre_int)
# plot(test_prewet)
# 
# AIC(test_null)
# AIC(test_pre)
# AIC(test_wet)
# AIC(test_pre_wet)
# AIC(test_pre_int)
# AIC(test_wet_pre_int)
# AIC(test_wet_int_pre_int)
# AIC(test_prewet)
# 
# summary(test_null)
# summary(test_pre)
# summary(test_wet)
# summary(test_pre_wet)
# summary(test_pre_int)
# summary(test_wet_pre_int)
# summary(test_wet_int_pre_int)
# summary(test_prewet)

