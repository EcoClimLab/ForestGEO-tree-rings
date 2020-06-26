# ---compare regular dendro analysis to climwin's one ---####

# clear environment ####
rm(list = ls())

# load libraries ####
library(caTools) # for runmean
library(bootRes)
library(dplR) # for read.rwl
library(climwin)

# source function ####

source("https://raw.githubusercontent.com/SCBI-ForestGEO/climate_sensitivity_cores/master/scripts/0-My_dplR_functions.R")

# path to data *** TO BE EDITED **** ####
## because the dendro repo is private I can't look into the url to find what folders/species we want to pull in... so I have to use absolute path... not ideal I can't find a vetter solution for now.
path_to_COFECHA <- "c:/Users/HerrmannV/Dropbox (Smithsonian)/GitHub/EcoClimLab/ForestGEO_dendro/data_processed/COFECHA/"

# decide what sites and species to run ####
sites_species <-  c("SCBI_qual", "SCBI_litu", "Zofin_ABAL", "CedarBreaks_PSME", "ScottyCreek_PIMA")[c(1,2,4)]
# sites <- unique(strsplit(sites_species, "_")[[1]][1])

paths_to_rwl<- list.files(path_to_COFECHA, pattern = "rwl", full.names = T)
paths_to_tab <-  list.files(path_to_COFECHA, recursive = T,  pattern = "tabs.txt", full.names = T)
paths_to_out <-  list.files(path_to_COFECHA, recursive = T,  pattern = "out.txt", full.names = T)




# Define parameters and variables ####
## variables_to_show ?
variables_to_show <- c("pet", "tmx")

## saving or not saving outputs ? ####
save.plots <- TRUE
save.result.table <- TRUE

## Define full.time.frame.end.years **** TO BE EDITED ****####
full.time.frame.end.years <- list(SCBI = 2009,
                                  Zofin = NULL,
                                  CedarBreaks = 2009, #using SCBI's for now but need to be changed
                                  ScotyyCreek = NULL)

## Define how to run it regarding the starting year ####
type.of.start.date <- c("1901_2009") # , "1920_1949", "1950_1979", "1980_2009"


## Define sss threshold ####
sss.threshold = 0.75

## Define start and end month for anlysis ####
start <- -4#-4 # April of previous year
end <- 8 # August of current year
reference_date <- c(30, 8) # refday in slidingwin
window_range <- c(16,0) #c(16, 0) #range in slidingwin 16 is april (-4 in mdcc)

# Load climate data ####
## climate data
all_Clim <- read.csv("processed_data/Climate_data_all_sites.csv")

### add year column
all_Clim$year <- as.numeric(format(as.Date(all_Clim$Date, format = "%d/%m/%Y"), "%Y"))
### add month column
all_Clim$month <- as.numeric(format(as.Date(all_Clim$Date, format = "%d/%m/%Y"), "%m"))


# load individual chronologies ####
for(ssp in sites_species) {
  x <- read.csv(paste0("processed_data/core_data_with_best_climate_signal/log_core_measurement/", strsplit(ssp, "_")[[1]][1], ".csv"), stringsAsFactors = F)
  x <- x[x$species_code %in% strsplit(ssp, "_")[[1]][2], c("residuals", "coreID", "Year", "Date")]
  assign(paste0(ssp, "_ind_chron"), x)
}

# Run analyses ####

dir.create(paste0("results/formal_comparison/"), showWarnings = F)

## get to SSS for each species ####
all_sss <- NULL

sd_coreres <- NULL # will store SD of the detrended chronologie
mean_core_raw_per_species <- NULL # will store the mean radius increment per indviduals

for(ssp in sites_species) {
  f_raw <- grep(ssp, paths_to_rwl, value = T)
  f_tab <- grep(ssp, paths_to_tab, value = T)
  f_out <- grep(ssp, paths_to_out, value = T)
  
  # get the raw data
  core_raw <- read.rwl(f_raw)
  mean_core_raw_per_species <- c(mean_core_raw_per_species, mean(apply(core_raw, 2, mean, na.rm = T)))
  
  # get the detrended data
  core <- read.table(f_tab, sep = "\t", h = T)
  core <- data.frame(res = core$res,  samp.depth = core$num, row.names = core$year)
  
  # output the SD of the detrended chronologies (see issue # 63 on GitHub)
  sd_coreres <- rbind(sd_coreres, data.frame(Species = ssp, SD = round(sd(core$res), 2)))
  
  # get the Subsample Signal Strength (sss as function of the number of trees in sample, the last one appearing in the "xxx_drop.rxl_out.txt files)
  
  sss <- readLines(f_out)
  sss <- sss[grep("sss", sss)]
  
  sss <- sss[grep("  sss:   ", sss)[c(rep(FALSE, 3*length(seq(grep("  sss:   ", sss)))/4), rep(TRUE, 1*length(seq(grep("  sss:   ", sss)))/4))]] # keep only last rows that have sss: in them
  
  sss <- sub("  sss:   ", "", sss)
  sss <- as.numeric(unlist(strsplit(sss, " " ))) # keep only numbers and store them as a vector
  
  sss <- data.frame(Species = ssp, "Num_of_trees" = 1:length(sss), sss)
  
  Year_to_Num_of_trees <- apply(core_raw, 1, function(x) sum(!is.na(x)))
  Year_to_Num_of_trees <- data.frame(Species = ssp, Year = as.numeric(names(Year_to_Num_of_trees)), Num_of_trees= Year_to_Num_of_trees)
  
  match(Year_to_Num_of_trees$Num_of_trees, sss$Num_of_trees)
  
  Year_to_Num_of_trees$sss <- NA
  for(i in 1:nrow(Year_to_Num_of_trees)) {
    
    Year_to_Num_of_trees$sss[i] <- rev(sss[sss$Num_of_trees <= Year_to_Num_of_trees$Num_of_trees[i],]$sss)[1]
    
  }
  
  sss <- Year_to_Num_of_trees
  
  assign(ssp, core)
  assign(paste0(ssp, "_sss"), sss)
  
  all_sss <- rbind(all_sss, sss)
  
}

# save SSS for all species 
write.csv(all_sss, file = paste0("results/formal_comparison/SSS_as_a_function_of_the_number_of_trees_in_sample.csv"), row.names = F)

# save sd_coreres for all species
write.csv(sd_coreres, file = paste0("results/formal_comparison/SD_of_each_detrended_chornologies.csv"), row.names = F)


# save mean radius increment 
## see: https://github.com/SCBI-ForestGEO/climate_sensitivity_cores/issues/62
write.csv(data.frame(Species = sites_species, mean_rad_inc = mean_core_raw_per_species), file = paste0("results/formal_comparison/mean_radius_increment.csv"), row.names = F)


## Define start and end year for analysis for each species ####

start.years.sss <- NULL # species specific
for(ssp in sites_species) {
  sss <- get(paste0(ssp, "_sss"))
  start.years.sss <- c(start.years.sss, sss[sss$sss >= sss.threshold, ]$Year[1])
}

# full.time.frame.end.year = full.time.frame.end.years[[site]]  # for now and later for sd.clim_fill_time_frame (in 0-My_dplR_functions.R)

## Plot SSS for the the decided threshold ####

if(save.plots) tiff(paste0("results/formal_comparison/SSS_as_a_function_of_the_number_of_trees_in_sample.tiff"), res = 150, width = 169, height = 169, units = "mm", pointsize = 10)

op <- par(mfrow = c(2, 1), oma = c(5, 5, 2, 0), mar = c(0, 0, 0, 1))

cols <- data.frame(col = rainbow(length(sites_species)), row.names = sites_species, stringsAsFactors = F)

years <- NULL
for(sp in levels(all_sss$Species)){
  x = all_sss[all_sss$Species %in% sp,]
  year <- x$Year[x$sss > sss.threshold][1]
  years <- c(years, year)
}

plot.nb <- 1

for(sp in levels(all_sss$Species)){
  full.time.frame.end.year <- full.time.frame.end.years[[strsplit(sp, "_")[[1]][1]]] 
  
  x <- all_sss[all_sss$Species %in% sp,]
  x <- x[x$Year <= full.time.frame.end.year,] 
  # n.core <- x$Num_of_trees[x$sss > sss.threshold][1]
  
  if(plot.nb %in% 1) {
    plot(Num_of_trees ~ Year, data = x, type = "l", col = cols[sp,], xlim = c(min(all_sss$Year), full.time.frame.end.year), ylim = range(all_sss$Num_of_trees), lwd = 2, log = "y", las = 1, ylab = "", xaxt = "n")
    mtext(side= 2 , "log(No. cores)", line = 3)
    axis(1, labels = F, tcl = 0.5)
    axis(1, labels = F, tcl = -0.5)
    mtext("a)", side = 1, line = -1, adj = 0.01, font = 2)
    
  } else {
    lines(Num_of_trees ~ Year, data = x, col = cols[sp,], lwd = 2)
  }
  
  plot.nb <- plot.nb +1
}

abline(v =years,  col = cols$col, lty = 2)
legend("topleft", col = cols$col, lty = 1, bty = "n", legend = paste(levels(all_sss$Species), years, sep = " - "), lwd = 2, cex = 0.8)

plot.nb <- 1

for(sp in levels(all_sss$Species)){
  
  full.time.frame.end.year <- full.time.frame.end.years[[strsplit(sp, "_")[[1]][1]]] 
  
  x <- all_sss[all_sss$Species %in% sp,]
  x <- x[x$Year <= full.time.frame.end.year,] 
  
  year <- x$Year[x$sss > sss.threshold][1]
  
  if(plot.nb %in% 1) {
    plot(sss ~ Year, data = x, type = "l", col = cols[sp,], xlim = c(min(all_sss$Year), full.time.frame.end.year), lwd = 2, las = 1, xaxt = "n")
    abline(v = year, lty = 2, col = cols[sp,])
    abline(h = 0.75, lty = 3)
    axis(1, labels = T, tcl = 0.5)
    axis(1, labels = F, tcl = -0.5)
    mtext(side= 2 , "sss", line = 3)
    mtext("b)", side = 1, line = -1, adj = 0.01, font = 2)
    
  } else {
    lines(sss ~ Year, data = x, col = cols[sp,], lwd = 2)
    abline(v = x$Year[x$sss > sss.threshold][1], lty = 2, col = cols[sp,])
  }
  plot.nb <- plot.nb +1
}

title(paste("SSS threshold =", sss.threshold), outer = T)
par(op)

if(save.plots) dev.off()
par(op)


## Run the 3 types of analysis to compare ####
mean_and_std_of_clim <- NULL
for(type.start in type.of.start.date) {
  
  print(type.start)
  
  dir.create(paste0("results/formal_comparison/", type.start, "/figures"), recursive = T, showWarnings = F)
  dir.create(paste0("results/formal_comparison/", type.start, "/tables"), recursive = T, showWarnings = F)
  
  
  for(f in sites_species) {
    print(f)
    
    full.time.frame.end.year <- full.time.frame.end.years[[strsplit(f, "_")[[1]][1]]] 
    
    if(type.start %in% "1901_2009") {
      start.years <- start.years.sss
      end.year <- full.time.frame.end.year
    }
    
    if(type.start %in% "1920_1949") {
      start.years <- ifelse(start.years.sss > 1920, start.years.sss, 1920)
      end.year <- 1949
    }
    
    if(type.start %in% "1950_1979") {
      start.years <- ifelse(start.years.sss > 1950, start.years.sss, 1950)
      end.year <- 1979
    }
    
    if(type.start %in% "1980_2009") {
      start.years <- ifelse(start.years.sss > 1980, start.years.sss, 1980)
      end.year <- 2009
    }
    
    # species_to_drop <- "none"
    # 
    # if(type.start %in% c( "1920_1949", "1950_1979", "1980_2009")) {
    #   # Species should be dropped when they have SSS<0.75 for any part of the period (which means we'll miss a lot in 1910-1939).
    #   species_to_drop <- site_sps[which(start.years != min(start.years))]
    #   start.years <- start.years[!site_sps %in% species_to_drop]
    #   
    #   
    # }
    
    
    # load species chronology data ####
    core <- get(f)
    
    #load individual chronology data ####
    ind_chron <- get(paste0(f, "_ind_chron"))
    
    # load climate data for corresponding site  ####
    clim <- all_Clim[all_Clim$site %in% strsplit(f, "_")[[1]][1], -which(names(all_Clim) %in% c("sites.sitename", "site"))] # , "Date"
    
    ### crop last year to full.time.frame.end.year
    clim <- clim[clim$year <= full.time.frame.end.year, ]
    
    ## mean and std of climate variables ####
    ## see https://github.com/SCBI-ForestGEO/climate_sensitivity_cores/issues/41
    
    columns_to_remove <- which(names(clim) %in% c("Date", "year", "month"))
    mean_and_std_of_clim <- rbind(mean_and_std_of_clim,
                                  data.frame(climate.data = "CRU", start.year = max(min(start.years), min(clim$year)), end.year,
                                             variable = colnames(clim[clim$year >=  min(start.years) & clim$year <= end.year, -columns_to_remove]),
                                             do.call(rbind, lapply(seq_along(clim[clim$year >=  min(start.years) & clim$year <= end.year, -columns_to_remove]),
                                                                   function(i) {
                                                                     
                                                                     X <- clim[clim$year >=  min(start.years) & clim$year <= end.year, -columns_to_remove][,i]
                                                                     v <- names(clim[clim$year >=  min(start.years) & clim$year <= end.year, -columns_to_remove])[i]
                                                                     X.year <- clim[clim$year >=  min(start.years) & clim$year <= end.year, ]$year
                                                                     
                                                                     X.month <- clim[clim$year >=  min(start.years) & clim$year <= end.year, ]$month
                                                                     
                                                                     temp.month <- tapply(X, X.month, function(x) return(data.frame(mean = mean(x), sd = sd(x))))
                                                                     names(temp.month) <- month.abb[as.numeric(names(temp.month))]
                                                                     temp.month <- as.data.frame(do.call(cbind, temp.month))
                                                                     
                                                                     temp.annual <- tapply(X, X.year, function(x) {
                                                                       if(v %in% c("pre", "PCP", "pet_sum", "wet")) return(sum(x))
                                                                       if(v %in% "PETminusPRE")  return(sum(x[x>0]))
                                                                       if(!v %in% c("pre","PCP", "pet_sum", "wet", "PETminusPRE"))  return(mean(x))
                                                                     })
                                                                     
                                                                     
                                                                     temp.annual <- data.frame(Annual.mean = mean(temp.annual), Annual.sd = sd(temp.annual))
                                                                     
                                                                     
                                                                     return(cbind(temp.month, temp.annual))
                                                                     
                                                                   }))))
    
    ## trim measurement years ####
    ## remove years of core measurement that are before climate record (+ first few first months to be able to look at window before measurement)
    core <- core[as.numeric(rownames(core)) >= min(as.numeric(clim$year))+  window_range[1]/12, ]
    ind_chron <- ind_chron[ind_chron$Year>= min(as.numeric(clim$year))+  window_range[1]/12, ]
    
    ## remove years that are after climate record
    core <- core[as.numeric(rownames(core)) <= max(as.numeric(clim$year)), ]
    ind_chron <- ind_chron[ind_chron$Year <= max(as.numeric(clim$year)), ]
    
    start.year <- max(min(clim$year), min(as.numeric(rownames(core))))# max(min(clim$year), start.years[which(site_sps[!site_sps %in% species_to_drop] %in% f)])
    
    # run analysis for each variable
    for (v in variables_to_show) { # names(clim)[-columns_to_remove]
      print(v)
      
      ## traditional ####
      corr.dcc.output <- my.dcc(core, clim[, c("year", "month", v)], method = "correlation", start = ifelse(v %in% "frs", start.frs, start), end = ifelse(v %in% "frs", end.frs, end), timespan = c(start.year, end.year), ci = 0.05, ci2 = 0.002)
      
      # resp.dcc.output.v <- my.dcc(core, clim[, c("year", "month", v)], method = "response", start = ifelse(v %in% "frs", start.frs, start), end = ifelse(v %in% "frs", end.frs, end), timespan = c(start.year, end.year), ci = 0.05, ci2 = 0.002)
      
      ## climwin on species chronologies ####
      baseline_sp = "lm(res ~ 1, data = core)"
      # mean_scaled_x <- function(x) mean(scale(x)) # this to average the scaled variables fir each window tested, instead of using just the scaled data.
      climwin.output_sp <- slidingwin( baseline = eval(parse(text = baseline_sp)),
                                    xvar = as.list(clim[v]), 
                                    type = "absolute", 
                                    range = window_range,
                                    stat = c("mean"),
                                    func = "lin", # c("lin","quad")
                                    refday = reference_date,
                                    cinterval = "month",
                                    cdate = clim$Date, bdate =  paste0("15/06/", rownames(core))) 
      
      climwin.response_sp <- climwin.output_sp[[1]]$Dataset
      climwin.response_sp <- climwin.response_sp[climwin.response_sp$WindowOpen == climwin.response_sp$WindowClose, ]
      climwin.response_sp <- climwin.response_sp[order(climwin.response_sp$WindowOpen, decreasing = T), c(2,3,4,5)]
      
   
      
      ## plot climwin species chronology model ####
      if(save.plots)  png(paste0('results/formal_comparison/', type.start, '/figures/climwin_sp_', f, "_", paste((data.frame(lapply(climwin.output_sp$combos[1, c(2, 5, 7, 8)], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.png'),
                          width = 10,
                          height =8,
                          units = "in",
                          res = 300)
      plotall(dataset = climwin.output_sp[[1]]$Dataset, 
              bestmodel = climwin.output_sp[[1]]$BestModel,
              bestmodeldata = climwin.output_sp[[1]]$BestModelData,
              title=paste(data.frame(lapply(climwin.output_sp$combos[1,], as.character), stringsAsFactors=FALSE), collapse = "_"))
      
      
      if(save.plots) dev.off()
      
      ## climwin on individual chronologies ####
      baseline_ind = "lmer(residuals ~ 1 + (1|coreID), data = ind_chron)"

      climwin.output_ind <- slidingwin( baseline = eval(parse(text = baseline_ind)),
                                       xvar = as.list(clim[v]), 
                                       type = "absolute", 
                                       range = window_range,
                                       stat = c("mean"),
                                       func = "lin", # c("lin","quad")
                                       refday = reference_date,
                                       cinterval = "month",
                                       cdate = clim$Date, bdate =  ind_chron$Date) 
     
      
      climwin.response_ind <- climwin.output_ind[[1]]$Dataset
      climwin.response_ind <- climwin.response_ind[climwin.response_ind$WindowOpen == climwin.response_ind$WindowClose, ]
      climwin.response_ind <- climwin.response_ind[order(climwin.response_ind$WindowOpen, decreasing = T), c(2,3,4,5)]
      
      ## get beta coeficient of both climwin results but for the  best window with individual chronologies ####
      best_window_ind <- climwin.output_ind$combos[1, c("WindowOpen", "WindowClose")]
      beta_best_window <- c(ind = fixef(climwin.output_ind[[1]]$BestModel)[[2]],
                            sp = climwin.output_sp[[1]]$Dataset[climwin.output_sp[[1]]$Dataset$WindowOpen %in% best_window_ind$WindowOpen & climwin.output_sp[[1]]$Dataset$WindowClose %in% best_window_ind$WindowClose,]$ModelBeta)
      
      ## plot climwin individual chronology model ####
      
     
      if(save.plots)  png(paste0('results/formal_comparison/', type.start, '/figures/climwin_ind_', f, "_", paste((data.frame(lapply(climwin.output_ind$combos[1, c(2, 5, 7, 8)], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.png'),
                          width = 10,
                          height =8,
                          units = "in",
                          res = 300)
      plotall(dataset = climwin.output_ind[[1]]$Dataset, 
              bestmodel = climwin.output_ind[[1]]$BestModel,
              bestmodeldata = climwin.output_ind[[1]]$BestModelData,
              title=paste(data.frame(lapply(climwin.output_ind$combos[1,], as.character), stringsAsFactors=FALSE), collapse = "_"))
      
      
      if(save.plots) dev.off()
      
      
      
      
      ## plot barplots ####
      
      ### prepare ploting device ####
      if(save.plots)  png(paste0('results/formal_comparison/', type.start, '/figures/climwin_vs_dcc_', f, "_",v, '.png'),
                          width = 10,
                          height =8,
                          units = "in",
                          res = 300)
      
      par(mfrow = c(1, 3), mar = c(4,4,1,0), oma = c(0,0,1,1))
      
      colors_plot <- c("black", "grey40", "grey70")
     ### prepare data to plot ####
      
       barplot_matrix <- rbind(corr.dcc.output$chg_rad_inc_clim, 
                              climwin.response_sp$ModelBeta,
                              climwin.response_ind$ModelBeta
      )
      barplot_ci <- list(
        lower = rbind(corr.dcc.output$chg_rad_inc_clim_ci.lower,
                      climwin.response_sp$ModelBeta - (1.96*climwin.response_sp$Std.Error),
                      climwin.response_ind$ModelBeta - (1.96*climwin.response_ind$Std.Error)),
        upper = rbind(corr.dcc.output$chg_rad_inc_clim_ci.upper,
                      climwin.response_sp$ModelBeta + (1.96*climwin.response_sp$Std.Error),
                      climwin.response_ind$ModelBeta + (1.96*climwin.response_ind$Std.Error)))
      barplot_border <- colors_plot[c(1,2,3)]
      barplot_col <- rbind(ifelse(corr.dcc.output$significant, colors_plot[1], "transparent"), 
                           ifelse((climwin.response_sp$ModelBeta - (1.96*climwin.response_sp$Std.Error)) <0 & (climwin.response_sp$ModelBeta + (1.96*climwin.response_sp$Std.Error)) >0, "transparent",  colors_plot[2]),
                           ifelse((climwin.response_ind$ModelBeta - (1.96*climwin.response_ind$Std.Error)) <0 & (climwin.response_ind$ModelBeta + (1.96*climwin.response_ind$Std.Error)) >0, "transparent", colors_plot[3])
      )
     
      ### plot ###
      b <- barplot(barplot_matrix, col = barplot_col, border = barplot_border , beside = T, ylim = c(min(barplot_ci$lower), max(barplot_ci$upper)))
      
      segments(b, barplot_ci$lower, b, barplot_ci$upper, col = "grey")
      
      
      ylim_plot <- range(c(climwin.response_sp$ModelBeta, climwin.response_ind$ModelBeta))
      plot(corr.dcc.output$chg_rad_inc_clim, climwin.response_sp$ModelBeta, xlab = "dcc", ylab = "climwin", ylim = ylim_plot)
      abline(0, 1)
      points(corr.dcc.output$chg_rad_inc_clim, climwin.response_ind$ModelBeta, xlab = "dcc", pch = 16)
      
      legend("topleft", bty = "n",
             pch = c(1,16),
             legend = c("climwin chron",
                        "climwin ind"))
      
      # barplot best month climwin ind
      
      barplot(beta_best_window, col = colors_plot[c(2, 3)])
      
      if(save.plots) dev.off()
      
      
     
    
    }# for (v in names(clim)[-c(1:2)])
    
    
  } # for(f in sites_species)
  
  ### clean and save results###
  
  corr.dcc.output$variable <- sapply(strsplit(row.names(corr.dcc.output), "\\."), function(x) x[1])
  corr.dcc.output$month <- sapply(strsplit(row.names(corr.dcc.output), "\\."), function(x) paste(x[2], x[3], sep ="."))
  # all.dcc.output$month <- gsub("[0-9]", "",   all.dcc.output$month)
  
  # resp.dcc.output.v$variable <- sapply(strsplit(row.names(resp.dcc.output.v), "\\."), function(x) x[1])
  # resp.dcc.output.v$month <- sapply(strsplit(row.names(resp.dcc.output.v), "\\."), function(x) paste(x[2], x[3], sep ="."))
  
  if(save.result.table) {
    write.csv(corr.dcc.output, file = paste0("results/formal_comparison/", type.start, "/tables/correlation_with_" , "CRU_climate_data.csv"), row.names = F)
    # write.csv(resp.dcc.output.v, file = paste0("results/formal_comparison/", type.start, "/tables/response_to_" , "CRU_climate_data.csv"), row.names = F)
  } 
  
} # for(type.start in type.of.start.date) 


# save mean_and_std_of_clim ####
# dir.create(paste0("results/formal_comparison/", site, "/climate"), showWarnings = F, recursive = T)
# write.csv(mean_and_std_of_clim, file = paste0("results/formal_comparison/", site, "/climate/mean_and_std_of_climate_variables.csv"), row.names = F)







