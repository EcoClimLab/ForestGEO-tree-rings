# ---run regular dendro analysis for all sites ---####

# clear environment ####
rm(list = ls())

# load libraries ####
library(caTools) # for runmean
library(bootRes)
library(dplR) # for read.rwl

# source function ####

source("https://raw.githubusercontent.com/SCBI-ForestGEO/climate_sensitivity_cores/master/scripts/0-My_dplR_functions.R")

# path to data *** TO BE EDITED **** ####
## because the dendro repo is private I can't look into the url to find what folders/species we want to pull in... so I have to use absolute path... not ideal I can't find a vetter solution for now.
path_to_COFECHA <- "c:/Users/HerrmannV/Dropbox (Smithsonian)/GitHub/EcoClimLab/ForestGEO_dendro/data_processed/COFECHA/"

# find out sites and species we can run ####
sites_species <- list.dirs(path_to_COFECHA, recursive = F, full.names = F)
sites <- unique(sapply(strsplit(sites_species, "_"), "[", 1))

paths_to_rwl<- list.files(path_to_COFECHA, pattern = "rwl", full.names = T) 
paths_to_tab <-  list.files(path_to_COFECHA, recursive = T,  pattern = "tabs.txt", full.names = T)
paths_to_out <-  list.files(path_to_COFECHA, recursive = T,  pattern = "out.txt", full.names = T)



# Define parameters and variables ####

## saving or not saving outputs ? ####
save.plots <- TRUE
save.result.table <- TRUE

## Define order of the species in the  plots, based on ANPP contribution####
# ANPP_contribution <- read.csv(text=getURL("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/summary_data/ANPP_total_and_by_species.csv"), header=T) 

# SPECIES_IN_ORDER <- toupper(ANPP_contribution$species[ ANPP_contribution$species %in% c("litu", "qual", "quru", "quve", "qupr", "fram", "cagl", "caco", "cato", "juni", "fagr", "caovl", "pist", "frni")])
# SPECIES_IN_ORDER <- gsub("CAOVL", "CAOV", SPECIES_IN_ORDER)

## Define sets of methods to run ####

methods.to.run <- c("correlation") # c("correlation", "response", "moving_correlation")

## Define full.time.frame.end.years ####
full.time.frame.end.years <- list(SCBI = 2009, CedarBreaks = 2009)
## Define how to run it regarding the starting year ####
type.of.start.date <- c("1901_2009", "1920_1949", "1950_1979", "1980_2009")


## Define sss threshold ####
sss.threshold = 0.75

## Define start and end month for anlaysis ####
start <- -4 # April of previous year
end <- 8 # August of current year

# start.frs <- -10 # october of previous year (for freeze days variable only - otherwise error because all 0 in other months)
# end.frs <- 5 # may of current year (for freeze days variable only)


# Load climate data ####
## climate data
all_Clim <- read.csv("processed_data/Climate_data_all_sites.csv")
### add year column
all_Clim$year <- as.numeric(format(as.Date(all_Clim$Date, format = "%d/%m/%Y"), "%Y"))
### add year column
all_Clim$month <- as.numeric(format(as.Date(all_Clim$Date, format = "%d/%m/%Y"), "%m"))

# Load and prepare core data ####
for(site in sites){
  
  dir.create(paste0("results/traditional_analysis/", site), showWarnings = F)
  
  # filenames <- list.dirs("data/cores/", full.names = F, recursive = F  )
  # filenames <- filenames[!grepl("[a-z]", filenames)] # keep only all caps names


all_sss <- NULL

sd_coreres <- NULL # will store SD of the detrended chronologie
mean_core_raw_per_species <- NULL # will store the mean radisu increment per indviduals

site_sps <- grep(site, sites_species, value = T)
for(ssp in site_sps) {
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

write.csv(all_sss, file = paste0("results/traditional_analysis/", site, "/SSS_as_a_function_of_the_number_of_trees_in_sample.csv"), row.names = F)

# save sd_coreres for all species
write.csv(sd_coreres, file = paste0("results/traditional_analysis/", site, "/SD_of_each_detrended_chornologies.csv"), row.names = F)


# save mean radius increment 
## see: https://github.com/SCBI-ForestGEO/climate_sensitivity_cores/issues/62
write.csv(data.frame(Species = site_sps, mean_rad_inc = mean_core_raw_per_species), file = paste0("results/traditional_analysis/", site, "/mean_radius_increment.csv"), row.names = F)


## Define start and end year for analysis, common to all species and one for each species ####

start.years.sss <- NULL # species specific
for(ssp in site_sps) {
  sss <- get(paste0(ssp, "_sss"))
  start.years.sss <- c(start.years.sss, sss[sss$sss >= sss.threshold, ]$Year[1])
}

full.time.frame.end.year = full.time.frame.end.years[[site]]  # for now and later for sd.clim_fill_time_frame (in 0-My_dplR_functions.R)

# Plot SSS for the the decided threshold ####

if(save.plots) png(paste0("results/traditional_analysis/", site, "/SSS_as_a_function_of_the_number_of_trees_in_sample.png"), res = 150, width = 169, height = 169, units = "mm", pointsize = 10)

op <- par(mfrow = c(2, 1), oma = c(5, 5, 2, 0), mar = c(0, 0, 0, 1))

cols <- data.frame(col = rainbow(length(site_sps)), row.names = site_sps, stringsAsFactors = F)

years <- NULL
for(sp in levels(all_sss$Species)){
  x = all_sss[all_sss$Species %in% sp,]
  year <- x$Year[x$sss > sss.threshold][1]
  years <- c(years, year)
}

plot.nb <- 1

for(sp in levels(all_sss$Species)){
  
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

# Run analysis ####
mean_and_std_of_clim <- NULL

  ## Load climate data for site  ####
  
  clim <- all_Clim[all_Clim$site %in% site, -which(names(all_Clim) %in% c("sites.sitename", "site", "Date"))]

  ### crop last year to full.time.frame.end.year
  clim <- clim[clim$year <= full.time.frame.end.year, ]
  
  ### get a moving average and sd of climate varibales, by month (for moving correlation)
  
  clim.moving.avg <- NULL
  clim.moving.sd <- NULL
  
  for(mth.i in (unique(clim$month))) {
    mth <- tolower(month.abb[mth.i])
    x.clim <- clim[clim$month %in% mth.i, ]
    x.clim.ma <- apply(x.clim, 2, runmean, k = 25, endrule = "NA", align = "center")
    x.clim.msd <- apply(x.clim, 2, runsd, k = 25, endrule = "NA", align = "center")
    
    rownames(x.clim.ma) <- paste(x.clim$year-12, x.clim$year + 12, sep = "-")
    rownames(x.clim.msd) <- paste(x.clim$year-12, x.clim$year + 12, sep = "-")
    
    clim.moving.avg[[mth]] <- x.clim.ma
    clim.moving.sd[[mth]] <- x.clim.msd
  }
  
  
  ## Run analysis on core data ####
  
  for(type.start in type.of.start.date) {
    
    print(type.start)
    
    dir.create(paste0("results/traditional_analysis/", site, "/", type.start, "/figures"), recursive = T, showWarnings = F)
    dir.create(paste0("results/traditional_analysis/", site, "/", type.start, "/tables"), recursive = T, showWarnings = F)
    
    
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
    
    species_to_drop <- "none"
    
    if(type.start %in% c( "1920_1949", "1950_1979", "1980_2009")) {
      # Species should be dropped when they have SSS<0.75 for any part of the period (which means we'll miss a lot in 1910-1939).
      species_to_drop <- site_sps[which(start.years != min(start.years))]
      start.years <- start.years[!site_sps %in% species_to_drop]
      # site_sps <- site_sps[!site_sps %in% species_to_drop]
      # SPECIES_IN_ORDER <- SPECIES_IN_ORDER[!SPECIES_IN_ORDER %in% gsub("CAOVL", "CAOV", species_to_drop)]
    }
    
    
    ## mean and std of climate variables ####
    ## see https://github.com/SCBI-ForestGEO/climate_sensitivity_cores/issues/41
    
    columns_to_remove <- which(names(clim) %in% c("year", "month"))
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
    
    ## Run analysis on core data ####
    
    for(method.to.run in methods.to.run) {
      
      print(method.to.run)
      
      dir.create(paste0("results/traditional_analysis/", site, "/", type.start, "/tables/monthly_", method.to.run ), showWarnings = F)
      
      
      all.dcc.output <- NULL
      
      ### run analysis ###
      for(f in site_sps[!site_sps %in% species_to_drop]) {
        print(f)
        
        core <- get(f)
        
        core <- core[rownames(core) %in% clim$year, ] # trim to use only years for which with have clim data 
        
        start.year <- max(min(clim$year), start.years[which(site_sps[!site_sps %in% species_to_drop] %in% f)])
        
        dcc.output <- NULL
        
        for (v in names(clim)[-columns_to_remove]) {
          print(v)
          
          if(method.to.run %in% c("correlation", "response")) {
            dcc.output <- rbind(dcc.output, my.dcc(core, clim[, c("year", "month", v)], method = method.to.run, start = ifelse(v %in% "frs", start.frs, start), end = ifelse(v %in% "frs", end.frs, end), timespan = c(start.year, end.year), ci = 0.05, ci2 = 0.002))
          }
          
          if(method.to.run %in% "moving_correlation" & type.start %in% "1901_2009") {
            all.dcc.output[[f]][[v]] <- bootRes::mdcc(core, clim[, c("year", "month", v)], method = "corr", start = ifelse(v %in% "frs", start.frs, start), end = ifelse(v %in% "frs", end.frs, end), timespan = c(start.year, end.year), win.size = 25, win.offset = 1, startlast = T,  boot = TRUE, ci = 0.05)
          }
          
        }# for (v in names(clim)[-c(1:2)])
        
        if(method.to.run %in% c("correlation", "response")) {
          all.dcc.output <- rbind(all.dcc.output, data.frame(cbind(Species = f, dcc.output)))
        }
        
      } # for(f in site_sps)
      
      ### clean and save results###
      if(method.to.run %in% c("correlation", "response")) {  
        all.dcc.output$variable <- sapply(strsplit(row.names(all.dcc.output), "\\."), function(x) x[1])
        all.dcc.output$month <- sapply(strsplit(row.names(all.dcc.output), "\\."), function(x) paste(x[2], x[3], sep ="."))
        all.dcc.output$month <- gsub("[0-9]", "",   all.dcc.output$month)
        
        if(save.result.table) write.csv(all.dcc.output, file = paste0("results/traditional_analysis/", site, "/", type.start, "/tables/monthly_", method.to.run, "/", method.to.run, ifelse(grepl("corr", method.to.run), "_with_", "_to_") , "CRU_climate_data.csv"), row.names = F)
      }
    
      
      ## Plot results ####
      
      if(method.to.run %in% c("correlation")) {
        for(v in names(clim)[-columns_to_remove]) {
          print(v)
          
          X <- all.dcc.output[all.dcc.output$variable %in% v, ]
          
          x <- data.frame(reshape(X[, c("month", "Species", "coef")], idvar = "month", timevar = "Species", direction = "wide"))
          rownames(x) <- ifelse(grepl("curr",  rownames(x)), toupper(rownames(x)), tolower( rownames(x)))
          rownames(x) <- gsub(".*curr.|.*prev.", "",   rownames(x), ignore.case = T)
          
          x.sig <- reshape(X[, c("month", "Species", "significant")], idvar = "month", timevar = "Species", direction = "wide")
          x.sig2 <- reshape(X[, c("month", "Species", "significant2")], idvar = "month", timevar = "Species", direction = "wide")
          
          colnames(x) <- gsub("coef.", "", colnames(x))
          colnames(x.sig) <- gsub("significant.", "", colnames(x.sig))
          colnames(x.sig2) <- gsub("significant2.", "", colnames(x.sig2))
          
          x <- x[, -1]
          x.sig <- x.sig[, -1]
          x.sig2 <- x.sig2[, -1]
          
          # x <- x[, rev(SPECIES_IN_ORDER[!SPECIES_IN_ORDER %in% gsub("CAOVL", "CAOV", species_to_drop)])]
          # x.sig <- x.sig[, rev(SPECIES_IN_ORDER[!SPECIES_IN_ORDER %in% gsub("CAOVL", "CAOV", species_to_drop)])]
          # x.sig2 <- x.sig2[, rev(SPECIES_IN_ORDER[!SPECIES_IN_ORDER %in% gsub("CAOVL", "CAOV", species_to_drop)])]
          
          if(save.plots)  {
            dir.create(paste0("results/traditional_analysis/", site, "/", type.start, "/figures/monthly_", method.to.run), showWarnings = F)
            dir.create(paste0("results/traditional_analysis/", site, "/", type.start, "/figures/monthly_", method.to.run), showWarnings = F)
            png(paste0("results/traditional_analysis/", site, "/", type.start, "/figures/monthly_", method.to.run, "/", v, ".png"), res = 150, width = 169, height = 169, units = "mm", pointsize = 10)
          }
          
          v <-  toupper(v)
          v <- gsub("PDSI_PREWHITEN" , "PDSI", v)
          
          my.dccplot(x = as.data.frame(t(x)), sig = as.data.frame(t(x.sig)), sig2 = as.data.frame(t(x.sig2)),  main = ifelse(v %in% "PETminusPRE", "PET-PRE", v), method = method.to.run)
          
          if(save.plots) dev.off()
        } #   for(v in names(clim)[-columns_to_remove])
      } # if(method.to.run %in% c("correlation") 
    
    } #  for(method.to.run in methods.to.run)
  } # for(type.start in type.of.start.date) 
  
  
  
} # for(site in sites)


# save mean_and_std_of_clim ####
dir.create(paste0("results/traditional_analysis/", site, "/climate"), showWarnings = F, recursive = T)
write.csv(mean_and_std_of_clim, file = paste0("results/traditional_analysis/", site, "/climate/mean_and_std_of_climate_variables.csv"), row.names = F)




