# clear environment ####
rm(list = ls())

# load libraries ####
library(bootRes)
library(dplR) # for read.rwl
library(climwin)

# source function ####

source("https://raw.githubusercontent.com/SCBI-ForestGEO/climate_sensitivity_cores/master/scripts/0-My_dplR_functions.R")

path_to_sp_res_chrons <- "C:/Users/world/Documents/GitHub/growth_phenology/Data/tree_rings/Harvard/" # replace to the path where Neil's chronlogies are


#species <-  c("vector_species_you_have", "as they are names in Neils files")
species <-  c("ACRU", "BEAL", "QURU", "TSCA")

climate_variables <- c("tmx", "tmn")

## Define start and end year for analysis for each species ####

#start.years.sss <-c(Species1 = 2010,
#                    Species2 = 2010) # enter the name in species and the corresponding year the analysis should start at # these dates may have been given by email by Neil Pederson when he created the psecies chrionologies, it is a date at which sss passes a certain threshold sss (know what the threshold is .75, or .8?)

start.years.sss <-c(ACRU = 1930,#
                    BEAL = 1952,
                    QURU = 1898,
                    TSCA = 1930) # enter the name in species and the corresponding year the analysis should start at # these dates may have been given by email by Neil Pederson when he created the psecies chrionologies, it is a date at which sss passes a certain threshold sss (know what the threshold is .75, or .8?)

#full.time.frame.end.years <- c(Species1 = 2010,
#                               Species2 = 2010) # enter the name in species and the corresponding year the analysis should stop at
full.time.frame.end.years <- c(ACRU = 2013,#
                               BEAL = 2012,
                               QURU = 2013,
                               TSCA = 2013) # enter the name in species and the corresponding year the analysis should stop at

## Define start and end month for anlysis ####
start <- 1#January of current year
end <- 8 # August of current year

# Load climate data ####

## load it from here: https://github.com/forestgeo/Climate/tree/master/Climate_Data/CRU/CRU_v4_04
## needs to be in different format: column for Date, year and then one column per climate variable

# something like this should do
for(clim_v in climate_variables) {
  print(clim_v)
             x <- read.csv(paste0("https://raw.githubusercontent.com/forestgeo/Climate/master/Climate_Data/CRU/CRU_v4_04/", clim_v,  ".1901.2019-ForestGEO_sites-6-03.csv"))


  ### subset for the sites we care about
  x <- droplevels(x[x$sites.sitename %in% "Harvard_Forest", ])

  ### reshape to long format
  x_long <- reshape(x,
                    times = names(x)[-1], timevar = "Date",
                    varying = list(names(x)[-1]), direction = "long", v.names = clim_v)

  ### format date
  x_long$Date <- gsub("X", "", x_long$Date)
  x_long$Date <- as.Date(x_long$Date , format = "%Y.%m.%d")#changed format to work with Harvard data


  ### combine all variables in one
  if(clim_v == climate_variables[1]) all_Clim <- x_long[, c(1:3)]
  else all_Clim <- merge(all_Clim, x_long[, c(1:3)], by = c("sites.sitename", "Date"), all = T)

}

### add year column
all_Clim$year <- as.numeric(format(as.Date(all_Clim$Date, format = "%d/%m/%Y"), "%Y"))
### add month column
all_Clim$month <- as.numeric(format(as.Date(all_Clim$Date, format = "%d/%m/%Y"), "%m"))

# load species chronologies ####
for(ssp in species) {
  #the files I had were .txt files - converted to csv to work with read.csv (the only read function I know)
  x <- read.csv(paste0(path_to_sp_res_chrons,"TP_", ssp, ".rwl_tabs.csv"), stringsAsFactors = F, row.names = 1) # chanhe accorfingly to the names of the files you have (from Neil)

  assign(ssp, x)
}

# Run analyses ####

# output the SD of the detrended chronologies and mean_core_raw_per_species (if you want)
sd_cores <- NULL # will store SD of the detrended chronologie
# par(mfrow = c(1,4))
for(ssp in species) {
  core <- get(ssp)

  sd_cores <- rbind(sd_cores, data.frame(Species = ssp, SD = round(sd(core$res), 2)))
  # plot(core$res ~ rownames(core), xlim = c(start.years.sss[ssp], full.time.frame.end.years[ssp]), main = paste(ssp, "\n", length(unique(get(paste0(ssp, "_ind_chron"))$coreID))))
}

# save sd_coreres for all species
write.csv(sd_cores, file = paste0("results/SD_of_each_detrended_chornologies.csv"), row.names = F)






## Run analysis to compare ####
all.dcc.output <- NULL#
corr.dcc.output <- NULL#
for(f in species) {
  print(f)

  end.year <- full.time.frame.end.years[f]
  start.year <- start.years.sss[f]


  # load species chronology data ####
  core <- get(f)

  # load climate data for corresponding site (not necessary since you have only one site, but renaming to clim so that the rest works)  ####
  clim <- all_Clim

  ### crop last year to full.time.frame.end.year
  clim <- clim[clim$year <= end.year, ]


  # trim measurement years ####
  ## remove years of core measurement that are before climate record (+ first few first months to be able to look at window before measurement)
 #Wasn't sure what window_range was meant to be, so just removed it. I assume this will make the first year's correlation less reliable but since there is a lot more years it shouldn't have a big impact?
   core <- core[as.numeric(rownames(core)) >= (min(as.numeric(clim$year))),]#window_range[1]/12), ]


  ## remove years that are after climate record
  core <- core[as.numeric(rownames(core)) <= max(as.numeric(clim$year)), ]

  start.year <- max(min(clim$year), start.year)# max(min(clim$year), start.years[which(site_sps[!site_sps %in% species_to_drop] %in% f)])

  # run analysis for each variable
  for (v in  climate_variables) {
    print(v)


    corr.dcc.output <- my.dcc(chrono = core["res"], clim = clim[, c("year", "month", v)], method = "correlation", start = start, end =  end, timespan = c(start.year, end.year), ci = 0.05, ci2 = 0.002)
    all.dcc.output <- rbind(all.dcc.output, data.frame(cbind(Species = substr(f, 1, 4), corr.dcc.output)))#

  }

                                      ### plot ####

    # you should know ploting function

  }


all.dcc.output$variable <- substr(paste(row.names(all.dcc.output)), 1, 3)#get variable from row name
all.dcc.output$month <- substr(paste(row.names(all.dcc.output)), 5, 12)#get month from row name

write.csv(all.dcc.output, file = "results/Harvard_Forest_core_corr.csv", row.names = FALSE)
#} # for(f in species)

#############################################
##Copy/Paste this section from other script##
#############################################
save.plots = TRUE
for(v in climate_variables) {
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
#  x.sig <- x.sig[, rev(SPECIES_IN_ORDER[!SPECIES_IN_ORDER %in% gsub("CAOVL", "CAOV", species_to_drop)])]
#  x.sig2 <- x.sig2[, rev(SPECIES_IN_ORDER[!SPECIES_IN_ORDER %in% gsub("CAOVL", "CAOV", species_to_drop)])]

 # if(save.plots)  {
#    dir.create(paste0("results/", type.start, "/figures/monthly_", method.to.run), showWarnings = F)
#    dir.create(paste0("results/", type.start, "/figures/monthly_", method.to.run, "/", c), showWarnings = F)
#    tiff(paste0("results/", type.start, "/figures/monthly_", method.to.run, "/", c, "/", v, ".tif"), res = 150, width = 169, height = 169, units = "mm", pointsize = 10)
#  }

  v <-  toupper(v)
  v <- gsub("PDSI_PREWHITEN" , "PDSI", v)
  tiff(paste0("results/", "monthly_", "correlation", "harvard", v, ".tif"), res = 150, width = 169, height = 169, units = "mm", pointsize = 10)

  my.dccplot(x = as.data.frame(t(x)), sig = as.data.frame(t(x.sig)), sig2 = as.data.frame(t(x.sig2)),  main = ifelse(v %in% "PETminusPRE", "PET-PRE", v), method = "correlation")

  if(save.plots) dev.off()
}

