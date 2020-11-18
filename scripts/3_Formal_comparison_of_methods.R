# ---compare regular dendro analysis to climwin's one ---####

# clear environment ####
rm(list = ls())

# load libraries ####
library(caTools) # for runmean
library(bootRes)
library(dplR) # for read.rwl
library(climwin)
library(lme4)


# source function ####

source("https://raw.githubusercontent.com/SCBI-ForestGEO/climate_sensitivity_cores/master/scripts/0-My_dplR_functions.R")

# path to data *** TO BE EDITED **** ####
## because the dendro repo is private I can't look into the url to find what folders/species we want to pull in... so I have to use absolute path... not ideal I can't find a vetter solution for now.
path_to_sp_res_chrons <- "c:/Users/HerrmannV/Dropbox (Smithsonian)/GitHub/EcoClimLab/ForestGEO_dendro/data_processed/sp_chronologies/0_CRNs/"
path_to_rwls <- "c:/Users/HerrmannV/Dropbox (Smithsonian)/GitHub/EcoClimLab/ForestGEO_dendro/data_processed/sp_chronologies/"

# decide what sites and species to run ####
sites_species <-  c("SCBI_litu", "ScottyCreek_PIMA", "Zofin_ABAL", "CedarBreaks_PSME")


# Define parameters and variables ####
## variables_to_show ?
variables_to_show <- c("pet", "tmx", "wet", "pre")

## saving or not saving outputs ? ####
save.plots <- TRUE
save.result.table <- TRUE

## Define full.time.frame.end.years ####
full.time.frame.end.years <- c(SCBI_litu = 2010,
                               Zofin_ABAL = 2010,
                               CedarBreaks_PSME = 2007, #using SCBI's for now but need to be changed
                               ScottyCreek_PIMA = 2013)

## Define start and end year for analysis for each species ####

start.years.sss <- c(SCBI_litu = 1919, # 11 trees
                     CedarBreaks_PSME = 1840, # 3 trees
                     Zofin_ABAL = 1676, # 7 trees (serious dating errors according to Neil)
                     ScottyCreek_PIMA = 1854 #13 trees
) # these dates were given by email by Neil Pederson when he created the psecies chrionologies for us... his threshold sss was .8


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


# load species chronologies ####
for(ssp in sites_species) {
  x <- read.csv(paste0(path_to_sp_res_chrons,ssp, "_col.csv"), stringsAsFactors = F, row.names = 1)
  assign(ssp, x)
}

# load individual chronologies (same as rwl files but in different format) ####
# par(mfrow = c(1, 4))
for(ssp in sites_species) {
  x <- read.csv(paste0("processed_data/core_data_with_best_climate_signal/log_core_measurement/", strsplit(ssp, "_")[[1]][1], ".csv"), stringsAsFactors = F)
  # hist(unique(x[ , c("coreID", "dbh_at_coring")])$dbh_at_coring, breaks = 20, main = "Histogram of dbh at coring", xlab = "")
  # summary(unique(x[ , c("coreID", "dbh_at_coring")])$dbh_at_coring)
  # 
  # hist(tapply(x$core_measurement, x$coreID, function(y) sum(y*2/10)), main = "Histogram of sum core measurements", xlab = "")
  # summary(tapply(x$core_measurement, x$coreID, function(y) sum(y*2/10)))
  # A <- tapply(x$dbh_at_coring, x$coreID, max)
  # if(ssp %in% "ScottyCreek_PIMA") x <- x[x$coreID %in% names(A)[order(A, decreasing = T)[1:50]], ]
  x <- x[grepl(strsplit(ssp, "_")[[1]][2], x$species_code, ignore.case = T), c("residuals", "coreID", "Year", "Date")]
  
  
  # plot(residuals ~ Year, data = x, ylim = c(-5, 5), main = paste(ssp, "\n", length(unique(x$coreID))))
  
  assign(paste0(ssp, "_ind_chron"), x)
}


# Run analyses ####

dir.create(paste0("results/formal_comparison/"), showWarnings = F)


# output the SD of the detrended chronologies and mean_core_raw_per_species
sd_cores <- NULL # will store SD of the detrended chronologie
# par(mfrow = c(1,4))
for(ssp in sites_species) {
  core <- get(ssp)
  
  sd_cores <- rbind(sd_cores, data.frame(Species = ssp, SD = round(sd(core$res), 2)))
  # plot(core$res ~ rownames(core), xlim = c(start.years.sss[ssp], full.time.frame.end.years[ssp]), main = paste(ssp, "\n", length(unique(get(paste0(ssp, "_ind_chron"))$coreID))))
}

# save sd_coreres for all species
write.csv(sd_cores, file = paste0("results/formal_comparison/SD_of_each_detrended_chornologies.csv"), row.names = F)






## Run the 3 types of analysis to compare ####

dir.create(paste0("results/formal_comparison/figures"), recursive = T, showWarnings = F)
file.remove(list.files(paste0("results/formal_comparison/figures"), full.names = T))

dir.create(paste0("results/formal_comparison/tables"), recursive = T, showWarnings = F)
file.remove(list.files(paste0("results/formal_comparison/tables"), full.names = T))


for(f in sites_species) {
  print(f)
  
  end.year <- full.time.frame.end.years[f] 
  start.year <- start.years.sss[f]
  
  
  # load species chronology data ####
  core <- get(f)
  
  #load individual chronology data ####
  ind_chron <- get(paste0(f, "_ind_chron"))
  
  # load climate data for corresponding site  ####
  clim <- all_Clim[all_Clim$site %in% strsplit(f, "_")[[1]][1], -which(names(all_Clim) %in% c("sites.sitename", "site"))] # , "Date"
  
  ### crop last year to full.time.frame.end.year
  clim <- clim[clim$year <= end.year, ]
  
  
  # trim measurement years ####
  ## remove years of core measurement that are before climate record (+ first few first months to be able to look at window before measurement)
  core <- core[as.numeric(rownames(core)) >= (min(as.numeric(clim$year)) +  window_range[1]/12), ]
  ind_chron <- ind_chron[ind_chron$Year>= (min(as.numeric(clim$year)) +  window_range[1]/12), ]
  
  ## remove years that are after climate record
  core <- core[as.numeric(rownames(core)) <= max(as.numeric(clim$year)), ]
  ind_chron <- ind_chron[ind_chron$Year <= max(as.numeric(clim$year)), ]
  
  start.year <- max(min(clim$year), start.year, min(ind_chron$Year))# max(min(clim$year), start.years[which(site_sps[!site_sps %in% species_to_drop] %in% f)])
  
  # run analysis for each variable
  for (v in  switch(f, "CedarBreaks_PSME" = c("tmx", "pre"), variables_to_show)) { # names(clim)[-columns_to_remove]
    print(v)
    
    ## traditional ####
    boolFalse<-F
    while(boolFalse==F) {
      tryCatch({
        corr.dcc.output <- tryCatch(my.dcc(chrono = core["res"], clim = clim[, c("year", "month", v)], method = "correlation", start = start, end =  end, timespan = c(start.year, end.year), ci = 0.05, ci2 = 0.002))
        boolFalse<-T
      },error=function(e){
      },finally={})
    } # doing a trycatch because "ScottyCreek_PIMA" failes often....
    
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
    if(save.plots)  png(paste0('results/formal_comparison/figures/climwin_sp_', f, "_", paste((data.frame(lapply(climwin.output_sp$combos[1, c(2, 5, 7, 8)], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.png'),
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
    beta_best_window <- rbind(sp = climwin.output_sp[[1]]$Dataset[climwin.output_sp[[1]]$Dataset$WindowOpen %in% best_window_ind$WindowOpen & climwin.output_sp[[1]]$Dataset$WindowClose %in% best_window_ind$WindowClose,],
                              ind = climwin.output_ind[[1]]$Dataset[1,])
    
    ## plot climwin individual chronology model ####
    
    
    if(save.plots)  png(paste0('results/formal_comparison/figures/climwin_ind_', f, "_", paste((data.frame(lapply(climwin.output_ind$combos[1, c(2, 5, 7, 8)], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.png'),
                        width = 10,
                        height =8,
                        units = "in",
                        res = 300)
    plotall(dataset = climwin.output_ind[[1]]$Dataset, 
            bestmodel = climwin.output_ind[[1]]$BestModel,
            bestmodeldata = climwin.output_ind[[1]]$BestModelData,
            title=paste(data.frame(lapply(climwin.output_ind$combos[1,], as.character), stringsAsFactors=FALSE), collapse = "_"))
    
    
    if(save.plots) dev.off()
    
    
    
    
    ## comparison plots ####
    
    ### prepare ploting device ####
    if(save.plots)  png(paste0('results/formal_comparison/figures/climwin_vs_dcc_', f, "_",v, '.png'),
                        width = 8,
                        height =4,
                        units = "in",
                        res = 300, 
                        pointsize = 14)
    
    par(mfrow = c(1, 3), mar = c(4,4,4,0), oma = c(0,0,1,1), mgp = c(2, 1, 0))
    
    colors_plot <- c("black", "grey40", "grey70")
    ### prepare data and plot ####
    
    barplot_matrix <- rbind(traditional_sp = corr.dcc.output$chg_rad_inc_clim, 
                            # climwin.response_sp$ModelBeta,
                            climwin_ind = climwin.response_ind$ModelBeta
    )
    
    barplot_ci <- list(
      lower = rbind(corr.dcc.output$chg_rad_inc_clim_ci.lower,
                    # climwin.response_sp$ModelBeta - (1.96*climwin.response_sp$Std.Error),
                    climwin.response_ind$ModelBeta - (1.96*climwin.response_ind$Std.Error)),
      upper = rbind(corr.dcc.output$chg_rad_inc_clim_ci.upper,
                    # climwin.response_sp$ModelBeta + (1.96*climwin.response_sp$Std.Error),
                    climwin.response_ind$ModelBeta + (1.96*climwin.response_ind$Std.Error)))
    
    barplot_border <- colors_plot[c(1,3)] #[c(1,2,3)]
    
    barplot_density <- rbind(ifelse(corr.dcc.output$significant, NA, 20), 
                             # ifelse((climwin.response_sp$ModelBeta - (1.96*climwin.response_sp$Std.Error)) <0 & (climwin.response_sp$ModelBeta + (1.96*climwin.response_sp$Std.Error)) >0, "transparent",  colors_plot[2]),
                             ifelse((climwin.response_ind$ModelBeta - (1.96*climwin.response_ind$Std.Error)) <0 & (climwin.response_ind$ModelBeta + (1.96*climwin.response_ind$Std.Error)) >0, 20, NA)
    )
    
    x.axis.labels <- reference_date[2] - c(window_range[1]:window_range[2])
    x.axis.labels_prev <- x.axis.labels <=0
    x.axis.labels <- ifelse(x.axis.labels_prev, rev(1:12)[abs(x.axis.labels) +1], x.axis.labels)
    
    x.axis.labels <- paste0(ifelse(x.axis.labels_prev, "p.", "c."),substr(month.abb, 1, 1)[x.axis.labels])
    
    ### plot ####
    b <- barplot(barplot_matrix, col = barplot_border, density = barplot_density, border = barplot_border, beside = T, ylim = c(min(barplot_ci$lower), max(barplot_ci$upper)), main = "Month-by-month", ylab = "Beta coeficient", legend.text = T, args.legend = list(fill = barplot_border, density = NULL, legend =  c("Traditional (sp)", "Climwin (ind)"), x = "top", bty = "n", border = "transparent", cex = .8))
    
    text(x = b[c(TRUE, FALSE)], y = par()$usr[3]*(1+0.19), labels =  x.axis.labels, srt = 45, cex = 0.8, xpd = NA)
    
    segments(b, barplot_ci$lower, b, barplot_ci$upper, col = barplot_border)
    
    # plot sp vs individual climwin only
    ylim_plot <- range(c(climwin.response_ind$ModelBeta - 1.96 * climwin.response_ind$Std.Error, climwin.response_ind$ModelBeta + 1.96 * climwin.response_ind$Std.Error))
    xlim_plot <- range(c(climwin.response_sp$ModelBeta - 1.96 * climwin.response_sp$Std.Error, climwin.response_sp$ModelBeta + 1.96 * climwin.response_sp$Std.Error))
    
    plot(x = climwin.response_sp$ModelBeta, y = climwin.response_ind$ModelBeta, xlab = "Species Chronology", ylab = "Individuals in Mixed-effects Model", ylim = ylim_plot, xlim = xlim_plot, main = "Comparison Beta coefficients\n1-month window Climwin", pch = 16, xpd = NA)
    segments(x0 = climwin.response_sp$ModelBeta - 1.96*climwin.response_sp$Std.Error, 
             x1 = climwin.response_sp$ModelBeta + 1.96*climwin.response_sp$Std.Error,
             y0 = climwin.response_ind$ModelBeta,
             y1 = climwin.response_ind$ModelBeta,
             col = "grey") # x axis CI
    segments(x0 = climwin.response_sp$ModelBeta,
             x1 = climwin.response_sp$ModelBeta,
             y0 = climwin.response_ind$ModelBeta - 1.96*climwin.response_ind$Std.Error,
             y1 = climwin.response_ind$ModelBeta + 1.96*climwin.response_ind$Std.Error,
             col = "grey") # y axis CI
    abline(0, 1, lty = 3)
    # text(x = climwin.response_sp$ModelBeta, y = climwin.response_ind$ModelBeta, labels = x.axis.labels, cex = 0.8, pos = 4, col = "grey")
    
    
    # points(x = corr.dcc.output$chg_rad_inc_clim, y = climwin.response_ind$ModelBeta, xlab = "dcc")
    # 
    # legend("topleft", bty = "n",
    #        pch = c(16,1),
    #        legend = c("climwin sp",
    #                   "climwin ind"))
    
    # plot best month climwin ind
    
    time_window <- reference_date[2] - as.numeric(best_window_ind[, c("WindowOpen", "WindowClose")])
    time_window_prev <- time_window <= 0
    time_window <- ifelse(time_window_prev, rev(1:12)[abs(time_window) +1], time_window)
    
    time_window_text <- paste(paste0(ifelse(time_window_prev, "p.", "c."),month.abb[time_window]), collapse = "-") #  paste0("\nfrom ", paste(paste0(ifelse(time_window_prev, "prev. ", "curr. "), month.abb[time_window]), collapse = "\nto "))
    
    ylim_plot <- range(c(beta_best_window$ModelBeta - 1.96 * beta_best_window$Std.Error, beta_best_window$ModelBeta + 1.96 * beta_best_window$Std.Error))
    
    
    plot(beta_best_window$ModelBeta,
         xlim = c(0,3), 
         ylim = ylim_plot,
         col = colors_plot[c(1, 3)], 
         main = paste("Climwin best window (ind)\n", time_window_text), 
         ylab = "Beta coeficient Climwin", 
         xlab = "", 
         xaxt = "n",
         pch = 16,
         bty = "n")
    segments(x0 = c(1,2), x1 = c(1,2),
             y0 = beta_best_window$ModelBeta - 1.96 * beta_best_window$Std.Error, 
             y1 = beta_best_window$ModelBeta + 1.96 * beta_best_window$Std.Error,
             col = colors_plot[c(1, 3)])
    
    
    axis(side = 1,
         at = c(1, 2), tick = F,
         labels = rownames(beta_best_window)
         
         
    )
    
    
    
    
    
    ### dev.off() ####
    
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
  write.csv(corr.dcc.output, file = paste0("results/formal_comparison/tables/correlation_with_" , "CRU_climate_data.csv"), row.names = F)
} 

# delete core_data_with_best_climate_signal ####
# file.remove("processed_data/core_data_with_best_climate_signal/")



# compile plots ####
library(png)
imgs <- list.files(paste("results/formal_comparison/figures/"), full.names = T)
for(ssp in sites_species) {
  
  for(v in  switch(f, "CedarBreaks_PSME" = c("tmx", "pre"), variables_to_show)) {
    
    
    # prepare file
    png(paste0("results/formal_comparison/", ssp, "_", v, ".png" ), width = 10, height = 8, units = "in", res = 300)
    
    layout(matrix(c(1,2,3,4,5,6), nrow = 2, byrow = F))
    
    # get the plots we need
    img_sp <- readPNG(imgs[grepl("sp", imgs) & grepl(paste(ssp), imgs) & grepl(paste(v), imgs)])
    img_ind <- readPNG(imgs[grepl("ind", imgs) & grepl(paste(ssp), imgs) & grepl(paste(v), imgs)])
    img_vs <- readPNG(imgs[grepl("vs", imgs) & grepl(paste(ssp), imgs) & grepl(paste(v), imgs)])
    
    # plot ####
    
    img1 <- img_sp[1:(nrow(img_sp)/2), (ncol(img_sp)/3):((ncol(img_sp)/3)*2),]
    img2 <- img_ind[1:(nrow(img_ind)/2), (ncol(img_ind)/3):((ncol(img_ind)/3)*2),] 
    img3 <- img_sp[(nrow(img_sp)/2):nrow(img_sp), (ncol(img_sp)/3):((ncol(img_sp)/3)*2),]
    img4 <- img_ind[(nrow(img_ind)/2):nrow(img_ind), (ncol(img_ind)/3):((ncol(img_ind)/3)*2),]
    img5 <- img_vs[, ((ncol(img_vs)/3)*2):ncol(img_vs),]
    img6 <- img_vs[, (ncol(img_vs)/3):((ncol(img_vs)/3)*2),]
    
    par(mar = c(0,0,0,0))
    for(i in 1:6) {
      plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
      
      rasterImage( get(paste0("img",i)) , xleft = 0, xright = 100,
                   ybottom = 0, ytop = 100)
    }
    
    
    
    # dev.off() ####
    
    dev.off()
  }
}


