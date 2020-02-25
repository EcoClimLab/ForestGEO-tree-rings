
# clear environment ####
rm(list = ls())

# load libraries ####
library(climwin)
library(lme4)

# set parameters ####
core_type <- c("CSV") # options are "detrended_chronologies" or "CSV"
reference_date <- c(30, 8) # refday in slidingwin
window_range <- c(15, 0) #range in slidingwin

if(core_type == "CSV") {
  baseline = "lmer(core_measurement ~ 1 + dbh + (1 | sp) + (1 | coreID), data = Biol)"
  baseline_by_sp <- "lmer(core_measurement ~ 1 + dbh + (1 | coreID), data = Biol[Biol$sp %in% sp,])"
}

if(core_type == "detrended_chronologies") {
  baseline = "lmer(res ~ 1 + (1 | sp), data = Biol)"
  baseline_by_sp <- "lm(res ~ 1 , data = Biol[Biol$sp %in% sp,])"
  
}

species <- c("CACO", "CAGL", "CAOVL", "CATO", "FAGR", "FRAM", "FRNI", "JUNI", "LITU", "PIST", "QUAL", "QUPR", "QURU", "QUVE")


# load data ####

## core data ####

if(core_type == "CSV") {
  cores <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_cores/cross-dated_cores_CSVformat/all_core_chronologies.csv")
  
  dbh_2008 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem1.csv", stringsAsFactors = FALSE)
  
  dbh_2013 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem2.csv", stringsAsFactors = FALSE)

  dbh_2018 <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_main_census/data/census-csv-files/scbi.stem3.csv", stringsAsFactors = FALSE)
}

if(core_type == "detrended_chronologies") {
  cores <- NULL
  for(f in species) {
    # get the raw data
    # core_raw <- dplR:: read.rwl(paste0("https://raw.githubusercontent.com/SCBI-ForestGEO/climate_sensitivity_cores/master/data/cores/", f,"/", tolower(f), "_drop.rwl"))
    # mean_core_raw_per_species <- c(mean_core_raw_per_species, mean(apply(core_raw, 2, mean, na.rm = T)))

    # get the detrended data
    core <- read.table(paste0("https://raw.githubusercontent.com/SCBI-ForestGEO/climate_sensitivity_cores/master/data/cores/", f,"/ARSTANfiles/", tolower(f), "_drop.rwl_tabs.txt"), sep = "\t", h = T)
    core <- data.frame(res = core$res,  samp.depth = core$num, row.names = core$year)
    
    ## can do more things, look at orginal script in https://github.com/SCBI-ForestGEO/climate_sensitivity_cores/blob/master/scripts/1-Calculate_and_plot_correlations_and_responses_between_tree-ring_chronologies_and_climate_variables.R##
    
    cores <- rbind(cores, data.frame(sp = f, core))
  }
}

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

if(core_type == "CSV") {
  
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

  
}


if(core_type == "detrended_chronologies") {
  Biol <- cores
  Biol$Date <- paste0("15/06/", rownames(Biol))
}

### remove years that are before climate record (+ first few first months to be able to look at window before measurement)
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) >= min(as.numeric(substr(Clim$Date, 7, 10)))+  window_range[1]/12, ]

### remove years that are after climate record
Biol <- Biol[as.numeric(as.numeric(substr(Biol$Date, 7, 10))) <= max(as.numeric(substr(Clim$Date, 7, 10))), ]


## plot core measurements splines against dbh and year ####

if(core_type == "CSV") {
  # plot against dbh, one page per tag
  for(t in unique(Biol$tag)) {
    x <- Biol[Biol$tag %in% t,]
    plot(core_measurement~ dbh,
         ylim = c(0,10),
         data= x, main = paste(x$tag[1], x$status.at.coring[1], sep = " - "))
    fit <- smooth.spline(x$dbh, x$core_measurement,df=4)
    lines(fit, col = "red")
  }
  
  # plot against year, one page per tag
  for(t in unique(Biol$tag)) {
    x <- Biol[Biol$tag %in% t,]
    plot(core_measurement~ Year,
         ylim = c(0,10),
         data= x, main = paste(x$tag[1], x$status.at.coring[1], sep = " - "))
    fit <- smooth.spline(x$Year, x$core_measurement,df=4)
    lines(fit, col = "red")
  }
  
  
  # plot against dbh, one page per sp
  for (sp in unique(Biol$sp)) {
    X <-   Biol[Biol$sp %in% sp,]
    X <- X[X$core_measurement < 15, ]
    X <- X[order(X$status.at.coring), ]
    for(t in unique(X$tag)) {
      x <- X[X$tag %in% t,]
      pch = 16 #ifelse(x$status.at.coring[1] == "alive", 16, 1)
      col = ifelse(x$status.at.coring[1] == "alive", 1, 2)
      
      if(t %in%  unique(X$tag)[1])  {
        plot(core_measurement~ dbh,
             xlim = c(min(X$dbh), max(X$dbh))
             ,ylim = c(min(X$core_measurement), max(X$core_measurement)),
             pch = pch,
             col = rgb(t(col2rgb(col)), alpha = 50, maxColorValue = 255),
             data= x, main = sp[1]) 
      } else {
        points(core_measurement~ dbh, data = x, pch = pch, col =  rgb(t(col2rgb(col)), alpha = 50, maxColorValue = 255))
      } 
      fit <- smooth.spline(x$dbh, x$core_measurement,df=4)
      lines(fit, col = col)
    }
    
    abline(lm(core_measurement ~ dbh, data = X), lwd = 2)
    lines(smooth.spline(X$dbh, X$core_measurement,df = 10), lwd = 2, col = "darkgreen")
    
    dev.print(tiff, paste0('results/explorations/', sp, "_core_meas_DBH_spline", '.tif'),
              width = 10,
              height =8,
              units = "in",
              res = 300)
  }
  
  # plot against year, one page per sp
  for (sp in unique(Biol$sp)) {
    X <-   Biol[Biol$sp %in% sp,]
    X <- X[X$core_measurement < 15, ]
    X <- X[order(X$status.at.coring), ]
    for(t in unique(X$tag)) {
      x <- X[X$tag %in% t,]
      pch = 16 #ifelse(x$status.at.coring[1] == "alive", 16, 1)
      col = ifelse(x$status.at.coring[1] == "alive", 1, 2)
      
      if(t %in%  unique(X$tag)[1])  {
        plot(core_measurement~ Year,
             xlim = c(min(X$Year), max(X$Year))
             ,ylim = c(min(X$core_measurement), max(X$core_measurement)),
             pch = pch,
             col = rgb(t(col2rgb(col)), alpha = 50, maxColorValue = 255),
             data= x, main = sp[1]) 
      } else {
        points(core_measurement~ Year, data = x, pch = pch, col =  rgb(t(col2rgb(col)), alpha = 50, maxColorValue = 255))
      } 
      fit <- smooth.spline(x$Year, x$core_measurement,df=4)
      lines(fit, col = col)
    }
    
    abline(lm(core_measurement ~ Year, data = X), lwd = 2)
    lines(smooth.spline(X$Year, X$core_measurement,df=10), lwd = 3, col = "darkgreen")
    
    
    dev.print(tiff, paste0('results/explorations/', sp, "_core_meas_YEAR_spline", '.tif'),
              width = 10,
              height =8,
              units = "in",
              res = 300)
    
  }
  
}




# run slidingwin on all species with sp as random effect ####

results <- slidingwin( baseline = eval(parse(text = baseline)),
                       xvar =list(cld = Clim$cld, 
                                  # dtr = Clim$dtr, 
                                  # frs = Clim$frs,# removed wet because model failed to converge 
                                  pet = Clim$pet, 
                                  pre = Clim$pre, 
                                  tmn = Clim$tmn, 
                                  tmp = Clim$tmp, 
                                  tmx = Clim$tmx, 
                                  vap = Clim$vap, 
                                  wet = Clim$wet
                                  ),
                       type = "absolute", 
                       range = window_range,
                       stat = c("mean"),
                       func = c("lin","quad"),
                       refday = reference_date,
                       cinterval = "month",
                       cdate = Clim$Date, bdate = Biol$Date) 

results$combos
best_mod_first_step <- which.min(results$combos$DeltaAICc)
results$combos[best_mod_first_step,]
results[[best_mod_first_step]][[1]]

randomized1 <- randwin(repeats = 100,     
                       baseline = eval(parse(text = baseline)),
                       xvar = list(Clim[,as.character(results$combos$climate[best_mod_first_step])]),
                       type = results$combos$type[best_mod_first_step], 
                       range = window_range,
                       stat = results$combos$stat[best_mod_first_step],
                       func = results$combos$func[best_mod_first_step],
                       refday = reference_date,
                       cinterval = "month",
                       cdate = Clim$Date, bdate = Biol$Date,
                       window= "sliding")


output <- results
pvalue(datasetrand = randomized1[[1]], dataset = output[[best_mod_first_step]]$Dataset, metric = "AIC", sample.size = nrow(Biol))


plotall(datasetrand = randomized1[[1]],
        dataset = output[[best_mod_first_step]]$Dataset, 
        bestmodel = output[[best_mod_first_step]]$BestModel,
        bestmodeldata = output[[best_mod_first_step]]$BestModelData,
        title=paste((data.frame(lapply(output$combos[best_mod_first_step,], as.character), stringsAsFactors=FALSE)), collapse = "_"))

## save outputs ####
dev.print(tiff, paste0('results/ALL_species_mixed_model_', paste((data.frame(lapply(output$combos[best_mod_first_step,], as.character), stringsAsFactors=FALSE)), collapse = "_"), '.tif'),
          width = 10,
          height =8,
          units = "in",
          res = 300)


save(results, best_mod_first_step, file = paste0('results/ALL_species_mixed_model_', ifelse(core_type == "CSV", "core_meas", "res"), ".RData"))

write.csv(results$combos[order(results$combos$DeltaAICc),], 
          file = paste0('results/ALL_species_mixed_model_', ifelse(core_type == "CSV", "core_meas_", "res_"), "combos.csv"))


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


# try gam model ####
plot(core_measurement~ dbh,
     data= Biol[Biol$sp %in% "qual",])

lines(smooth.spline(Biol[Biol$sp %in% "qual",]$dbh,Biol[Biol$sp %in% "qual",]$core_measurement,df=4), col = "red")

set.seed(0)
dat <- gamSim(1,n=200)
b<-gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat)
plot(b,pages=1)
gam.check(b,pch=19,cex=.3)

set.seed(0)
n<-200;sig2<-4
x1 <- runif(n, 0, 1);x2 <- runif(n, 0, 1);x3 <- runif(n, 0, 1)
fac<-c(rep(1,n/2),rep(2,n/2)) # create factor
fac.1<-rep(0,n)+(fac==1);fac.2<-1-fac.1 # and dummy variables
fac<-as.factor(fac)
f1 <-  exp(2 * x1) - 3.75887
f2 <-  0.2 * x1^11 * (10 * (1 - x1))^6 + 10 * (10 * x1)^3 * (1 - x1)^10
f<-f1*fac.1+f2*fac.2+x2
e <- rnorm(n, 0, sqrt(abs(sig2)))
y <- f + e
# NOTE: smooths will be centered, so need to include fac in model....
b<-gam(y~fac+s(x1,by=fac)+x2) 
plot(b,pages=1)

null.space.dimension(2,0)

n <- 100; set.seed(2)
x <- runif(n); y <- x + x^2*.2 + rnorm(n) *.1


## is smooth significantly different from straight line?
g1 <- gam(y~s(x,m=c(2,0))+x,method="REML") 
summary(g1)## not quite
plot(g1, pages = 1)

## is smooth significatly different from zero?
g2 <- gam(y~s(x),method="REML")
summary(g2) ## yes!
plot(g2, pages = 1)

## 
g3 <- gam(y~s(x, m=c(2,0)),method="REML")
summary(g3)
plot(g3)

## Fool bam(...,discrete=TRUE) into (strange) nested
## model fit...
set.seed(2) ## simulate some data... 
dat <- gamSim(1,n=400,dist="normal",scale=2)
dat$x1a <- dat$x1 ## copy x1 so bam allows 2 copies of x1
## Following removes identifiability problem, by removing
## linear terms from second smooth, and then re-inserting
## the one that was not a duplicate (x2)...
b <- bam(y~s(x0,x1)+s(x1a,x2,m=c(2,0))+x2,data=dat,discrete=TRUE)
summary(b)
plot(b, pages = 1)

test <- bam(core_measurement~ 1 + s(dbh, k = 4) + s(Year, by = tag),
              data= Biol[Biol$sp %in% "qual",])
plot(test, pages = 10)

for(sp in unique(Biol$sp)) {
  test <- bam(core_measurement~ 1 + dbh + s(Year, by = tag),
              data= Biol[Biol$sp %in% sp,])
  summary(test)
  plot(test, pages = 10, main = sp)
}



test <- gam(core_measurement~ 1 + s(dbh, k = -1),
            data= Biol[Biol$sp %in% "qual",])
plot(test, pages = 1)


test <- gam(core_measurement~ 1 + s(dbh, k = 4) + s(Year, by = tag),
            data= Biol[Biol$sp %in% "qual",])
plot(test, pages = 100)
gam.check(test)



test <- gamm4(core_measurement~ 1 + dbh + s(Year, k=4) , random=list(coreID=~1),
             data= Biol[Biol$sp %in% "qual",])
test <- gamm(core_measurement~ 1 + s(dbh, k=4) + s(Year, by=coreID), random=list(coreID=~1),
             data= Biol[Biol$sp %in% "qual",])


# other stuff ####
# sp = "QUAL"
# Biol$QUAL_signal <- NA
# Biol[Biol$sp %in% sp, ]$QUAL_signal <- output[[best_mod_first_step]]$BestModelData$climate
# 
# print(sp)
# results <- slidingwin( baseline = lm(res~1 + QUAL_signal + I(QUAL_signal^2), data = Biol[Biol$sp %in% sp,]),
#                        xvar =list(
#                                   tmx = Clim$tmx
#                        ),
#                        type = "absolute", 
#                        range = window_range,
#                        stat = c("mean"),
#                        func = c("lin","quad"),
#                        refday = reference_date,
#                        cinterval = "month",
#                        cdate = Clim$Date, bdate = Biol$Date[Biol$sp %in% sp]) 
# 
# results$combos
# best_mod_first_step <- which.min(results$combos$DeltaAICc)
# results$combos[best_mod_first_step,]
# results[[best_mod_first_step]][[1]]
# 
# best_mod_first_step_by_sp <- rbind(best_mod_first_step_by_sp, data.frame(sp = sp, results$combos[best_mod_first_step,]))
# 
# randomized1 <- randwin(repeats = 100,     
#                        baseline = lm(res~1 + QUAL_signal + I(QUAL_signal^2), data = Biol[Biol$sp %in% sp,]),
#                        xvar = list(Clim[,as.character(results$combos$climate[best_mod_first_step])]),
#                        type = results$combos$type[best_mod_first_step], 
#                        range = window_range,
#                        stat = results$combos$stat[best_mod_first_step],
#                        func = results$combos$func[best_mod_first_step],
#                        refday = reference_date,
#                        cinterval = "month",
#                        cdate = Clim$Date, bdate = Biol$Date[Biol$sp %in% sp],
#                        window= "sliding")
# 
# 
# output <- results
# pvalue(datasetrand = randomized1[[1]], dataset = output[[best_mod_first_step]]$Dataset, metric = "AIC", sample.size = nrow(Biol))
# 
# 
# plotall(datasetrand = randomized1[[1]],
#         dataset = output[[best_mod_first_step]]$Dataset, 
#         bestmodel = output[[best_mod_first_step]]$BestModel,
#         bestmodeldata = output[[best_mod_first_step]]$BestModelData,
#         title = paste(sp, as.character(results$combos$climate[best_mod_first_step]), sep = " - ")
# )
# 
# 
# Biol$drt <- output[[best_mod_first_step]]$BestModelData$climate
# plot(core_measurement ~ drt, data = Biol[Biol$coreID %in% Biol$coreID[1], ])
# 
# for(coreID in levels(Biol$coreID)) {
#   plot(core_measurement ~ drt, data = Biol[Biol$coreID %in% coreID, ])
#   points(predict(lm(core_measurement ~ poly(drt, 2), data = Biol[Biol$coreID %in% coreID,])) ~ Biol[Biol$coreID %in% coreID,]$drt, pch = 16)
#   
# }



