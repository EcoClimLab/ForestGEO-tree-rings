# climwin figures together to compare core increment, agb increment and BAI 

# clear environment ####
rm(list = ls())

# load libraries ####
library(png)
sites <- list.files("results/log_core_measurement", pattern = "env.RData")
sites <- gsub("_all_env.RData", "", sites)

# make plots ####
for(solution in c("/", "detrend_climate/", "old_records_only/", "young_records_only/")[]) {
  # clear folder of old plots ####
  dir.create(paste0("results/", solution, "climwin_plots_combined/"), recursive = T)
  file.remove(list.files(paste0("results/", solution, "climwin_plots_combined/"), full.names = T))
  
  # prepare site names ####
  
  
  # prepare order_plots ####
  order_plots <- list(log_core_measurement = expression(RW~"- all cores"), 
                      log_core_measurement_dbh= expression(RW~"- trees with DBH"), 
                      log_BAI_dbh = "BAI",
                      log_agb_inc_dbh = "AGB")
  
  # prepare variable groupe ####
  clim_var_group <- list(c("pre", "wet"),
                         c("tmp", "tmn", "tmx", "pet")#,
                         # c("dtr", "cld", "pet")
  )
  
  # prepare and save figures ####
  for (site in switch(solution, "/" = sites, c("ScottyCreek", "NewMexico", "SCBI"))) {
    imgs <- list.files(paste0("results/", solution, c("log_core_measurement", "log_core_measurement_dbh", "log_agb_inc_dbh", "log_BAI_dbh"), "/", site), pattern = "climwin", full.names = T)
    
    variables <- regmatches(imgs, regexpr("climwin_\\D{3}_", imgs))
    variables <- unique(gsub("climwin|_", "", variables))
    
    for(g in clim_var_group) {
      
      variables_in_group <- variables[variables %in% g]
      imgs_v <- imgs[unlist(sapply(variables_in_group, grep, imgs))]  
      nb_plots <- length(imgs_v)
      
      
      #order the plots
      imgs_v <- imgs_v[match(sapply(paste0(names(order_plots), "/"), grep, imgs_v), c(1,2,3,4))]
      imgs_v <- imgs_v[!is.na(imgs_v)]
      
      # prepare file
      png(paste0("results/", solution, "climwin_plots_combined/", site, "_", paste(variables_in_group, collapse = "_"), ".png" ), width = 10, height = 8, units = "in", res = 300)
      
      layout(matrix(1:(nb_plots*2), nrow = nb_plots, byrow = T), widths = c(4, 1))
      
      # plot
      
      for(i in seq_along(imgs_v)){
        y_lab = order_plots[[strsplit(imgs_v[i], "/")[[1]][3]]]
        v = gsub("^.*climwin_|_quad.*$", "", imgs_v[i])
        img <- readPNG(imgs_v[i])
        img1 <- img[1:(nrow(img)/2), ,] # keep only first half of plot
        img2 <- img[c((nrow(img)/2)+1):nrow(img), c((ncol(img1)/4): (2*ncol(img1)/4)),] # keep only plot 2 of second row
        
        
        par(mar = c(0,3,0,0))
        plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
        
        rasterImage( img1 , xleft = 0, xright = 100,
                     ybottom = 0, ytop = 100)
        
        mtext(side = 2, text = v)
        mtext(side = 2, text = y_lab, line = -2)
        
        par(mar = c(0,0,0,0))
        plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
        
        rasterImage( img2 , xleft = 0, xright = 100,
                     ybottom = 0, ytop = 100)
        
      }
      
      
      dev.off()
    }
  }
  
} # for(i in c("", "with_detrended_climate", "old_records_only"))

