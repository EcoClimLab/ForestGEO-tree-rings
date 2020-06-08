# climwin figures together to compare core increment, agb increment and BAI 

# clear environment ####
rm(list = ls())

# load libraries ####
library(png)

# prepare site names ####
sites <- list.files("results", pattern = ".RData")
sites <- gsub("_all_env.RData", "", sites)

# preapre and save figures ####
for (site in sites) {
  imgs <- list.files(paste0("results/", c("log_core_measurement", "log_core_measurement_dbh", "log_agb_inc_dbh", "log_BAI_dbh"), "/", site), pattern = "climwin", full.names = T)
  
  variables <- regmatches(imgs, regexpr("climwin_\\D{3}_", imgs))
  variables <- unique(gsub("climwin|_", "", variables))
  
  for(v in variables) {
    
    imgs_v <- imgs[grepl(v, imgs)]                 
    nb_plots <- length(imgs_v)
   
    
    
    png(paste0("results/climwin_plots_combined/", site, "_",v, ".png" ), width = 10, height = 8, units = "in", res = 300)
    layout(matrix(1:(nb_plots*2), nrow = nb_plots, byrow = T), widths = c(4, 1))
   
    
    for(i in seq_along(imgs_v)){
      y_lab = gsub("log_|_inc", "", strsplit(imgs_v[i], "/")[[1]][2])
      img <- readPNG(imgs_v[i])
      img1 <- img[1:(nrow(img)/2), ,] # keep only first half of plot
      img2 <- img[c((nrow(img)/2)+1):nrow(img), c((ncol(img1)/4): (2*ncol(img1)/4)),] # keep only plot 2 of second row
      # img1[, (ncol(img1)+1):  (ncol(img1)+ncol(img1)/4), ] <- 0
      # 
      # img1[,c((ncol(img1)/4): (ncol(img1) - ncol(img1)/4)) , ] <- img1[,c( (2 *ncol(img1)/4): ncol(img1)) , ] #replace plot 2 and 3 by 3 and 4
      # img1[,c(  (ncol(img1) - ncol(img1)/4): ncol(img1)) , ] <- img[c((nrow(img)/2)+1):nrow(img), c((ncol(img1)/4): (2*ncol(img1)/4)),] # replace plot 4 by plot 6 in orginal image
      
      par(mar = c(0,3,0,0))
      plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
      
      rasterImage( img1 , xleft = 0, xright = 100,
                  ybottom = 0, ytop = 100)
      
      mtext(side = 2, text = y_lab)
      
      par(mar = c(0,0,0,0))
      plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
      
      rasterImage( img2 , xleft = 0, xright = 100,
                   ybottom = 0, ytop = 100)
      
    }
    
    dev.off()
    

  
  
  }
  
}