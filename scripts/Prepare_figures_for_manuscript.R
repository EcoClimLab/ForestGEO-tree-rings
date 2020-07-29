# --- create figures for manuscript --- #

# clear environment ####
rm(list = ls())

# load libraries ####
library(png)
library(gridExtra)
library(ggplot2)

# prepare site list and site coordinates ####
sites <- list.dirs("results/log_core_measurement_dbh", full.names = F, recursive = F)


sites_coords <- read.csv("https://raw.githubusercontent.com/forestgeo/Site-Data/master/ForestGEO_site_data.csv")
setdiff(sites, sites_coords$Site.name)
sites.sitenames <- c(BCI = "Barro Colorado Island", 
                     CedarBreaks = "Utah Forest Dynamics Plot",
                     HarvardForest = "Harvard Forest",
                     LillyDickey = "Lilly Dickey Woods",
                     SCBI = "Smithsonian Conservation Biology Institute",
                     ScottyCreek = "Scotty Creek",
                     Zofin = "Zofin",
                     HKK = "Huai Kha Khaeng",
                     NewMexico = "New_Mexico")[sites]
sites_coords <- sites_coords[sites_coords$Site.name %in% sites.sitenames, c("Site.name", "Latitude", "Longitude")]

# add coordinates of New Mexico site
sites_coords <- rbind(sites_coords, data.frame(Site.name = "New_Mexico", Latitude = 35.738376, Longitude = -105.838154 ))

# add site column 
sites_coords$site <- names(sites.sitenames)[match(sites_coords$Site.name, sites.sitenames)]

# order sites by latitude
sites <- sites_coords$site[order(sites_coords$Latitude)]
n_sites <- length(sites)

# DBH response at each sites and for each response ####
what_to_show <- c("log_core_measurement_dbh" = expression(Delta*r), "log_BAI_dbh" = expression(BAI), "log_agb_inc_dbh" = expression(Delta*AGB))


g_legend <- function(){
  a.gplot <- ggplot(data = p$data, aes(x = varying_x, y = expfit)) + geom_line(aes(group = species, col = species)) +  geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = species), alpha=0.25) + theme(legend.title=element_blank())+ guides(col=guide_legend(ncol=2), bg=guide_legend(ncol=2))
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  }

all_plots <- list()
for(site in sites){
  for(what in names(what_to_show)) {
    temp_env <- new.env()
    load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
    
    p <- get("p_dbh", envir = temp_env)
    
    # remove title and xlab of p
    p$labels$title <- NULL
    p$labels$x <- NULL
  
    
    # save into all_plots
    all_plots[[paste0(site, what)]] <- p
    } # for what in ...
  
  # add a plot for the legend
  all_plots[[paste0(site, "leg")]] <- g_legend()
} # for site in sites

png("doc/manuscript/tables_figures/DBH_responses.png", width = 8, height = 10, res = 300, units = "in")

grid.arrange(arrangeGrob(grobs = all_plots, ncol = 4, vp= grid::viewport(width=0.95, height=0.95)))


grid::grid.text(sites, x = unit(0.01, "npc"), y = unit(rev(cumsum(c(1/n_sites/2, rep(1/n_sites, n_sites-1)))), "npc"), rot = 90)

grid::grid.text(what_to_show, x = unit(cumsum(c(1/4/2, rep(1/4, 2))), "npc"), y = unit(.99,  "npc"))

grid::grid.text("dbh (cm)", x = unit(.5, "npc"), y = unit(0.015,  "npc"))

dev.off()

# # preapre and save figures ####
# for (site in sites) {
#   imgs <- list.files(paste0("results/", c("log_core_measurement", "log_core_measurement_dbh", "log_agb_inc_dbh", "log_BAI_dbh"), "/", site), pattern = "climwin", full.names = T)
#   
#   variables <- regmatches(imgs, regexpr("climwin_\\D{3}_", imgs))
#   variables <- unique(gsub("climwin|_", "", variables))
#   
#   for(v in variables) {
#     
#     imgs_v <- imgs[grepl(v, imgs)]                 
#     nb_plots <- length(imgs_v)
#     
#     
#     
#     png(paste0("results/climwin_plots_combined/", site, "_",v, ".png" ), width = 10, height = 8, units = "in", res = 300)
#     layout(matrix(1:(nb_plots*2), nrow = nb_plots, byrow = T), widths = c(4, 1))
#     
#     
#     for(i in seq_along(imgs_v)){
#       y_lab = gsub("log_|_inc", "", strsplit(imgs_v[i], "/")[[1]][2])
#       img <- readPNG(imgs_v[i])
#       img1 <- img[1:(nrow(img)/2), ,] # keep only first half of plot
#       img2 <- img[c((nrow(img)/2)+1):nrow(img), c((ncol(img1)/4): (2*ncol(img1)/4)),] # keep only plot 2 of second row
#       # img1[, (ncol(img1)+1):  (ncol(img1)+ncol(img1)/4), ] <- 0
#       # 
#       # img1[,c((ncol(img1)/4): (ncol(img1) - ncol(img1)/4)) , ] <- img1[,c( (2 *ncol(img1)/4): ncol(img1)) , ] #replace plot 2 and 3 by 3 and 4
#       # img1[,c(  (ncol(img1) - ncol(img1)/4): ncol(img1)) , ] <- img[c((nrow(img)/2)+1):nrow(img), c((ncol(img1)/4): (2*ncol(img1)/4)),] # replace plot 4 by plot 6 in orginal image
#       
#       par(mar = c(0,3,0,0))
#       plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
#       
#       rasterImage( img1 , xleft = 0, xright = 100,
#                    ybottom = 0, ytop = 100)
#       
#       mtext(side = 2, text = y_lab)
#       
#       par(mar = c(0,0,0,0))
#       plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
#       
#       rasterImage( img2 , xleft = 0, xright = 100,
#                    ybottom = 0, ytop = 100)
#       
#     }
#     
#     dev.off()
#     
#     
#     
#     
#   }
#   
# }