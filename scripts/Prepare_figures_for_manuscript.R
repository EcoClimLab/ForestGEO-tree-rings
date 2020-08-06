# --- create figures for manuscript --- #

# clear environment ####
rm(list = ls())

# load libraries ####
library(png)
library(gridExtra)
library(ggplot2)
library(grid)

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

# prepare function ####

g_legend <- function(x = "pt"){
  p <- get(x, temp_env)
  a.gplot <- ggplot(data = p, aes(x = varying_x, y = expfit)) + geom_line(aes(group = species, col = species)) +  geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = species), alpha=0.25) + theme(legend.title=element_blank())+ guides(col=guide_legend(ncol=2), bg=guide_legend(ncol=2))
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
# load and calculate y limits ####
load("results/ylims_for_GLS_plots.RData")
ylim_p <- lapply(ylim_p, lapply, range)
ylim_p[["log_agb_inc_dbh"]] <- lapply(ylim_p[["log_agb_inc_dbh"]], function(x) x*1000) # convert agb to kg
# DBH response at each sites and for each response ####
what_to_show <- c("log_core_measurement_dbh" = expression(Delta*r~(mm)), "log_BAI_dbh" = expression(BAI~(cm^2)), "log_agb_inc_dbh" = expression(Delta*AGB~(kg)))



all_plots <- list()
for(site in sites){
  for(what in names(what_to_show)) {
    temp_env <- new.env()
    load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
    
    p <- get("p_dbh", envir = temp_env)
    
    # remove title and xlab of p
    p$labels$title <- NULL
    p$labels$x <- NULL
    
    # change ylim to scale across sites 
    p <- p + ylim(ylim_p[[what]][["dbh"]])
    
    # if p is AGB, convert to kg
    if(what == "log_agb_inc_dbh") p$data[c("expfit", "lwr", "upr")] <-   p$data[c("expfit", "lwr", "upr")]*1000
    
    # save into all_plots
    all_plots[[paste0(site, what)]] <- p
  } # for what in ...
  
  # add a plot for the legend
  all_plots[[paste0(site, "leg")]] <- g_legend()
} # for site in sites

png("doc/manuscript/tables_figures/DBH_responses.png", width = 8, height = 10, res = 300, units = "in")

grid.arrange(arrangeGrob(grobs = all_plots, ncol = 4, vp= grid::viewport(width=0.95, height=0.95)))


grid::grid.text(sites, x = unit(0.01, "npc"), y = unit(rev(cumsum(c(1/n_sites/2, rep(1/n_sites, n_sites-1)))), "npc"), rot = 90)

grid::grid.text(what_to_show, x = unit(cumsum(c(.05 +.9/4/2, rep(.9/4, 2))), "npc"), y = unit(.99,  "npc"))

grid::grid.text("dbh (cm)", x = unit(.5, "npc"), y = unit(0.015,  "npc"))

dev.off()

# Pre, Temp and cld groups ####

all_plots <- list()
for(site in sites){
  what = "log_core_measurement"
  
  temp_env <- new.env()
  load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
  
  existing_plots <- ls(temp_env)[grepl("^p_", ls(temp_env))]
  
  # get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
  clim_var_group <- get("clim_var_group"  , temp_env)
  
  existing_plots <- existing_plots[match(c(1,2,3), sapply(gsub("p_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))]
  
  assign("leg", g_legend(), envir = temp_env)
  existing_plots <- c(existing_plots, "leg")
  
  all_plots[[paste0(site, what)]] <- grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.rect(gp=gpar(col="white")) else get(x, temp_env)}), ncol = 4)))
  
} # for what in ...

png("doc/manuscript/tables_figures/pre_temp_cld_groups.png", width = 8, height = 10, res = 300, units = "in")

grid.arrange(grobs = all_plots, vp= grid::viewport(width=0.95, height=0.95), ncol = 1)


grid::grid.text(sites, x = unit(0.01, "npc"), y = unit(rev(cumsum(c(1/n_sites/2, rep(1/n_sites, n_sites-1)))), "npc"), rot = 90)

grid::grid.text(c("Precipiation group", "Temperature group", "Cloud group"), x = unit(cumsum(c(.05 +.9/4/2, rep(.9/4, 2))), "npc"), y = unit(.99,  "npc"))


dev.off()

# comparison with quilt ####
site = "SCBI"
v = "pre"
what = "log_core_measurement"

png("doc/manuscript/tables_figures/quilt_comparison.png", width = 10, height = 4 , units = "in", res = 300)

layout(matrix(c(1, 1, 2, 3, 5, 5, 4, 6, 6, 7, 7, 7), nrow = 3), heights = c(3,2,1), widths = c(2.2, 1, 1, 2))
par(mar = c(0,0,0,0))

## a) quilt ####
img <- readPNG(paste0("results/traditional_analysis/", site, "/1901_2009/figures/monthly_correlation/", v, ".png"))
img1 <- img[50:(nrow(img)-150), 1:800, ] # keep only quilt half of plot
img2 <- img[,800:(ncol(img)-60), ] # keep only legend


par(mar = c(0,0,0,0))
plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")

rasterImage( img1 , xleft = 0, xright = 100,
             ybottom = 0, ytop = 100)

par(mar = c(0,0,0,0))
plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")

rasterImage( matrix(as.vector(as.raster(img2)), ncol = nrow(img2))[ncol(img2):1,], xleft = 10, xright = 110,
             ybottom = 5, ytop = 100)

## b,c,d,e) climwin ####
img <- readPNG(paste0("results/", what, "/", site, "/climwin_", v, "_quad_14_0_SCBI.png"))

img1 <- img[80:(nrow(img)/2), c((2*ncol(img)/4): (3*ncol(img)/4)) ,] # beta linear
img2 <- img[80:(nrow(img)/2), c((3*ncol(img)/4): (4*ncol(img)/4)) ,] # beta quadratic
img3 <- img[80:(nrow(img)/2), c(1: (1*ncol(img)/4)) ,] # AIC
img4 <- img[c((nrow(img)/2)+1):nrow(img), c((1*ncol(img)/4): (2*ncol(img)/4)) ,] # response

for(i in 1:4) {
  plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
  
  rasterImage( get(paste0("img", i)) , xleft = 0, xright = 100,
               ybottom = 0, ytop = 100)
}

## f) response curves ####
temp_env <- new.env()
load(paste0("results/", what, "/", site, "/env.RData"), temp_env)

p <- get(paste0("p_", v), temp_env)
p <- p + labs(y = expression(Delta*r~(mm)))
plot.new()              ## suggested by @Josh
vps <- gridBase::baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
vp1 <- plotViewport(c(1.8,1,0,1)) ## create new vp with margins, you play with this values 

print(p,vp = vp1) 
## dev.off() ####
dev.off()