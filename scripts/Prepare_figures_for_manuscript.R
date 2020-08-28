# --- create figures for manuscript --- #

# clear environment ####
rm(list = ls())

# load libraries ####
library(png)
library(gridExtra)
library(ggplot2)
library(grid)

# prepare site list and site coordinates ####
sites <- list.dirs("results/log_core_measurement", full.names = F, recursive = F)


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

# sites with dbh
sites_with_dbh <- sites[-grep("CedarBreaks", sites)]

# give site abbrevationtion in paper
sites_abb <- list(BCI  = "BCI",
                  HKK = "HKK",
                  NewMexico = "LT",
                  CedarBreaks = "CB",
                  SCBI = "SCBI",
                  LillyDickey = "LDW",
                  HarvardForest = "HF",
                  NB = "NB",
                  Zofin = "ZOF",
                  ScottyCreek = "SC")

# prepare function ####

g_legend <- function(x = "pt"){
  p <- get(x, temp_env)
  a.gplot <- ggplot(data = p, aes(x = varying_x, y = expfit)) + geom_line(aes(group = species, col = species)) +  geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = species), alpha=0.25) + theme(legend.title=element_blank())+ guides(col=guide_legend(ncol=2), bg=guide_legend(ncol=2))
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# DBH response at each sites and for each response ####
what_to_show <- c("log_core_measurement_dbh" = expression(Delta*r~(mm)), "log_BAI_dbh" = expression(BAI~(cm^2)), "log_agb_inc_dbh" = expression(Delta*AGB~(kg)))

n_sites <- length(sites_with_dbh)

all_plots <- list()
for(site in sites_with_dbh){
  for(what in names(what_to_show)) {
    temp_env <- new.env()
    load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
    
    p <- get("p_dbh", envir = temp_env)
    
    # remove title and xlab of p
    p$labels$title <- NULL
    p$labels$x <- NULL
    
    # change ylim to scale across sites 
    p <- p + ylim(range(get("ylim_p", temp_env)) * ifelse(what == "log_agb_inc_dbh", 1000, 1)) # ylim(ylim_p[[what]][["dbh"]])
    
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


grid::grid.text(sites_abb[sites_with_dbh], x = unit(0.025, "npc"), y = unit(rev(cumsum(c(.95/n_sites/2, rep(.95/n_sites, n_sites-1)))) + .035, "npc"), rot = 90)

grid::grid.text(what_to_show, x = unit(cumsum(c(.05 +.9/4/2, rep(.9/4, 2))), "npc"), y = unit(.99,  "npc"))

grid::grid.text("dbh (cm)", x = unit(.5, "npc"), y = unit(0.015,  "npc"))

dev.off()

# Pre and Temp groups ####


what = "log_core_measurement"

n_sites <- length(sites)

all_plots <- list()
for(site in sites){
 
  
  temp_env <- new.env()
  load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
  
  existing_plots <- ls(temp_env)[grepl("^p_", ls(temp_env))]
  
  # get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
  clim_var_group <- get("clim_var_group"  , temp_env)
  
  existing_plots <- existing_plots[match(c(1,2), sapply(gsub("p_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))]
  
  assign("leg", g_legend(), envir = temp_env)
  existing_plots <- c(existing_plots, "leg")
  
  all_plots[[paste0(site, what)]] <- grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.rect(gp=gpar(col="white")) else get(x, temp_env)}), ncol = 3)))
  
} # for(site in sites)

png("doc/manuscript/tables_figures/pre_temp_groups.png", width = 8, height = 10, res = 300, units = "in")

grid.arrange(grobs = all_plots, vp= grid::viewport(width=0.95, height=0.95), ncol = 1)


grid::grid.text(sites_abb[sites], x = unit(0.025, "npc"), y = unit(rev(cumsum(c(.95/n_sites/2, rep(.95/n_sites, n_sites-1))) + 0.05), "npc"), rot = 90)

grid::grid.text(c("Precipiation group", "Temperature group"), x = unit(cumsum(c(.05 +.9/3/2, rep(.9/3, 1))), "npc"), y = unit(.99,  "npc"))


dev.off()

# comparison with quilt ####
site = "SCBI"
v = "pet"
what = "log_core_measurement"

png("doc/manuscript/tables_figures/quilt_comparison.png", width = 10, height = 5, units = "in", res = 300)

layout(matrix(c(1, 1, 2, 3, 5, 5, 4, 6, 6, 7, 7, 7), nrow = 3), heights = c(3,2,1), widths = c(2.2, 1, 1, 2))
par(oma = c(1,0,2,0))

## a) quilt ####
img <- readPNG(paste0("results/traditional_analysis/", site, "/1901_2009/figures/monthly_correlation/", v, ".png"))
img1 <- img[50:(nrow(img)-150), 1:800, ] # keep only quilt half of plot
img2 <- img[300:(nrow(img)-150),800:(ncol(img)-60), ] # keep only legend


par(mar = c(3,0,0,0), mgp = c(3,0,0))
plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")

rasterImage( img1 , xleft = 0, xright = 100,
             ybottom = 0, ytop = 100)

segments(x0 = c(28, 63), y0 = rep(97, 2), x1 = c(62, 93), y1 = rep(97, 2), lwd = 2, col = c("grey", "black"))
text(x = c(mean(c(28, 62)), mean(c(63, 93))),
     y = 100,
     labels = c("previous year", "current year"),
     col = c("grey", "black"))

mtext("a)", side = 3, adj = 0.1, line = -2, cex = .8)
axis(1, at = seq(27, 94.5, length.out = 17), labels = rev(c(0:15, "")), cex.axis = .65, line = -1, col.ticks = "white", hadj = 1.5, tcl=-.1
      )
mtext("months prior to current August", side = 1, cex = .6, adj = 0.7, line = 0)

title(expression(bold(underline("Traditional analysis"))), xpd = NA, line = .5, adj  = .7)

par(mar = c(0,5,1,2))
plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")

rasterImage( matrix(as.vector(as.raster(img2)), ncol = nrow(img2))[ncol(img2):1,], xleft = 10, xright = 110,
             ybottom = 0, ytop = 100)
text(x = 0, y = 50, labels = expression(bold("Correlation")), xpd = NA)
legend(x = -15, y = 180, pch = c(21, 24), legend = c("0.05", "0.0002"), xpd = NA, bty = "n", title = expression(bold("Significance")))

## b,c,d,e) climwin ####
img <- readPNG(list.files(paste0("results/", what, "/", site), pattern = v, full.names = T))

img1 <- img[80:(nrow(img)/2), c((2*ncol(img)/4): (3*ncol(img)/4)) ,] # beta linear
img2 <- img[80:(nrow(img)/2), c((3*ncol(img)/4): (4*ncol(img)/4)) ,] # beta quadratic
img3 <- img[80:(nrow(img)/2), c(1: (1*ncol(img)/4)) ,] # AIC
img4 <- img[c((nrow(img)/2)+1):nrow(img), c((1*ncol(img)/4+1): (2*ncol(img)/4)) ,] # response

for(i in 1:4) {
  par(mar = c(0,0,0,0))
  plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")
  
  rasterImage( get(paste0("img", i)) , xleft = 0, xright = 100,
               ybottom = 0, ytop = 100)
  
  if(i == 3) {
    rect(0, 90, 100, 100, col = "white", border = "white")
    mtext(expression(Delta*AIC), side = 3, line = -1, cex = .7)
    mtext("(compared to null model)", side = 3, line = -1.7, cex = .7, adj = 1)
  }
  mtext(paste0(letters[i+1], ")"), side = 3, adj = 0.1, line = -2, cex = .8)
  
  
}

title(expression(bold(underline("Variable identification in"~ bolditalic("climwin") ~" step"))), xpd = NA, line = .5, adj = .55, outer = T)

## f) response curves ####
temp_env <- new.env()
load(paste0("results/", what, "/", site, "/env.RData"), temp_env)

p <- get(paste0("p_", v), temp_env)
p <- p + labs(y = expression(Delta*r~(mm)))
p <- p + theme(legend.position="right", legend.text = element_text(size = 8)) # add legend
p$layers[c(1,2)] <- NULL # remove vertical line and shading
plot.new()              ## suggested by @Josh
vps <- gridBase::baseViewports()
pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
vp1 <- plotViewport(c(2,0,2,0)) ## create new vp with margins, you play with this values 

print(p, vp = vp1)
# windowsFonts(Times=windowsFont("TT Times New Roman"))
grid.text(paste0(letters[6], ")"),x = 0.1, y = .90, gp=gpar(fontsize=9.8, fontfamily=""))

title(expression(bold(underline("GLS output"))), xpd = NA, line = .5, adj = .85,  outer = T)
# dev.off()####
dev.off()


# create composite image of all models for each sites + show case 2 sites ####
what_to_show <- c("log_core_measurement_dbh" = expression(Delta*r~(mm)), "log_BAI_dbh" = expression(BAI~(cm^2)), "log_agb_inc_dbh" = expression(Delta*AGB~(kg)))

sites_to_show_case <- c("SCBI", "NewMexico")
show_case <- list()
for(site in sites_with_dbh){
  all_plots <- list()
  for( what in names(what_to_show)) {
  
  temp_env <- new.env()
  load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
  
  existing_plots <- ls(temp_env)[grepl("^p_", ls(temp_env))]
  
  # get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
  clim_var_group <- get("clim_var_group"  , temp_env)
  
  existing_plots <- existing_plots[c(1, match(c(1,2), sapply(gsub("p_", "", existing_plots), function(x) grep(x,  clim_var_group)[1])))]
  
  # change ylim to scale across sites 
  lapply(existing_plots, function(x) assign(x, get(x, temp_env) + ylim(range(get("ylim_p", temp_env))* ifelse(what == "log_agb_inc_dbh", 1000, 1)), envir = temp_env))
 
  # if p is AGB, convert to kg
  if(what == "log_agb_inc_dbh")   lapply(existing_plots, function(x) {
    p <- get(x, temp_env)
    p$data[c("expfit", "lwr", "upr")] <-   p$data[c("expfit", "lwr", "upr")]*1000
    assign(x, p, envir = temp_env)
  })
  
 # get the legend
  
  assign("leg", g_legend(), envir = temp_env)
  existing_plots <- c(existing_plots, "leg")
  
  all_plots[[what]] <- grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.rect(gp=gpar(col="white")) else get(x, temp_env)}), ncol = 4)))
  
  if (site %in% sites_to_show_case) show_case[[paste0(site, what)]] <- grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots[-4], function(x)  {if(is.na(x)) grid.rect(gp=gpar(col="white")) else get(x, temp_env)}), ncol = 3)))
  }  

png(paste0("results/composite_plots/", site, ".png"), width = 8, height = 10, res = 300, units = "in")

grid.arrange(grobs = all_plots, vp= grid::viewport(width=0.95, height=0.95), ncol = 1)

grid::grid.text(what_to_show,  x = unit(0.01, "npc"), y = unit(rev(cumsum(c(1/length(what_to_show)/2, rep(1/length(what_to_show), length(what_to_show)-1)))), "npc"), rot = 90)

dev.off()


if(site %in% sites_to_show_case){
  # save the plot
  show_case[[paste0(site, "leg")]] <- get("leg", temp_env)
}
}

png(paste0("doc/manuscript/tables_figures/show_case_response_plots.png"), width = 10, height = 8, res = 300, units = "in")

grid.arrange(grobs = show_case, vp= grid::viewport(width=0.95, height=0.95), layout_matrix = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4))

grid::grid.text(what_to_show,  x = unit(0.015, "npc"), y = unit(rev(cumsum(c(1/length(what_to_show)/2, rep(1/length(what_to_show), length(what_to_show)-1))))*3/4+1/4, "npc"), rot = 90)

grid::grid.text(sites_abb[sites_to_show_case[order(match(sites_to_show_case, sites))]], x = unit(cumsum(c(.05 +.9/2/2, rep(.9/2, 1))), "npc"), y = unit(.99,  "npc"))
dev.off()
