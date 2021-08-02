# give site abbrevationtion in paper
sites_abb <- list(BCI  = "BCNM",
                  HKK = "HKK",
                  NewMexico = "LT",
                  CedarBreaks = "CB",
                  SCBI = "SCBI",
                  LillyDickey = "LDW",
                  HarvardForest = "HF",
                  # Nebraska = "NE",
                  Niobara = "NIO",
                  Hansley = "NE",
                  Zofin = "ZOF",
                  ScottyCreek = "SC")


# order sites by average temperature
MAT_order <- read.csv("doc/manuscript/tables_figures/sites.csv")

sites <- names(sites_abb[match(MAT_order[,1], sites_abb)])

# load all legends ####
load("results/all_legends.Rdata")

# sites with dbh
sites_with_dbh <- sites #[-grep("CedarBreaks", sites)]

# standardize variable names ####
v_names <- list(tmn = "expression(T[min]~",
                tmx = "expression(T[max]~",
                tmp = "expression(T[mean]~",
                pet = "expression(PET~",
                pre = "expression(PPT~",
                wet = "expression(PDF~",
                dbh = "expression(DBH~")

with_Year = TRUE
site = "NewMexico"
what = "log_core_measurement_dbh"
temp_env <- new.env()
all_plots <- list()

load(paste0('results/with_Year/', what, "/", site, "/env.RData"), envir = temp_env) # to get the Year plot

existing_plots <- ls(temp_env)[grepl("^p_", ls(temp_env))]

# get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
clim_var_group <- get("clim_var_group"  , temp_env)

if(with_Year) existing_plots <- existing_plots[c(1, length(existing_plots), na.omit(match(c(1,2), sapply(gsub("p_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))))]
if(!with_Year) existing_plots <- existing_plots[c(1, na.omit(match(c(1,2), sapply(gsub("p_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))))]

# change ylim to scale across sites 
lapply(existing_plots, function(x) assign(x, get(x, temp_env) + ylim(range(get("ylim_p", temp_env))* ifelse(what == "log_agb_inc_dbh", 1000, 1)), envir = temp_env))

# if p is AGB, convert to kg
if(what == "log_agb_inc_dbh")   lapply(existing_plots, function(x) {
  p <- get(x, temp_env)
  p$data[c("expfit", "lwr", "upr")] <-   p$data[c("expfit", "lwr", "upr")]*1000
  p$layers[[1]]$aes_params$ymin  <- p$layers[[1]]$aes_params$ymin * 1000
  p$layers[[1]]$aes_params$ymax  <- p$layers[[1]]$aes_params$ymax * 1000
  
  assign(x, p, envir = temp_env)
})

# standardize variable names
lapply(existing_plots[switch(as.character(with_Year), "TRUE" = -2, "FALSE" = c(1:length(existing_plots)))], function(x) { # -2 is to not do it for Year
  p <- get(x, temp_env)
  p$labels$x <- eval(parse(text = gsub(" |  ", "~", gsub("-1", "\\^-1", paste0(gsub(substr(p$labels$x, 1, 4), v_names[tolower(substr(p$labels$x, 1, 3))], p$labels$x), ")")))))
  assign(x, p, temp_env)
})

# 
# # keep only PIPO
# p$layers[[1]]$aes_params$linetype <- p$layers[[1]]$aes_params$linetype[p$data$species_code == "PIPO"]
# p$data <- p$data[p$data$species_code == "PIPO", ]
# p
# 
# # add data points 
# 
temp_all_env <- new.env()
load(paste0('results/with_Year/', what, "/", site, "_all_env.RData"), envir = temp_all_env) # to get the Year plot
ls(temp_all_env)

dp <-  get("PIPO_best_model", temp_all_env)$data
# p + geom_point(data = dp, mapping = aes(x = Year, y = log_core_measurement))

lapply(existing_plots, function(x) { # -2 is to not do it for Year
  p <- get(x, temp_env)
  
  # keep only PIPO
  # p$layers[[grep("linetype", sapply(lapply(p$layers, function(x) x$aes_params), names))]]$aes_params$linetype <- p$layers[[grep("linetype", sapply(lapply(p$layers, function(x) x$aes_params), names))]]$aes_params$linetype[p$data$species_code == "PIPO"]
  p$data <- p$data[p$data$species_code == "PIPO", ]
  # p
  
  # add data points 
  
  dp$x <- dp[, gsub("p_", "",  x)]
  
  p <- p + geom_point(data = dp, mapping = aes(x = x, y = exp(log_core_measurement)), #, col = as.numeric(dp$coreID) , 
                      alpha= 0.1)+ ylim(range(c(p$data$upr, exp(dp$log_core_measurement))))
  
  # p <- p + scale_y_continuous(trans = "log", labels = scales::scientific)
  
  # remove rectangle
  p$layers[grepl("Rect|Vline", sapply(p$layers, function(x) class(x$geom)[[1]]))] <- NULL
  
  # reorder layers
  p$layers <- p$layers[rev(1:length(p$layers))]
  
  
  # save
  assign(x, p, temp_env)
})



all_plots[[what]] <- do.call(cowplot::plot_grid, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.text(label = "no significant\nmain effect") else get(x, temp_env)}), ncol = ifelse(with_Year, 4, 3), list(align = "hv")))






png(paste0("doc/manuscript/tables_figures/schematic_figure_bottom_part.png"), width = 10, height = 3, res = 300, units = "in")
all_plots[[what]] 
grid::grid.text("RW (mm)",  x = unit(0.015, "npc"), y = 0.5, rot = 90)
dev.off()
# get the legend ####
# all_Biol <- get("all_Biol", temp_all_env)
# species_colors <- get("species_colors", temp_all_env)
# names(species_colors) <- substr(names(species_colors), 1,4)
# 
# x <- all_Biol[[site]]
# x <- x[!is.na(x$dbh),]
# # x <- x[x$species_code %in% summary_data$species_code[summary_data$site %in% site],]
# x <- droplevels(x[!duplicated(x$species_code) & x$species_code %in% "PIPO",])
# 
# a.gplot <- ggplot(x, aes(x = Year, y = dbh)) + geom_line(aes(group = paste(genus, species), col = paste(genus, species))) +  geom_ribbon(aes(ymin=min(dbh), ymax=max(dbh), bg = paste(genus, species)), alpha=0.25) + labs(col = sites_abb[[site]], bg = sites_abb[[site]])+ theme(legend.background = element_blank(), legend.box.background =element_blank(),  legend.justification = "left", legend.margin=margin(0,0,0,0), legend.box.margin=margin(0,0,0,10)) 
# 
#  leg <- ggplot_gtable(ggplot_build(a.gplot))$grobs[[which(sapply( ggplot_gtable(ggplot_build(a.gplot))$grobs, function(x) x$name) == "guide-box")]]
# 
#  
#  # make plot
# grid.arrange(
#   all_plots[[what]] ,   leg,
#   ncol = 2,  widths = c(6,1))
# 
# grid::grid.text("RW (mm)",  x = unit(0.015, "npc"), y = 0.5, rot = 90)
# 

