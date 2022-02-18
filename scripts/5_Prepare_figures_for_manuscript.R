# --- create figures for manuscript --- #

# clear environment ####
rm(list = ls())

# load libraries ####
library(png)
library(gridExtra)
library(ggplot2)
library(grid)
library(tiff)

# prepare site list and order by average temperature ####
sites <- list.dirs("results/log_core_measurement", full.names = F, recursive = F)
sites <- sites[!sites%in%"Hansley"]

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


# prepare function ####

# g_legend <- function(x = "pt"){
#   p <- get(x, temp_env)
#   a.gplot <- ggplot(data = p, aes(x = varying_x, y = expfit)) + geom_line(aes(group = species, col = species)) +  geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = species), alpha=0.25)+
#     scale_color_manual(values = species_colors) + scale_fill_manual(values = species_colors) + theme(legend.title=element_blank(),  legend.background = element_blank(), legend.box.background =element_blank())+ guides(col=guide_legend(ncol=2), bg=guide_legend(ncol=2))
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)
# }


# DBH response at each sites and for each response ####
what_to_show <- c("log_core_measurement_dbh" = expression(RW~(mm)), "log_BAI_dbh" = expression(BAI~(cm^2)), "log_agb_inc_dbh" = expression(Delta*AGB~(kg)))

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
    
    
    # get the species colors
    species_colors <- get("species_colors", temp_env)
  } # for what in ...
  
  # add a plot for the legend
  # all_plots[[paste0(site, "leg")]] <- all_legends[[site]]#g_legend()
} # for site in sites


all_legends_left <- all_legends[sites_with_dbh[seq(1, length(sites_with_dbh), by = 2)]]
all_legends_right <- all_legends[sites_with_dbh[seq(2, length(sites_with_dbh), by = 2)]]

leg_lengths_left <- unlist(lapply(all_legends_left, function(x) nrow(x$grobs[[1]])))
leg_lengths_right <- unlist(lapply(all_legends_right, function(x) nrow(x$grobs[[1]])))

layout_matrix_left <- matrix(rep( 1:length(leg_lengths_left), leg_lengths_left), ncol = 1)
layout_matrix_right <- matrix(rep( 1:length(leg_lengths_right), leg_lengths_right), ncol = 1)


# png("doc/manuscript/tables_figures/species_legends.png", width =4, height = 10, res = 300, units = "in")
# grid.arrange(
#              arrangeGrob( grobs =all_legends_left, layout_matrix = layout_matrix_left),
#              arrangeGrob( grobs =all_legends_right, layout_matrix = layout_matrix_right), ncol = 2)
# dev.off()


png("doc/manuscript/tables_figures/DBH_responses.png", width = 10, height = 10, res = 300, units = "in")

grid.arrange(arrangeGrob(grobs = all_plots, ncol = 3, vp= grid::viewport(width=0.95, height=0.95)), 
arrangeGrob( grobs =all_legends_left, layout_matrix = layout_matrix_left, size="first"),
arrangeGrob( grobs =all_legends_right, layout_matrix = layout_matrix_right), ncol = 3,  widths = c(4,1,1))

grid::grid.text(sites_abb[sites_with_dbh], x = unit(0.025, "npc"), y = unit(rev(cumsum(c(.95/n_sites/2, rep(.95/n_sites, n_sites-1)))) + .035, "npc"), rot = 90)

grid::grid.text(what_to_show, x = unit(cumsum(c(.05 +.9/4.5/2, rep(.9/4.5, 2))), "npc"), y = unit(.99,  "npc"))

grid::grid.text("DBH (cm)", x = unit(.35, "npc"), y = unit(0.015,  "npc"))

dev.off()

# Year response at each sites and for each response ####
what_to_show <- c("log_core_measurement_dbh" = expression(RW~(mm)), "log_BAI_dbh" = expression(BAI~(cm^2)), "log_agb_inc_dbh" = expression(Delta*AGB~(kg)))

n_sites <- length(sites_with_dbh)

all_plots <- list()
for(site in sites_with_dbh){
  for(what in names(what_to_show)) {
    temp_env <- new.env()
    load(paste0('results/with_Year/', what, "/", site, "/env.RData"), envir = temp_env) 
    
    p <- get("p_Year", envir = temp_env)
    
    # remove title and xlab of p
    p$labels$title <- NULL
    p$labels$x <- NULL
    p$labels$xintercept <- NULL
    
    # change ylim to scale across sites 
    p <- p + ylim(range(get("ylim_p", temp_env)) * ifelse(what == "log_agb_inc_dbh", 1000, 1)) # ylim(ylim_p[[what]][["dbh"]])
    
    # if p is AGB, convert to kg
    if(what == "log_agb_inc_dbh") p$data[c("expfit", "lwr", "upr")] <-   p$data[c("expfit", "lwr", "upr")]*1000
    
    # save into all_plots
    all_plots[[paste0(site, what)]] <- p
    
    
    # get the species colors
    species_colors <- get("species_colors", temp_env)
  } # for what in ...
  
  # add a plot for the legend
  # all_plots[[paste0(site, "leg")]] <- all_legends[[site]]#g_legend()
} # for site in sites

png("doc/manuscript/tables_figures/Year_responses.png", width = 8, height = 10, res = 300, units = "in")

grid.arrange(arrangeGrob(grobs = all_plots, ncol = 3, vp= grid::viewport(width=0.95, height=0.95)))


grid::grid.text(sites_abb[sites_with_dbh], x = unit(0.025, "npc"), y = unit(rev(cumsum(c(.95/n_sites/2, rep(.95/n_sites, n_sites-1)))) + .035, "npc"), rot = 90)

grid::grid.text(what_to_show, x = unit(cumsum(c(.05 +.9/2.9/2, rep(.9/2.9, 2))), "npc"), y = unit(.99,  "npc"))

grid::grid.text("Year", x = unit(.5, "npc"), y = unit(0.015,  "npc"))

dev.off()

# Year response only for log_BAI_dbh response ####
what_to_show <- c("log_BAI_dbh" = expression(BAI~(cm^2)))

n_sites <- length(sites_with_dbh)

all_plots <- list()
xlim_p <- NULL
for(site in sites_with_dbh){
  for(what in names(what_to_show)) {
    temp_env <- new.env()
    load(paste0('results/with_Year/', what, "/", site, "/env.RData"), envir = temp_env) 
    
    p <- get("p_Year", envir = temp_env)
    
    # remove title and xlab of p
    p$labels$title <- NULL
    p$labels$x <- NULL
    p$labels$xintercept <- NULL
    
    # change ylim to scale across sites 
    p <- p + ylim(range(get("ylim_p", temp_env)) * ifelse(what == "log_agb_inc_dbh", 1000, 1)) # ylim(ylim_p[[what]][["dbh"]])
    
    # if p is AGB, convert to kg
    if(what == "log_agb_inc_dbh") p$data[c("expfit", "lwr", "upr")] <-   p$data[c("expfit", "lwr", "upr")]*1000
    
    # save into all_plots
    all_plots[[paste0(site, what)]] <- p
    
    # save x-range
    xlim_p[[paste0(site, what)]] <-  range(p$data$Year)
    
    # get the species colors
    species_colors <- get("species_colors", temp_env)
  } # for what in ...
  
  # add a plot for the legend

  # all_plots[[paste0(site, "leg")]] <- g_legend()

} # for site in sites

png("doc/manuscript/tables_figures/Year_responses_BAI_only.png", width = 10, height = 10, res = 300, units = "in")

grid.arrange(do.call(arrangeGrob, c(lapply(all_plots, function(x) if(!is.null(x$data)) x + xlim(range(xlim_p)) else x), ncol = 2)), vp= grid::viewport(width=0.98, height=0.98), 
             arrangeGrob( grobs =all_legends_left, layout_matrix = layout_matrix_left, size="first"),
             arrangeGrob( grobs =all_legends_right, layout_matrix = layout_matrix_right), ncol = 3,  widths = c(4,1.2,1))


grid::grid.text(sites_abb[sites_with_dbh], x = unit(c(0.2,0.5), "npc"), y = unit(rep(rev(cumsum(c(.95/(n_sites/2), rep(.95/(n_sites/2), (n_sites/2)-1)))),each = 2), "npc"))

grid::grid.text(what_to_show, x = unit(0.025, "npc"), y = unit(rev(cumsum(c(.95/(n_sites/2), rep(.95/(n_sites/2), (n_sites/2)-1))))-0.05, "npc"), rot = 90)

grid::grid.text("Year", x = unit(.35, "npc"), y = unit(0.015,  "npc"))

dev.off()

# Pre and Temp groups ####
what = c("log_core_measurement" = expression(RW~(mm)))

n_sites <- length(sites)

all_plots <- list()
for(site in sites){
 
  
  temp_env <- new.env()
  load(paste0('results/', names(what), "/", site, "/env.RData"), envir = temp_env) 
  
  existing_plots <- ls(temp_env)[grepl("^p_[^in]", ls(temp_env))]
  
  # get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
  clim_var_group <- get("clim_var_group"  , temp_env)
  
  existing_plots <- existing_plots[match(c(1,2), sapply(gsub("p_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))]
  
  # sandardize variable names
  lapply(existing_plots[!is.na(existing_plots)], function(x) {
    p <- get(x, temp_env)
    
    ylim_p <- get("ylim_p", temp_env)
    ylim_p <- ylim_p[!names(ylim_p) %in% "dbh"]
    
    p <- p + ylim(range(ylim_p))
    
    p$labels$x <- eval(parse(text = gsub(" |  ", "~", gsub("-1", "\\^-1", paste0(gsub(substr(p$labels$x, 1, 4), v_names[substr(p$labels$x, 1, 3)], p$labels$x), ")")))))
    p$theme$plot.background <-element_blank()
    assign(x, p, temp_env)
  })
  
  
  # get the species colors
  species_colors <- get("species_colors", temp_env)
  
  # add legend
  
  # assign("leg", g_legend(), envir = temp_env)
  # existing_plots <- c(existing_plots, "leg")
  
  all_plots[[paste0(site, what)]] <- grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.text(label = "no significant\nmain effect") else get(x, temp_env)}), ncol = 2)))
  
} # for(site in sites)

png("doc/manuscript/tables_figures/pre_temp_groups.png", width = 8.2, height = 8.2, res = 300, units = "in")

grid.arrange(grobs = all_plots, vp= grid::viewport(width=0.95, height=0.95), ncol = 2)


grid::grid.text(sites_abb[sites], x = unit(rep(c(0.032, 0.515), 2), "npc"), y = unit(rep(rev(cumsum(c(.95/(n_sites/2), rep(.95/(n_sites/2), n_sites/2-1)))+0.02), each = 2), "npc"))
grid::grid.text(what, x = unit(0.0265, "npc"), y = unit(rev(cumsum(c(.95/(n_sites/2), rep(.95/(n_sites/2), n_sites/2-1)))-0.05), "npc"), rot = 90)

grid::grid.text(c("Precipitation group", "Temperature group"), x = unit(cumsum(c(.05 +.9/3.9/2, rep(.9/3.9, 3))), "npc"), y = unit(.99,  "npc"))


dev.off()

# Pre and Temp groups with DBH interaction ####
what_to_show = c("log_core_measurement_dbh" = expression(RW~(mm)),
                 "log_BAI_dbh" = expression(BAI~(cm^2)), 
                 "log_agb_inc_dbh" = expression(Delta*AGB~(kg)))

n_sites <- length(sites)

for(what in names(what_to_show)) {
  
  all_plots <- list()
  for(site in sites){
    
    
    temp_env <- new.env()
    load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
    
    existing_plots <- ls(temp_env)[grepl("^p_int_", ls(temp_env))]
    
    # get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
    clim_var_group <- get("clim_var_group"  , temp_env)
    
    existing_plots <- existing_plots[match(c(1,2), sapply(gsub("p_int_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))]
    
    if(any(!is.na(existing_plots))) {
      # sandardize variable names
      lapply(existing_plots[!is.na(existing_plots)], function(x) {
        p <- get(x, temp_env) 
        
        ylim_p_int <- get("ylim_p_int", temp_env)
        ylim_p_int <- ylim_p_int[!names(ylim_p_int) %in% "dbh"]
   
        p <- p + ylim(range(ylim_p_int))
        p$labels$x <- eval(parse(text = gsub(" |  ", "~", gsub("-1", "\\^-1", paste0(gsub(substr(p$labels$x, 1, 4), v_names[substr(p$labels$x, 1, 3)], p$labels$x), ")")))))
        p$theme$plot.background <- element_blank()
        p$theme$legend.position <- "none"
        
        assign(x, p, temp_env)
        print(p)
      })
    }
    
    # get the species colors
    # species_colors <- get("species_colors", temp_env)
    
    # add legend
    
    # assign("leg", g_legend(), envir = temp_env)
    # existing_plots <- c(existing_plots, "leg")
    
    all_plots[[paste0(site, what)]] <- grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.text(label = "no significant\ninteractions")  else get(x, temp_env)}), ncol = 2)))
    
  } # for(site in sites)
  
  png(paste0("doc/manuscript/tables_figures/pre_temp_groups_dbh_interactions_", what, ".png"), width = 8.2, height = 8.2, res = 300, units = "in")
  
  grid.arrange(grobs = all_plots, vp= grid::viewport(width=0.95, height=0.95), ncol = 2)
  
  
  grid::grid.text(sites_abb[sites], x = unit(rep(c(0.032, 0.515), 2), "npc"), y = unit(rep(rev(cumsum(c(.95/(n_sites/2), rep(.95/(n_sites/2), n_sites/2-1)))+0.02), each = 2), "npc"))
  grid::grid.text(what_to_show[what], x = unit(0.0265, "npc"), y = unit(rev(cumsum(c(.95/(n_sites/2), rep(.95/(n_sites/2), n_sites/2-1)))-0.05), "npc"), rot = 90)
  
  grid::grid.text(c("Precipitation group", "Temperature group"), x = unit(cumsum(c(.05 +.9/3.9/2, rep(.9/3.9, 3))), "npc"), y = unit(.99,  "npc"))
  
  
  dev.off()
}


# Pre and Temp groups with DBH interaction show case ####
what_to_show = c("log_core_measurement_dbh" = expression(RW~(mm)))
sites_species_to_show = list(HKK = c("TOCI", "Toona ciliata"), 
                          LillyDickey = c("LITU", "Lirodendron tulipifera"),
                          CedarBreaks = c("PIPU", "Picea pungens"))

what = "log_core_measurement_dbh"

all_plots <- list()
rm(leg)
for(site in names(sites_species_to_show)){
  
  
  temp_env <- new.env()
  load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
  
  existing_plots <- ls(temp_env)[grepl("^p_int_", ls(temp_env))]
  
  # get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
  clim_var_group <- get("clim_var_group"  , temp_env)
  
  existing_plots <- existing_plots[match(c(1,2), sapply(gsub("p_int_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))]
  
  if(any(!is.na(existing_plots))) {
    # sandardize variable names
    lapply(existing_plots[!is.na(existing_plots)], function(x) {
      p <- get(x, temp_env) 
      
      ylim_p_int <- get("ylim_p_int", temp_env)
      ylim_p_int <- ylim_p_int[!names(ylim_p_int) %in% "dbh"]
      
      p <- p + ylim(range(ylim_p_int))
      
      p$labels$x <- eval(parse(text = gsub(" |  ", "~", gsub("-1", "\\^-1", paste0(gsub(substr(p$labels$x, 1, 4), v_names[substr(p$labels$x, 1, 3)], p$labels$x), ")")))))
      p$theme$plot.background <- element_blank()
      
      # keep only species we want
      p$data <- p$data[p$data$species_code %in% sites_species_to_show[[site]][1],]
      
      # remove species legend
      p <- p + guides(fill = FALSE, colour = FALSE)
      
      # remove title legend
      p <- p + theme(legend.title = element_blank())
      
      # change legend labels
      p <- p + scale_linetype_manual(values = c("solid", "dotted"), labels=c("max DBH" ,"min DBH"))
      
   
      
      # save legend
      leg <<- ggplotGrob(p)$grobs[[15]]
  
      #remove legend
      p$theme$legend.position <- "none"
      
      #save
      assign(x, p, temp_env)
    })
  }
  

  all_plots[[paste0(site, what)]] <- grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.text(label = "no significant\nmain effect")  else get(x, temp_env)}), ncol = 2)))
  
} # for(site in sites)

png(paste0("doc/manuscript/tables_figures/pre_temp_groups_dbh_interactions.png"), width = 8.2, height = 8.2, res = 300, units = "in")

grid.arrange(arrangeGrob(grobs = all_plots, vp= grid::viewport(width=0.95, height=0.95, y = 0.48), ncol = 1), arrangeGrob(leg), widths = c(8,1))


grid::grid.text(sites_abb[ names(sites_species_to_show)], x = unit(0.06, "npc"), y = unit(c(0.95, 0.63, 0.32), "npc"), just = "right")

grid::grid.text(sapply(sites_species_to_show, "[[", 2), x = unit(0.085, "npc"), y =  unit(c(0.95, 0.63, 0.32), "npc"), gp = gpar(fontface = "italic"), just = "left")

grid::grid.text(what_to_show[what], x = unit(0.0265, "npc"), y = unit(c(0.95, 0.63, 0.32)-0.12, "npc"), rot = 90)

grid::grid.text(c("Precipitation group", "Temperature group"), x = unit(c(0.25,0.65), "npc"), y = unit(.98,  "npc"))


dev.off()





# Pre and Temp groups with DBH interaction show case FOR NIDHI's PAPER####
what_to_show = c("log_core_measurement_dbh" = expression(RW~(mm)))
sites_species_to_show = list(HKK = c("TOCI", "Toona ciliata"), 
                             LillyDickey = c("LITU", "Lirodendron tulipifera"),
                             CedarBreaks = c("PIPU", "Picea pungens"))

what = "log_core_measurement_dbh"

all_plots <- list()
# rm(leg)
for(site in names(sites_species_to_show)){
  
  
  temp_env <- new.env()
  load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) 
  
  existing_plots <- ls(temp_env)[grepl("^p_int_", ls(temp_env))]
  
  # get variable in order or Precipitation, Temperature and cloud groups (but complicated because, pet is in both temp and dtr),....
  clim_var_group <- get("clim_var_group"  , temp_env)[2]
  
  existing_plots <- existing_plots[match(c(1,2), sapply(gsub("p_int_", "", existing_plots), function(x) grep(x,  clim_var_group)[1]))]
  
  if(any(!is.na(existing_plots))) {
    # sandardize variable names
    lapply(existing_plots[!is.na(existing_plots)], function(x) {
      p <- get(x, temp_env) 
      
      ylim_p_int <- get("ylim_p_int", temp_env)
      ylim_p_int <- ylim_p_int[!names(ylim_p_int) %in% "dbh"]
      
      # p <- p + ylim(range(ylim_p_int))
      
      p$labels$x <- eval(parse(text = gsub(" |  ", "~", gsub("-1", "\\^-1", paste0(gsub(substr(p$labels$x, 1, 4), v_names[substr(p$labels$x, 1, 3)], p$labels$x), ")")))))
      p$theme$plot.background <- element_blank()
      
      # keep only species we want
      p$data <- p$data[p$data$species_code %in% sites_species_to_show[[site]][1],]
      
      # remove species legend
      # p <- p + guides(fill = FALSE, colour = FALSE)
      
      # remove title legend
      # p <- p + theme(legend.title = element_blank())
      
      # remove vertical lines and shading
      p$layers[c(1,2)] <- NULL
      
   
      
      # change color to be based on DBH rather than species
      p$layers[[1]]$mapping$colour <- p$layers[[1]]$mapping$linetype
      p$layers[[2]]$mapping$fill <- p$layers[[1]]$mapping$colour 
      
      p <- p +  scale_color_manual(values = c("#EFA729", "#0466AC"), labels=c("max DBH" ,"min DBH"))
      p <- p +  scale_fill_manual(values = c("#EFA729", "#0466AC"), labels=c("max DBH" ,"min DBH"))
      
      # remove DBH as linetype
      p$layers[[1]]$mapping$linetype <- NULL
      
      # make line thicker and alpha less transparent
      p$layers[[1]]$aes_params$size = 1.2
      p$layers[[2]]$aes_params$alpha = 0.5
      # redo the plot to get legend (not sure why I have to do this)
  
      
      a <- ggplot(p$data[p$data$species_code == p$data$species_code[1],], aes(x = varying_x, y = expfit, group = DBH)) + 
        geom_line(aes(color = DBH), size=  1.2) + geom_ribbon(aes(ymin = lwr       , ymax = upr, fill = DBH), alpha = 0.5)  +  
        scale_color_manual(values = c("#EFA729", "#0466AC"), labels=c("max DBH" ,"min DBH")) +  
        scale_fill_manual(values = c("#EFA729", "#0466AC"), labels=c("max DBH" ,"min DBH")) + 
        guides(color = guide_legend(nrow = 1)) + 
        theme(legend.title = element_blank())
      
      # save legend
      leg <<- ggplotGrob(a)$grobs[[15]]
      
      #remove legend
      # p$theme$legend.position <- "none"
      
      #save
      assign(x, p, temp_env)
    })
  }
  
  
  all_plots[[paste0(site, what)]] <- get(existing_plots, temp_env) # grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.text(label = "no significant\nmain effect")  else get(x, temp_env)}), ncol = 2)))
  
} # for(site in sites)

png(paste0("doc/manuscript/tables_figures/pre_temp_groups_dbh_interactions_for_NIDHI.png"), width = 4.1, height = 8.2, res = 300, units = "in")

# grid.arrange(arrangeGrob(leg), arrangeGrob(grobs = all_plots, vp= grid::viewport(width=0.87, height=0.95, y = 0.48, x = 0.50), ncol = 1, size = "first"), heights  = c(.5, 8))

grid.arrange(arrangeGrob(leg), arrangeGrob(do.call(rbind, sapply( all_plots, ggplotGrob))), heights  = c(.5, 8), vp= grid::viewport(width=0.87, height=0.95, y = 0.48, x = 0.50))

grid::grid.text(c("(Huai Kha Khaeng, Tailand)", "(Lilly Dickey Woods, Indiana, USA)", "(Cedar Breaks, Utah, USA)"), x = unit(0.9, "npc"), y = unit(c(0.89, 0.58, 0.28)-0.025, "npc"), just = "right")

grid::grid.text(sapply(sites_species_to_show, "[[", 2), x = unit(0.9, "npc"), y =  unit(c(0.89, 0.58, 0.28), "npc"), gp = gpar(fontface = "italic"), just = "right")

grid::grid.text(expression(Ring~width ~ (mm)), x = unit(0.04, "npc"), y = unit( 0.5, "npc"), rot = 90,  gp=gpar(fontsize=18, fontface = "bold"))

# grid::grid.text(c("Precipitation group", "Temperature group"), x = unit(c(0.25,0.65), "npc"), y = unit(.98,  "npc"))


dev.off()




# comparison with quilt ####
site = "SCBI"
v = "pet"
what = "log_core_measurement"

png("doc/manuscript/tables_figures/quilt_comparison.png", width = 10, height = 5, units = "in", res = 300)

layout(matrix(c(1, 1, 2, 3, 5, 5, 4, 6, 6, 7, 7, 7), nrow = 3), heights = c(3,2,1), widths = c(2.2, 1, 1, 2))
par(oma = c(1,0,2,0))

## a) quilt ####
img <- readTIFF(paste0("C:/Users/HerrmannV/Dropbox (Smithsonian)/GitHub/climate_sensitivity_cores/results/1901_2009/figures/monthly_correlation/CRU_", site, "_1901_2016/", v, ".tif"), as.is = T)

img1 <- img[30:(nrow(img)-150), 1:800, ] # keep only quilt half of plot
img2 <- img[300:(nrow(img)-150),800:(ncol(img)-60), ] # keep only legend


par(mar = c(3,0,0,0), mgp = c(3,0,0))
plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")

rasterImage( img1 , xleft = 0, xright = 100,
             ybottom = 0, ytop = 100)

rect(xleft = rev(seq(27, 94.5, length.out = 17))[5]+1,
     ybottom = 0,
     xright = rev(seq(27, 94.5, length.out = 17))[2],
     ytop = 85, 
     lwd = 2) # add rectangle on May-Jul

# segments(x0 = c(28, 63), y0 = rep(97, 2), x1 = c(62, 93), y1 = rep(97, 2), lwd = 2, col = c("grey", "black"))
# text(x = c(mean(c(28, 62)), mean(c(63, 93))),
#      y = 100,
#      labels = c("previous year", "current year"),
#      col = c("grey", "black"))

mtext("a)", side = 3, adj = 0.1, line = -2, cex = .8)
axis(1, at = seq(27, 94.5, length.out = 17)+.2, labels = rev(c(0:15, "")), cex.axis = .65, line = -1, col.ticks = "white", hadj = 1.5, tcl=-.1
      )
mtext("months prior to current August", side = 1, cex = .6, adj = 0.7, line = 0)

title(expression(bold(underline("Traditional analysis"))), xpd = NA, line = .5, adj  = .7)

par(mar = c(0,5,1,2))
plot(0:100, 0:100, type = "n", axes = F, xlab = "", ylab = "")

rasterImage( matrix(as.vector(as.raster(img2)), ncol = nrow(img2))[ncol(img2):1,], xleft = 10, xright = 110,
             ybottom = 0, ytop = 100)
text(x = -5, y = 50, labels = expression(bold("Correlation")), xpd = NA)
legend(x = -18, y = 180, pch = c(21, 24), legend = c("0.05", "0.0002"), xpd = NA, bty = "n", title = expression(bold("Significance")))

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
p <- p + labs(y = expression(RW~(mm)))
p <- p + theme(legend.position="right", legend.text = element_text(size = 8)) # add legend
p$layers[c(1,2)] <- NULL # remove vertical line and shading
p$labels$x <- eval(parse(text = gsub(" |  ", "~", gsub("-1", "\\^-1", paste0(gsub(substr(p$labels$x, 1, 4), v_names[substr(p$labels$x, 1, 3)], p$labels$x), ")"))))) # change variable label
p <- p + scale_fill_discrete(labels = gsub(" \\(\\d*\\)", "", levels(p$data$species))) + scale_colour_discrete(labels = gsub(" \\(\\d*\\)", "", levels(p$data$species)))
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
what_to_show <- c("log_core_measurement_dbh" = expression(RW~(mm)), "log_BAI_dbh" = expression(BAI~(cm^2)), "log_agb_inc_dbh" = expression(Delta*AGB~(kg)))

sites_to_show_case <- c("SCBI", "NewMexico")


for(with_Year in c(FALSE, TRUE)) {
  show_case <- list()
  for(site in sites_with_dbh){
    all_plots <- list()

    for( what in names(what_to_show)) {
      
      temp_env <- new.env()
      if(with_Year) load(paste0('results/with_Year/', what, "/", site, "/env.RData"), envir = temp_env) # to get the Year plot
      if(!with_Year) load(paste0('results/', what, "/", site, "/env.RData"), envir = temp_env) # to get all other plots
      
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
      
      
      # get the species colors
      species_colors <- get("species_colors", temp_env)
      
      # get the legend
      
      # assign("leg", g_legend(), envir = temp_env)
      # existing_plots <- c(existing_plots, "leg")
      
      
      all_plots[[what]] <- do.call(cowplot::plot_grid, c(lapply(existing_plots, function(x)  {if(is.na(x)) grid.text(label = "no significant\nmain effect") else get(x, temp_env)}), ncol = ifelse(with_Year, 4, 3), list(align = "hv")))

      if (site %in% sites_to_show_case) show_case[[paste0(site, what)]] <-  all_plots[[what]]  #grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots[-4], function(x)  {if(is.na(x)) grid.rect(gp=gpar(col="white")) else get(x, temp_env)}), ncol = 3)))
    }  
    
    png(paste0("results/composite_plots/", site, ifelse(with_Year, "_with_Year", ""), ".png"), width = 10, height = 10, res = 300, units = "in")
    
    grid.arrange(
      arrangeGrob(grobs = all_plots, ncol = 1), 
      arrangeGrob( all_legends[[site]]),
      ncol = 2,  widths = c(4,1))
    
    # grid.arrange(grobs = all_plots, vp= grid::viewport(width=0.95, height=0.95), ncol = 1)
    
    grid::grid.text(what_to_show,  x = unit(0.01, "npc"), y = unit(rev(cumsum(c(1/length(what_to_show)/2, rep(1/length(what_to_show), length(what_to_show)-1)))), "npc"), rot = 90)
    
    dev.off()
    
    
    # if(site %in% sites_to_show_case){
    #   # save the plot
    #   show_case[[paste0(site, "leg")]] <- get("leg", temp_env)
    # }
  }



# png(paste0("doc/manuscript/tables_figures/show_case_response_plots", ifelse(with_Year, "_with_Year", ""), ".png"), width = 10, height = 8, res = 300, units = "in")

grid.arrange(
  arrangeGrob(grobs = show_case, layout_matrix = matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)), 
  arrangeGrob( grobs = all_legends_2_cols[sites_to_show_case], ncol = 2),
  ncol = 1,  heights = c(2,1))


# grid.arrange(grobs = show_case, vp= grid::viewport(width=0.95, height=0.95), layout_matrix = matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4))

grid::grid.text(what_to_show,  x = unit(0.015, "npc"), y = unit((rev(cumsum(c(1.3/length(what_to_show)/2, rep(1.3/2/length(what_to_show), length(what_to_show)-1)))))+0.25, "npc"), rot = 90)

grid::grid.text(sites_abb[sites_to_show_case[order(match(sites_to_show_case, sites))]], x = unit(cumsum(c(.05 +.9/2/2, rep(.9/2, 1))), "npc"), y = unit(.99,  "npc"))
# dev.off()

} # for(with_Year in c(FALSE, TRUE))