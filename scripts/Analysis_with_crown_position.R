# implement the plan which is:
# - use climwin to identify the strongest temperature related variable and moisture related variables signals (and time windoe) using mixed model with species and core ID (not nested) as random effects and the individual-levels residuals of core measurements spline.
# - using the best vriable, use a gam on the raw data, with a smoothing on dbh and on, Year at the individual level fo year. On model for each species.
# help about specifying gams https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.models.html


# clear environment ####
rm(list = ls())

# load libraries ####
library(climwin)
library(lme4)
library(mgcv)
library(splines)
library(grid)
library(MuMIn)
library(usdm)


# set parameters ####
core_type <- "CSV"
reference_date <- c(30, 8) # refday in slidingwin
window_range <- c(15, 0) #range in slidingwin

what = "log_core_measurement" # options: "log_core_measurement", "log_agb_inc"

# load data ####
## core data  (processed, which include climate data for best windows) ####

Biol <- read.csv(paste0("processed_data/core_data_with_best_climate_signal", ifelse(what %in% "log_agb_inc", "_AGB", ""), ".csv"))
best_results_combos <- readRDS(paste0("processed_data/best_results_combos", ifelse(what %in% "log_agb_inc", "_AGB", ""), ".rds"))

## crown position data ####
crown <- read.csv("https://raw.githubusercontent.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/master/tree_dimensions/tree_crowns/cored_dendroband_crown_position_data/dendro_cored_full.csv")


# prepare data ####
## keep only 1960 and later
Biol <- Biol[Biol$Year >= 1960, ]

## add crown position ####
if(!all(Biol$tag %in% crown$tag)) stop("Not all tags are in crown data set")
Biol$crown_position <- crown$crown.position[match(Biol$tag, crown$tag)]


## remove tags for which we do not have crown position ####
Biol <- Biol[!is.na(Biol$crown_position), ]

## combine crown position into 2 groups (DC and IS) ####
levels(Biol$crown_position) <- c("CD", "CD", "IS", "IS")

## edit species to be a species-canopy class combination ####
Biol$sp <- paste(Biol$sp, Biol$crown_position, sep = "_")

## remove  species-canopy combinations were we have less than 5 cores ###
sp_can_combi_to_remove <- table(Biol$sp[!duplicated(Biol$coreID)]) 
sp_can_combi_to_remove <- names(sp_can_combi_to_remove[sp_can_combi_to_remove<5])

Biol <- Biol[!Biol$sp %in% sp_can_combi_to_remove, ]


## run analysis ####

# After discussion, we like to get rid of collinearity by picking what variables made more sense biologically to us so here is the set we move on with:
variables_to_keep <- c("pre", "wet", "cld", "tmx", "tmn")
vif(Biol[, variables_to_keep])
vifstep(Biol[, variables_to_keep], th = 3) #--> pre and wet are correlated...
  
# now do a species by species gam ####
  
  
  # create gam formula
  
  full_model_formula <- switch(what, "log_core_measurement" =  paste("log_core_measurement ~ s(dbh, k = 3) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")),
                               log_agb_inc = paste("log_agb_inc ~ s(dbh, k = 3) + s(Year, bs ='re', by = tag) +", paste0("ns(", variables_to_keep, ", 2)", collapse = " + ")))

  ## identify what variables we should keep for each species, looking at the sum of AIC weights ####
  start_time <- Sys.time()
  
  sum_of_weights_for_each_term_by_sp <- NULL
  for(sp in unique(Biol$sp)) {
    print(sp)
    x <- Biol[Biol$sp %in% sp,]
    x$tag <- factor(x$tag)
    x$log_core_measurement <- log(x$core_measurement+0.1)
    x$log_agb_inc <-  log(x$agb_inc + 0.1)
    x <- x[, c("dbh", "Year", "tag", what, variables_to_keep)]
    
    x <- x[!is.na(x[, what]), ]
    fm1 <- gam(eval(parse(text = full_model_formula)), data = x, na.action = "na.fail") # using mixed model makes all variables important which I find suspicous. I feel that s(Year, by = tag) is enough to add some randomness by tag... + plus I am not even sure that actuallydoes what we want.
    dd <- dredge(fm1)
    dd$cw <- cumsum(dd$weight)
    
    sum_of_weights_for_each_term <- dd[, grepl(paste(c(variables_to_keep, "dbh", "Year"), collapse = "|"), names(dd))]
    sum_of_weights_for_each_term <- apply(sum_of_weights_for_each_term, 2, function(x) sum(dd$weight[!is.na(x)]))
    sum_of_weights_for_each_term
    sum_of_weights_for_each_term_by_sp <- rbind(sum_of_weights_for_each_term_by_sp,
                                                c(sum_of_weights_for_each_term))
    
    
    
    # get the results of the model that includes the variables that have sum of weight > 0.9
    
    if(length(names(sum_of_weights_for_each_term)[sum_of_weights_for_each_term > 0.9]) > 0) {
      best_model <- gam(eval(parse(text = paste( what,  "~", paste(names(sum_of_weights_for_each_term)[sum_of_weights_for_each_term > 0.9], collapse = " + ")))), data = x, na.action = "na.fail")
    } else {
      # if none of the variables make the cuttof, still keep dbh and year effect, but stop if they are not in the top model
      dbh_in_top_model <- sum(grepl("dbh", names(dd) [which(dd[1,]=="+")])) ==1
        
      if(!dbh_in_top_model) stop("no variables make the cut and dbh is not in best model...")

      best_model <- gam(eval(parse(text = paste( what,  "~", paste(grep("dbh|Year", names(dd) [which(dd[1,]=="+")], value = T), collapse = " + ")))), data = x, na.action = "na.fail")
      
    }
  
    
    
    # save results for individual species
    assign(paste0(sp, "_dd"), dd)
    assign(paste0(sp, "_best_model"), best_model)
    
    # remove x
    rm(x)
  }
  end_time <- Sys.time()
  
  
  (ellapsed_time <- difftime(end_time, start_time))
  
  rownames(sum_of_weights_for_each_term_by_sp) <- unique(Biol$sp)
  
  sum_of_weights_for_each_term_by_sp
  library(lattice)
  
  levelplot(t(sum_of_weights_for_each_term_by_sp), 
            scales=list(x=list(rot=45)), 
            xlab = "parameter", 
            ylab = "species",
            legend = list(top = list(fun = grid::textGrob("Sum of Weights", y=0, x=1.09))))
  
  # save the plot
  dev.print(png, paste0('results/Species_by_species_GAMS_on_raw_data_by_canopy/Sum_of_AICweights_', what, '.png'),
            width = 6,
            height =8,
            units = "in",
            res = 300)
  
  
  ## second, plot response curves for each species and variables ####
  
  for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
    
    print(sp)
    
    X <- Biol[Biol$sp %in% sp, ]
    X <- X[!is.na(X[, switch(what, log_core_measurement = "core_measurement", log_agb_inc = "agb_inc")]), ]

    best_model <- get(paste0(sp, "_best_model"))
    
    variables_to_look_at <- names((best_model$var.summary))[!names((best_model$var.summary)) %in% c("Year", "tag")]
    
    n_row <- ifelse(length(variables_to_look_at) <= 2, 1, 2)
    n_col <- ifelse(length(variables_to_look_at) ==2, 2, length(variables_to_look_at) %/% 2 +  length(variables_to_look_at)%% 2)
    
    if(length(variables_to_look_at) > 0) {
      par(mfrow=c(ifelse(n_row == 0, 1, n_row), n_col))
      
      
      for(v in variables_to_look_at) {
        
        X$varying_v <- X[, v]
        
        varying_x <- data.frame(floor(min(X$varying_v)): ceiling(max(X$varying_v))) 
        colnames(varying_x) <- v
        
        X$Y <- X[, switch(what, log_core_measurement = "core_measurement", log_agb_inc = "agb_inc")]
        plot(Y+0.1 ~ varying_v, data = X, log = "y", 
             pch = 16,
             # bg = rgb(0,0,0,0.2),
             col = rainbow(length(unique(X$tag)), s = 0.8, alpha = 0.2)[c(1:length(unique(X$tag)))[match(X$tag, unique(X$tag))]],
             main = paste0(sp[1], " - ", v, ifelse(v %in% best_results_combos$climate, paste0("\nfrom ",
                                                                                              paste(month.abb[reference_date[2] - as.numeric(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")])], collapse = " to ")), "")),
             xlab = v,
             ylab = switch(what, log_core_measurement = "core measurement (mm)", log_agb_inc = "AGB increment (Mg C)"),
             border = "grey")
        
        if(length(variables_to_look_at) > 1) {
          constant_variables <- variables_to_look_at[!variables_to_look_at %in% v]
          

            newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(X$", constant_variables, ")", collapse = ", "), ",  Year = median(X$Year), tag = factor(X$tag[1]))"))), varying_x)

        } else {
          newd <- data.frame(varying_x,  Year = median(X$Year), tag = factor(X$tag[1]))
        }
        
        ## prediction
        pt <- predict.gam(best_model, newd, type = "response", exclude =grep("Year", sapply(best_model$smooth, "[[", "label"), value = T), se.fit = T)
        
        ## add preditive line
        lines(exp(pt$fit) ~ varying_x[[1]], lwd = 2)
        lines(exp(pt$fit - 1.96 * pt$se.fit) ~ varying_x[[1]], lwd = 1, lty = 2)
        lines(exp(pt$fit + 1.96 * pt$se.fit) ~ varying_x[[1]], lwd = 1, lty = 2)
        
      }
      
    }
    
    
    
    # save plot
    dev.print(png, paste0('results/Species_by_species_GAMS_on_raw_data_by_canopy/GAM_results_raw_', sp, "_", what, ".png"),
              width = 8,
              height =8,
              units = "in",
              res = 300)
    
  }
  
  ## third,  plot response curves for each variable, with one curve per species ####
  #Create a custom color scale
  
  for(v in c("dbh", variables_to_keep)) {
    print(v)
    
    ## predictions
    pt <- NULL
    for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
      best_model <- get(paste0(sp, "_best_model"))
      
      varying_x <- data.frame(varying_x = seq(min(Biol[Biol$sp %in% sp, v], na.rm = T), max(Biol[Biol$sp %in% sp, v], na.rm = T), length.out = 100))
      
      colnames(varying_x) <- v
      
      constant_variables <- c("dbh", variables_to_keep)[!c("dbh", variables_to_keep) %in% v]
      
   
        newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(Biol$", constant_variables, ", na.rm = T)", collapse = ", "), ",  Year = median(Biol$Year, na.rm = T), tag = factor(Biol[Biol$sp %in% sp,]$tag[1]))"))), varying_x)
        
      if(v %in% names(best_model$var.summary)) {
        pt <- rbind(pt, data.frame(newd, variable = v, species = sp, varying_x = newd[, v], predict.gam(best_model, newd, type = "response", exclude =grep("Year", sapply(best_model$smooth, "[[", "label"), value = T), se.fit = T)))
      }
      
    } # ignore errors
    pt$species <- factor(pt$species, levels = rownames(sum_of_weights_for_each_term_by_sp))
    pt$expfit <- exp(pt$fit)
    pt$lwr <- exp(pt$fit - 1.96 * pt$se.fit)
    pt$upr <- exp(pt$fit + 1.96 * pt$se.fit)
    
    pt$crown_position <- factor(sapply(strsplit(as.character(pt$species), split = "_"), "[[", 2), levels = c("D", "C", "I", "S"))
    
    p <- ggplot(data = pt, aes(x = varying_x, y = expfit, group = species))
    if(!v %in% c("dbh")) p <- p + geom_rect(xmin = mean(Biol[, v], na.rm = T) - sd(Biol[, v], na.rm = T), ymin = min(pt$lwr, na.rm = T), xmax = mean(Biol[, v], na.rm = T) + sd(Biol[, v], na.rm = T), ymax = max(pt$upr), fill = "grey" , alpha=0.01) + geom_vline(xintercept = mean(Biol[, v]), col = "grey")
    
    p <- p + geom_line(aes(col = crown_position)) +
      # scale_x_continuous(trans= ifelse(v %in% "dbh", 'log','identity')) +
      labs(title = paste0(v, ifelse(v %in% best_results_combos$climate, paste0("\nfrom ",
                                                                               paste(month.abb[reference_date[2] - as.numeric(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")])], collapse = " to ")), "")),
           x = v,
           y = "") + #"core measurements") +
      geom_ribbon(aes(ymin=lwr, ymax=upr, col = NULL, bg = crown_position), alpha=0.25) + 
      scale_colour_hue(drop = F) + scale_fill_hue(drop = F) + 
      theme_classic()
    
    assign(paste0("p_", v), p +
             theme(legend.position="none"))
  }
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  
  
  grid.arrange(do.call(arrangeGrob, c(lapply(paste0("p_",  c("dbh", variables_to_keep)), function(x)  get(x)), ncol = 3)),
               g_legend(p),
               nrow = 1,
               widths = c(10, 1))
  
  grid.text(switch (what, log_core_measurement = "core measurement (mm)",
                    log_agb_inc = "AGB increment (Mg C)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90)
  
  
  
  # save plot
  dev.print(png, paste0('results/Species_by_species_GAMS_on_raw_data_by_canopy/ALL_variables_', what, '.png'),
            width = 8,
            height =8,
            units = "in",
            res = 300)

