rm(list = ls())
library(mgcv)
library(ggplot2)
library(gridExtra)

sites = c("BCI", "HKK", "SCBI", "LillyDickey", "HarvardForest", "Zofin", 
          "Niobara", "NewMexico", "CedarBreaks", "ScottyCreek")


what_to_show <- c("log_core_measurement", "log_core_measurement_dbh", "log_BAI_dbh", "log_agb_inc_dbh"
)



for(site in sites[]){
  print(site)
  
  for(what in what_to_show) {
    print(what)
    
    
    for(with_Year_or_CO2 in c("", "Year")) {
      print(with_Year_or_CO2)
    
    obj_to_keep <- c(ls(), "obj_to_keep")
    

    
    if(!file.exists( paste0("results/",  ifelse(with_Year_or_CO2 %in% "", "", paste0("with_", with_Year_or_CO2, "/")), what, "/", site, "_all_env.RData"))) next
    
    load( paste0("results/",  ifelse(with_Year_or_CO2 %in% "", "", paste0("with_", with_Year_or_CO2, "/")), what, "/", site, "_all_env.RData")) 

    
    
    ## Plot response curves for each variable, with one curve per species ####
    # first remove any object starting by p_
    rm(list = ls()[grepl("^p_", ls())])
    #Create a custom color scale
    ylim_p <- list()
    ylim_p_int <- list()
    for(v in c( variables_to_keep)) { 
      print(v)
      
      ## predictions
      pt <- NULL
      pt_int <- NULL
      for(sp in rownames(sum_of_weights_for_each_term_by_sp)) {
        best_model <- get(paste0(sp, "_best_model"))
        
        g <- names(which(sapply(clim_var_group, function(x) v %in% x)))
        
        if(exists(paste0(sp, "_best_model_int_", g))) {
          dbh_clim_int = T
          best_model_int <- get(paste0(sp, "_best_model_int_", g))
        } else {
          dbh_clim_int = F
        }
        
        
        x <- Biol[Biol$species_code %in% sp,]
        
        # keep only years with enough data if we look at CO2 or Year
        # if(!with_Year_or_CO2 %in% "") x <- x[x$Year %in% Species_Year_to_keep$Year[Species_Year_to_keep$species_code %in% sp],]
        
        varying_x <- data.frame(varying_x = seq(min(x[, v], na.rm = T), max(x[, v], na.rm = T), length.out = 100)) ; colnames(varying_x) <- v
        # constant_variables <- c("dbh", variables_to_keep)[!c("dbh", variables_to_keep) %in% v]
        constant_variables <- c(variables_to_keep)[!c(variables_to_keep) %in% v]
        
        newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(x$", constant_variables, ", na.rm = T)", collapse = ", "), ")"))), varying_x)  # newd <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(Biol$", constant_variables, ", na.rm = T)", collapse = ", "), ",  Year = median(Biol$Year, na.rm = T), coreID = factor(Biol[Biol$species_code %in% sp,]$coreID[1]))"))), varying_x)
        
        
        
        # if(v %in% names(best_model$var.summary)) {
        # pt <- rbind(pt, data.frame(newd, variable = v, species = sp, varying_x = newd[, v], do.call(rbind, lapply(1:nrow(newd), function(i) data.frame(predict(best_model, newd[i,], type = "link", level = 0, se.fit = T)))), draw = any(grepl(v,names(best_model$coefficients$fixed)))))
        
        
        pt <- rbind(pt, 
                    data.frame(newd, 
                               variable = v, 
                               species = sp, 
                               varying_x = newd[, v], 
                               data.frame(predict(best_model, newd, type = "link", level = 0, se.fit = T)), 
                               sigma = sigma(best_model), # get residual Standard error to correct prediction of mean when back-transforming
                               draw = any(grepl(v,names(best_model$coefficients$fixed))),
                               sig = switch(as.character(sum(summary(best_model)$tTable[grep(v, rownames(summary(best_model)$tTable)), "p-value"] < .05)), "2" = "solid", "1" = "twodash", "0" = "dotted")))
        
        
        if(dbh_clim_int) {
          constant_variables <- constant_variables[!constant_variables%in%"dbh"]
          newd_int <- cbind(eval(parse(text = paste0("data.frame(", paste0(constant_variables, " = median(x$", constant_variables, ", na.rm = T)", collapse = ", "), ", dbh = range(x$dbh))"))), varying_x[rep(1:nrow(varying_x), each = 2),, drop = F]) 
          
          
          pt_int <- rbind(pt_int, 
                          data.frame(newd_int, 
                                     variable = v, 
                                     species = sp, 
                                     varying_x = newd_int[, v], 
                                     data.frame(predict(best_model_int, newd_int, type = "link", level = 0, se.fit = T)), 
                                     sigma = sigma(best_model_int), # get residual Standard error to correct prediction of mean when back-transforming
                                     draw = any(grepl(v, names(best_model_int$coefficients$fixed))),
                                     # sig = switch(as.character(sum(summary(best_model_int)$tTable[grep(v, rownames(summary(best_model_int)$tTable)), "p-value"] < .05)), "2" = "solid", "1" = "twodash", "0" = "dotted"),
                                     DBH = factor(c("min", "max")[as.numeric(as.factor(newd_int$dbh))])
                          ))
        }
        # }
        
      } # ignore warnings
      
      # if(!is.null(pt)) {
      pt$species_code <- pt$species
      pt$species <- factor(pt$species, levels = rownames(sum_of_weights_for_each_term_by_sp))
      pt$species <- factor(paste0(pt$species, " (",tapply(Biol$coreID,  Biol$species_code, function(x) length(unique(x)))[as.character(pt$species)], ")"))
      pt$expfit <- exp(pt$fit + .5*pt$sigma^2) # exp(pt$fit) would give epected median.
      pt$lwr <- exp(pt$fit + .5*pt$sigma^2 - 1.96 * pt$se.fit)# exp(pt$fit - 1.96 * pt$se.fit)
      pt$upr <- exp(pt$fit + .5*pt$sigma^2 + 1.96 * pt$se.fit)
      
      species_colors <- all_species_colors[[site]]
      names(species_colors) <- pt$species[match(names(species_colors), pt$species_code)]
      species_colors <- species_colors[!is.na(names(species_colors))]
      
      
      time_window <- reference_date[2] - as.numeric(best_results_combos[best_results_combos$climate %in% v, c("WindowOpen", "WindowClose")])
      time_window_prev <- time_window <= 0
      time_window <- ifelse(time_window_prev, rev(1:12)[abs(time_window) +1], time_window)
      
      time_window_text <- paste(paste0(ifelse(time_window_prev, "p.", "c."),month.abb[time_window]), collapse = "-") #  paste0("\nfrom ", paste(paste0(ifelse(time_window_prev, "prev. ", "curr. "), month.abb[time_window]), collapse = "\nto "))
      
      
      p <- ggplot(data = pt[pt$draw,], aes(x = varying_x, y = expfit))
      
      # add gray rectangle when relevant
      if(!v %in% c("dbh", "Year")) p <- p + geom_rect(xmin = mean(Biol[, v], na.rm = T) - sd(Biol[, v], na.rm = T), ymin = min(pt$lwr), xmax = mean(Biol[, v], na.rm = T) + sd(Biol[, v], na.rm = T), ymax = max(pt$upr), fill = "grey" , alpha=0.01) + geom_vline(xintercept = mean(Biol[, v], na.rm = T), col = "grey")
      
      
      p <- p + geom_line(aes(col = species), linetype = as.character(p$data$sig)) +
        # scale_x_continuous(trans= ifelse(v %in% "dbh", 'log','identity')) +
        labs(#title = paste0(v, ifelse(v %in% best_results_combos$climate, time_window_text, "")),
          x = paste(v, ifelse(v %in% best_results_combos$climate, time_window_text, ""), variables_units[[v]]),
          y = "") + #"core measurements") +
        geom_ribbon(aes(ymin=lwr, ymax=upr, bg = species), alpha=0.25) +
        scale_color_manual(values = species_colors[unique(pt$species[pt$draw])]) + scale_fill_manual(values = species_colors[unique(pt$species[pt$draw])]) +
        # scale_colour_hue(drop = F) + scale_fill_hue(drop = F) + 
        theme_classic() +
        theme(text = element_text(size = 10))
      
      if(any(pt$draw) | v %in% "Year") {
        assign(paste0("p_", v), p +
                 theme(legend.position="none")) 
      }
      
      # save the ylimits ####
      ylim_p[[v]] <- range(pt$lwr, pt$upr)
      # if(!with_CO2) ylim_p[[what]][[v]][[site]] <- range(pt$lwr, pt$upr)
      
      
      if(!is.null(pt_int)) {
        pt_int$species_code <- pt_int$species
        pt_int$species <- factor(pt_int$species, levels = rownames(sum_of_weights_for_each_term_by_sp))
        pt_int$species <- factor(paste0(pt_int$species, " (",tapply(Biol$coreID,  Biol$species_code, function(x) length(unique(x)))[as.character(pt_int$species)], ")"))
        pt_int$expfit <- exp(pt_int$fit + .5*pt_int$sigma^2) # exp(pt_int$fit) would give epected median.
        pt_int$lwr <- exp(pt_int$fit + .5*pt_int$sigma^2 - 1.96 * pt_int$se.fit)# exp(pt_int$fit - 1.96 * pt_int$se.fit)
        pt_int$upr <- exp(pt_int$fit + .5*pt_int$sigma^2 + 1.96 * pt_int$se.fit)
        
        
        p_int <- ggplot(data = pt_int[pt_int$draw,], aes(x = varying_x, y = expfit, group = dbh))
        
        # add gray rectangle when relevant
        if(!v %in% c("dbh", "Year")) p_int <- p_int + geom_rect(xmin = mean(Biol[, v], na.rm = T) - sd(Biol[, v], na.rm = T), ymin = min(pt$lwr), xmax = mean(Biol[, v], na.rm = T) + sd(Biol[, v], na.rm = T), ymax = max(pt$upr), fill = "grey" , alpha=0.01) + geom_vline(xintercept = mean(Biol[, v], na.rm = T), col = "grey")
        
        
        p_int <- p_int + geom_line(aes(col = species, linetype = DBH)) +
          # scale_x_continuous(trans= ifelse(v %in% "dbh", 'log','identity')) +
          labs(#title = paste0(v, ifelse(v %in% best_results_combos$climate, time_window_text, "")),
            x = paste(v, ifelse(v %in% best_results_combos$climate, time_window_text, ""), variables_units[[v]]),
            y = "") + #"core measurements") +
          geom_ribbon(aes(ymin=lwr, ymax=upr, bg = species), alpha=0.25) +
          scale_color_manual(values = species_colors[unique(pt$species[pt$draw])]) + scale_fill_manual(values = species_colors[unique(pt$species[pt$draw])]) +
          # scale_colour_hue(drop = F) + scale_fill_hue(drop = F) + 
          theme_classic() +
          theme(text = element_text(size = 10))
        
        if(any(pt_int$draw) | v %in% "Year") {
          assign(paste0("p_int_", v), p_int) 
        }
        
        # save the ylimits ####
        ylim_p_int[[v]] <- range(pt_int$lwr, pt_int$upr)
        # if(!with_CO2) ylim_p[[what]][[v]][[site]] <- range(pt$lwr, pt$upr)
      }
      #  p$data <- pt_int
      #  p <- p + aes(group = factor(p$data$dbh))
      # p$layers[[3]]$aes_params$linetype <- as.character(p$data$sig)
      # p$layers[[3]]$aes_params$size <- p$data$linesize*0.5
      
      
    }
    
    g_legend <- function()
    {
      a.gplot <- ggplot(data = pt, aes(x = varying_x, y = expfit)) + geom_line(aes(group = species, col = species)) +  geom_ribbon(aes(ymin=lwr, ymax=upr, bg = species), alpha=0.25) +
        scale_color_manual(values = species_colors) + scale_fill_manual(values = species_colors)
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    
    # existing_plots <- paste0("p_",  c("dbh", variables_to_keep))
    existing_plots <- ls()[grepl("^p_[^in]", ls())]
    existing_plots_int <- ls()[grepl("^p_int_", ls())]
    
    grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x) {
      x <- get(x)
      x$layers[[1]]$aes_params$ymin <- range(ylim_p)[1]
      x$layers[[1]]$aes_params$ymax <- range(ylim_p)[2]
      x + ylim(range(ylim_p))
    }), ncol = ifelse(length(existing_plots) %in% 4, 2, length(existing_plots)))),
    g_legend(),
    nrow = 1,
    widths = c(10, 1))
    
    grid::grid.text(switch (gsub("_dbh", "", what), log_core_measurement = "core measurement (mm)",
                            log_agb_inc = "AGB increment (Mg C)", log_BAI = "BAI (cm2)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90, gp = grid::gpar(fontsize = 10))
    
    
    
    # save png plot ####
    png(paste0('results/', ifelse(with_Year_or_CO2 %in% "", "", paste0("with_", with_Year_or_CO2, "/")), ifelse(solution_to_global_trend %in% "none", "", paste0(solution_to_global_trend, "/")), what, "/", site, "/GLS_ALL_variables_responses_", site, '.png'),
        width = 10,
        height =8,
        units = "in",
        res = 300)
    
    grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots, function(x) {
      x <- get(x)
      x$layers[[1]]$aes_params$ymin <- range(ylim_p)[1]
      x$layers[[1]]$aes_params$ymax <- range(ylim_p)[2]
      x + ylim(range(ylim_p))
    }), ncol = ifelse(length(existing_plots) %in% 4, 2, length(existing_plots)))),
    g_legend(),
    nrow = 1,
    widths = c(10, 1))
    
    grid::grid.text(switch (gsub("_dbh", "", what), log_core_measurement = "core measurement (mm)",
                            log_agb_inc = "AGB increment (Mg C)",
                            log_BAI = "BAI (cm2)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90, gp = grid::gpar(fontsize = 10))
    dev.off()
    
    # save interaction png plot ####
    if(length(existing_plots_int)>0) {
      png(paste0('results/', what, "/", site, "/GLS_ALL_variables_responses_", site, '_dbh_interaction.png'),
          width = 10,
          height =8,
          units = "in",
          res = 300)
      
      grid.arrange(do.call(arrangeGrob, c(lapply(existing_plots_int, function(x) {
        x <- get(x)
        x$layers[[1]]$aes_params$ymin <- range(ylim_p_int)[1]
        x$layers[[1]]$aes_params$ymax <- range(ylim_p_int)[2]
        x + ylim(range(ylim_p_int))
      }), ncol = ifelse(length(existing_plots_int) %in% 4, 2, length(existing_plots_int)))),
      nrow = 1,
      widths = c(10, 1))
      
      grid::grid.text(switch (gsub("_dbh", "", what), log_core_measurement = "core measurement (mm)",
                              log_agb_inc = "AGB increment (Mg C)",
                              log_BAI = "BAI (cm2)"), x = unit(0.01, "npc"), y = unit(.51, "npc"), rot = 90, gp = grid::gpar(fontsize = 10))
      dev.off()
      
    }
    # save plots at this point to later fetch them  ####
    save(list = grep("^p_|pt|clim_var_group$|ylim_p|species_colors", ls(), value = T), file = paste0('results/', ifelse(with_Year_or_CO2 %in% "", "", paste0("with_", with_Year_or_CO2, "/")), ifelse(solution_to_global_trend %in% "none", "", paste0(solution_to_global_trend, "/")), what, "/", site, "/env.RData"))
    
    rm(list = ls()[!ls() %in% obj_to_keep])
    
    }  
  }
}
