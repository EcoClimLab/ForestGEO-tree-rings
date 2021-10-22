# --- create tables for manuscript --- #

# clear environment ####
rm(list = ls())

# load libraries ####
library(tidyverse)

# prepare sites abbrevaiations ####
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
                  # Hansley = "Hansley",
                  Zofin = "ZOF",
                  ScottyCreek = "SC")
# load data ####
species_list <- read.csv("https://raw.githubusercontent.com/EcoClimLab/ForestGEO_dendro/master/data/species%20list/sitespecies.csv?token=AEWDCIPU2N4SAZTFLW5S623BPRH3Y")
bark_allo <- read.csv("https://raw.githubusercontent.com/EcoClimLab/ForestGEO_dendro/master/data_processed/dbh_to_bark_allometries_table.csv?token=AEWDCILSHBGE6WBRU5XUFF3BPRH52")

load("results/log_core_measurement/BCI_all_env.RData")

# prepare species list ####
species_summary <- data.frame(species.code = species_list$species_code,
                              family = species_list$family,
                              latin.name = species_list$latin_name,
                              accepted.name =species_list$accepted_latin_name_with_authority,
                              'sites.sampled' = NA,
                              leaf.type = species_list$leaf_type,
                              leaf.phenology = species_list$leaf_pheno,
                              light.requirements = species_list$light_requirements,
                              bark.allometry = NA
)

## remove duplicates
species_summary <- species_summary[!duplicated(species_summary),]

## add list of sites

species_sites <- data.frame(species_code = unlist(sapply(all_Biol, function(x) as.character(unique(x$species_code)))))
species_sites <- species_sites[-grep("Hansley", row.names(species_sites)),,drop=FALSE] # remove Hansley if is there

species_sites$site <- unlist(sites_abb[gsub("\\d", "", row.names(species_sites))])

species_summary$sites.sampled <- tapply(species_sites$site, species_sites$species_code, paste0, collapse = ", ")[as.character(species_summary$species.code)]

## add bark allomatries
species_summary$bark.allometry <- tapply(ifelse(is.na(species_list$bark_species) | species_list$bark_species == "", "neglected", species_list$bark_latin), species_list$species_code, function(x) paste(unique(x), collapse = ", "))[as.character(species_summary$species.code)] #tapply(paste0(ifelse(is.na(species_list$bark_species) | species_list$bark_species == "", "neglected", paste(species_list$bark_latin, "in ")), species_list$bark_site), species_list$species_code, function(x) paste(unique(x), collapse = ", "))[as.character(species_summary$species.code)]


species_summary$bark.allometry <- gsub(", neglected.*$", "", species_summary$bark.allometry ) # remove repetition of neglected
species_summary$bark.allometry <- gsub(" ,", ",", species_summary$bark.allometry ) # remove space before comma

## change bark allometry to allometry ID
species_summary$bark.allometry <- ifelse(species_summary$bark.allometry %in% "neglected", "neglected", bark_allo$ID[match(paste0("$", gsub(" ", "$ $", species_summary$bark.allometry), "$"), bark_allo$species)])

## order by latin name
species_summary <- species_summary[order(species_summary$accepted.name), ]

## remove species that are not sampled in any site
species_summary <- species_summary[!is.na(species_summary$sites.sampled), ]

# change needleleaf to conifer

species_summary$leaf.type <- gsub("needleleaf", "conifer", species_summary$leaf.type)


## italicize latin names
species_summary$latin.name  <- paste0("$", gsub(" ", "$ $", species_summary$latin.name) ,  "$")

species_summary$accepted.name[!species_summary$accepted.name %in% ""] <- sapply(strsplit(species_summary$accepted.name[!species_summary$accepted.name %in% ""] , " "), function(x) {
 x[c(1,2)] <- paste0("$", x[c(1,2)] ,  "$")
x <-  paste(x , collapse = " ")
return(x)
})

## replace underscores by spaces

names(species_summary) <- gsub("\\.|_", " ", names(species_summary))


## save

write.csv(species_summary, "doc/manuscript/tables_figures/species.csv", row.names = F)

# prepapre summary cores ####
summary_cores <- read.csv("results/summary_cores_analyzed.csv")

names(summary_cores)

## get n trees all and n cores all
idx = !grepl("dbh", summary_cores$what)

A <- summary_cores[idx, ] %>%
  group_by(site, species_code) %>%
  summarize(n.trees.all = unique(n_trees),
            n.cores.all = unique(n_cores))


## get n trees all and n cores all
idx = grepl("dbh", summary_cores$what)

B <- summary_cores[idx, ] %>%
  group_by(site, species_code) %>%
  summarize(n.trees.dbh = unique(n_trees),
            n.cores.dbh = unique(n_cores),
            dbh.range.sampled = unique(paste(min_DBH_cores_sampled, max_DBH_cores_sampled, sep = "-")),
            dbh.range.reconstructed = unique(paste(min_DBH_cores_reconstructed , max_DBH_cores_reconstructed , sep = "-"))[1] # taking first value of unique because this is from the log_core_measuremnt_dbh and that has one mode value as the rest (since others are NA for the first value in the data)
            )

## get date range
C <- summary_cores %>%
  group_by(site, species_code) %>%
  summarize(date.range = unique(paste(start_year , end_year , sep = "-") )[1]) # taking first value of unique because this is from the log_core_measuremnt_dbh and that has one mode value as the rest (since others are NA for the first value in the data)


## get if Year was tested (looking at BAI only)
sp_analyzed_year <- list.files("results/with_Year/log_BAI_dbh/", recursive = T, pattern = "_model_comparisons.csv")

summary_cores$year_analyzed <- paste0(summary_cores$site, "/", summary_cores$species_code, sep = "_model_comparisons.csv") %in% sp_analyzed_year

idx = summary_cores$what %in% "log_BAI_dbh"

D <- summary_cores[idx,] %>%
  group_by(site, species_code) %>%
  summarize(year_analyzed = year_analyzed)

## put everything together
summary_samples <- left_join(left_join(left_join(A, B), C), D)

## replace n_trees_dbh by 0 if NA
summary_samples$n.trees.dbh[is.na(summary_samples$n.trees.dbh)] <- 0
summary_samples$n.cores.dbh[is.na(summary_samples$n.cores.dbh)] <- 0

## replace dbh.range.sampled NA-NA by "unknown"
summary_samples$dbh.range.sampled[summary_samples$dbh.range.sampled %in% "NA-NA"] <- "unknown"
  
## remove Hansley if there (should not starting  next run)
summary_samples <- droplevels(summary_samples[!summary_samples$site %in% "Hansley",])

## change site names to abbreviations
summary_samples$site <- unlist(sites_abb[as.character(summary_samples$site)])

## add asterix for dbh.range.reconstructed
names(summary_samples) <- gsub("dbh.range.reconstructed", "dbh.range.reconstructed\\*", names(summary_samples))

## replace underscores by spaces

names(summary_samples) <- gsub("\\.|_", " ", names(summary_samples))

## save
write.csv(summary_samples, "doc/manuscript/tables_figures/sampling_details.csv", row.names = F, quote = F)

# # sites and species combo ####
# site_sp_combo <- unique(do.call(rbind, all_Biol)[, c("site", "species_code")])
# site_sp_combo$site <- unlist(sites_abb[as.character(site_sp_combo$site)])
# 
# site_sp_combo <- site_sp_combo[order(site_sp_combo$site, site_sp_combo$species_code), ]
# 
# write.csv(site_sp_combo, "doc/manuscript/tables_figures/site_species_combos.csv", row.names = F)

