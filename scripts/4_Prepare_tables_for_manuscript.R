# --- create tables for manuscript --- #

# clear environment ####
rm(list = ls())

# load libraries ####
library(tidyverse)

# prepare sites abbrevaiations ####
# give site abbrevationtion in paper
sites_abb <- list(BCI  = "BCI",
                  HKK = "HKK",
                  NewMexico = "LT",
                  CedarBreaks = "CB",
                  SCBI = "SCBI",
                  LillyDickey = "LDW",
                  HarvardForest = "HF",
                  Nebraska = "NE",
                  Zofin = "ZOF",
                  ScottyCreek = "SC")
# load data ####
species_list <- read.csv("https://raw.githubusercontent.com/EcoClimLab/ForestGEO_dendro/master/data/species%20list/sitespecies.csv?token=AEWDCIKFARYAZGEDSDHZEOC7VPV64", stringsAsFactors = F)
# bark <- read.csv("https://raw.githubusercontent.com/EcoClimLab/ForestGEO_dendro/master/data/bark/bark_depth.csv?token=AEWDCIMNKUU2YC7ET5X5HZK7KJK3W")

load("results/BCI_all_env.RData")

# prepare species list ####
species_summary <- data.frame(species.code = species_list$species_code,
                              family = species_list$family,
                              latin.name = species_list$latin_name,
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

species_sites$site <- unlist(sites_abb[gsub("\\d", "", row.names(species_sites))])

species_summary$sites.sampled <- tapply(species_sites$site, species_sites$species_code, paste0, collapse = ", ")[as.character(species_summary$species.code)]

## add bark allomatries
species_summary$bark.allometry <- tapply(paste(ifelse(is.na(species_list$bark_species) | species_list$bark_species == "", "neglected", species_list$bark_species), "in", species_list$site), species_list$species_code, paste, collapse = ", ")[as.character(species_summary$species.code)]

## order by species code
species_summary <- species_summary[order(species_summary$species.code), ]

## remove species that are not sampled in any site
species_summary <- species_summary[!is.na(species_summary$sites.sampled), ]


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

# get date range
C <- summary_cores %>%
  group_by(site, species_code) %>%
  summarize(date.range = unique(paste(start_year , end_year , sep = "-") )[1]) # taking first value of unique because this is from the log_core_measuremnt_dbh and that has one mode value as the rest (since others are NA for the first value in the data)


## put everything together
summary_samples <- left_join(left_join(A, B), C)

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

