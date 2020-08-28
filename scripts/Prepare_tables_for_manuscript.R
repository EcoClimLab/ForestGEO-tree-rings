# --- create tables for manuscript --- #

# clear environment ####
rm(list = ls())

# load libraries ####

# prepare sites abbrevaiations ####
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
# load data ####
species_list <- read.csv("https://raw.githubusercontent.com/EcoClimLab/ForestGEO_dendro/master/data/species%20list/sitespecies.csv?token=AEWDCINETB6DDQ5LD37JT7C7KJIX2", stringsAsFactors = F)
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

## save

write.csv(species_summary, "doc/manuscript/tables_figures/species.csv", row.names = F)
