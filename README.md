# ForestGEO Climate Sensitivity

Repository for analysis of climate sensitivity in ForestGEO dendro data (Q1 in 2020 Scholarly Studies grant)

## Hypothesis
*From grant proposal, Hypothesis 1:*

Forest woody productivity is likely to decline under warmer and drier future climate conditions. Efforts to understand and model the inter-annual climate sensitivity of forest ecosystem productivity have focused primarily on gross primary production (GPP), yet our recent work indicates that there are potentially important disconnects between the climate responses of GPP and aboveground net primary productivity (ANPPwoody) that are not captured in ecosystem models. Specifically, whereas GPP may respond positively or negatively to both temperature and precipitation, tree growth and ANPPwoody appear to generally be greatest under relatively cool, wet conditions. Further decoupling is introduced through carbohydrate storage, which creates lags whereby growth is sensitive to climatic conditions of the previous growing season. Given these poorly understood disconnects between GPP and ANPPwoody, empirically characterizing the climate sensitivity of ANPPwoody across the world’s major forest biomes will be crucial to predicting responses of ANPPwoody and biomass to future shifts in climate conditions. As tree-rings commonly reveal negative growth responses to temperature and drought, we hypothesize that ANPPwoody will be reduced at all our sites, which are expected to experience increases in temperature and often precipitation, and increased drought. 

## Data

Data are stored in the [ForestGEO_dendro repository](https://github.com/EcoClimLab/ForestGEO_dendro), which will be private for the foreseeable future. Only data whose PIs have agreed to be made public should be placed in this repository.

SCBI data are public [here](https://github.com/SCBI-ForestGEO/SCBI-ForestGEO-Data/tree/master/tree_cores) and can be copied into this repository.

## Methods
*From grant proposal:*

**Study Sites & Data** This study will include 10 ForestGEO sites representing a range of forest types, climates, and climate trends. Tree cores have already been collected at all sites and, in the vast majority of cases, cross-dated following standard dendrochronological practices. We will complete cross-dating for a few species to maximize our coverage of ANPPwoody. Most sites, including the few where species with chronologies comprise <75% of ANPPwoody (i.e., tropical sites, where few species form annual rings), also have ≥10 years annual tree growth data from dendrometer bands. 

**Analyses** To address our hypothesis, we will use data from tree cores, usually supplemented with dendrometer bands, to characterize the climate sensitivity of ANPPwoody of the entire forest community for all sites. In most cases, the combination of these two data types will allow us to include species representing >90% of ANPPwoody. First, using existing R scripts developed in KAT’s lab, we will analyze tree-ring chronologies to identify the climate variables-months combinations with the strongest influence on the annual growth of each species. We will then combine tree-ring and dendrometer band data in mixed effects models to examine the joint and interactive effects of key climate variables with one another and with tree size, using model averaging (via the `MuMln` package in R) to combine top candidate models into the best predictive model. We will scale predicted species and individual-level responses to the ecosystem level using ForestGEO census data. Finally, we will use the individual-based statistical model (with interactive terms) to predict ANPPwoody under climatic conditions projected through 2100.
