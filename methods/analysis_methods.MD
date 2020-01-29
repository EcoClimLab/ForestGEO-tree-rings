# Analysis methods

## Climate
We are currently considering two completementary ways of looking at climate sensitivity:

1. The [climwin](https://cran.r-project.org/web/packages/climwin/index.html) R package, described in [this paper](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12590), is currently our best known option for characterizing climate effects. This works with species-level chronologies or individual-level analysis, and can integrate with most types of regression models that can be implemented in R. We have it working with *lmer* (see [issue #1](https://github.com/EcoClimLab/ForestGEO-climate-sensitivity/issues/1)).

2. The "quilt plots" in [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906) work well for visualizing the big-picture of multiple variables/ months / species, and is useful for complementing climwin. So far, we have only done this for species-level residual chronologies, not with individual-level analyses. However, there is no reason (other than potentially complex coding) why we couldn't make similar "quilt plots" for individual-level analyses.

## Species- versus individual-level analysis

Species-level chronologies provide the clearest picture of climate sensitivity, using methods optimized by dencrochronologists. However, this approach is limited in its ability to scale to ecosystem-level productivity, particularly in that it cannot represent interactions with tree size or other relevant factors. 

Here's the original plan from the proposal:

*First, using existing R scripts developed in KAT’s lab ([Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906)), we will analyze tree-ring chronologies to identify the climate variables-months combinations with the strongest influence on the annual growth of each species. We will then combine tree-ring and dendrometer band data in mixed effects models to examine the joint and interactive effects of key climate variables with one another and with tree size, using model averaging (MuMln in R) to combine top candidate models into the best predictive model.*

*Climwin* allows the option of skipping the first step. 

## Potential methods for individual tree growth analyses


## Scaling to the ecosystem
What others have done:
-[Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906) used tree-rings to charactize growth as a function of DBH, and assumed that the proportional response to climate was consistent across size classes. We want to move beyond that.
- [Teets et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.14120) used stratefied random sample based on community composition. They summed across individuals to get values proportional to ecosystem-level productivity. This approach does not work when sampling is not carefully stratefied relative to community composition.

