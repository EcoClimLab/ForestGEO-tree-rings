# Analysis methods

## Climate
We are currently considering two completementary ways of looking at climate sensitivity:

1. The [climwin](https://cran.r-project.org/web/packages/climwin/index.html) R package, described in [this paper](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12590), is currently our best known option for characterizing climate effects. This works with species-level chronologies or individual-level analysis, and can integrate with most types of regression models that can be implemented in R. We have it working with *lmer* (see [issue #1](https://github.com/EcoClimLab/ForestGEO-climate-sensitivity/issues/1)).

2. The "quilt plots" in [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906) work well for visualizing the big-picture of multiple variables/ months / species, and is useful for complementing climwin. So far, we have only done this for species-level residual chronologies, not with individual-level analyses. However, there is no reason (other than potentially complex coding) why we couldn't make similar "quilt plots" for individual-level analyses.


## Potential methods for individual tree growth analyses





