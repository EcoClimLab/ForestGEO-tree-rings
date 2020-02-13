# Analysis methods

## Overview

There are 3 basic steps:
1. Identify most important climatic drivers
2. Identify best multivariate statistical model accounting for climate and DBH
3. Scale to the ecosystem by combining statistical model with census data

## Step 1. Identify most important climatic drivers

[Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906), [GitHub repo with full results](https://github.com/SCBI-ForestGEO/climate_sensitivity_cores) analyzed climate sensitivity using traditional dendro methods and gave a reliable (as far as we know!) picture of climate sensitivity at SCBI. From our explorations of other analysis methods on the same data, it seems that pulling out a reliable climate signal *requires* indivdidual-level detrending of cores. Our attempt to do this with raw growth data, using climwin and a mixed effects model to account for DBH and year, resulted in weak or absent climate correlations that didn't align with Helcoski results. However, if we fit splines to individual core sequences and then analyze these data in *climwin*, we get results that are generally very consistent with Helcoski (except that they don't pull out an opposite signal from the previous growing season).

We have two completementary ways of looking at climate sensitivity:

1. The [climwin](https://cran.r-project.org/web/packages/climwin/index.html) R package, described in [this paper](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12590), is currently our best known option for identifying the most important climate drivers. This works with species-level chronologies or individual-level analysis, and can integrate with most types of regression models that can be implemented in R. We have it working with *lmer* (see [issue #1](https://github.com/EcoClimLab/ForestGEO-climate-sensitivity/issues/1)).

2. The "quilt plots" in [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906) work well for visualizing the big-picture of multiple variables/ months / species, and are a useful complement to *climwin*. So far, we have only done this for species-level residual chronologies, not with individual-level analyses. However, there is no reason (other than potentially complex coding) why we couldn't make similar "quilt plots" for individual-level analyses.

### Considerations
#### Species- versus individual-level analysis

Species-level chronologies provide the clearest picture of climate sensitivity, using methods optimized by dencrochronologists (e.g., Helcoski et al). We may consider whether it would be better to use this approach.

Individual-level residual analysis seems to work fine for pulling out the climate signal, and is our current approach.


## Step 2. Identify best multivariate statistical model accounting for climate and DBH

### Considerations:
- DBH must be accounted for in the model. [Peters et al. (2015)](https://onlinelibrary.wiley.com/doi/epdf/10.1111/gcb.12826) reviews some methods

- should we account for age too?

### Current method:

- Using *climwin* with with mixed effects model: `diameter.growth ~ [climate] + DBH + year + treeID (random effect)`

- Results are [here](https://github.com/EcoClimLab/ForestGEO-climate-sensitivity/tree/master/results)

- Some of the results align with [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906), but many don't (e.g., positive effect of DTR, where Helcoski found negative).

- Some potential clean-up that may make results more clear: remove outliers, ln-transform age and growth?, remove DTR as predictor variable.

### What others have done:
- [Foster et al. (2016)](https://github.com/EcoClimLab/ForestGEO-climate-sensitivity/blob/master/methods/references/Foster_et_al-2016-Global_Change_Biology.pdf) used an individuallevel Baysian heriarchical model to characterize to what extent whole-tree growth depends on climatic, tree, and stand structural variables.

- M. Ross Alexander and colleagues (in review-TREAT CONFIDENTIALLY) "used generalized additive models (GAMs) using the gam function in the mgcv R package to allow for non-linear climate and size effects without prescribing a priori functional forms (Polansky & Robbins 2013; Peters et al.2015, Simpson 2018, Wood 2003).  ....Year effects were estimated using a 4-knot thin plate regression spline fit by sampling plot to capture stand dynamics such as suppression or release events that may have influenced growth through time. Size effects were estimated for each species using 3-knot thin-plate regression splines fit to reconstructed DBH (Wood 2003)."

    Peters, R.L., Groenendijk, P., Vlam, M. & Zuidema, P.A. (2015). [Detecting long‐term growth trends using tree rings: a critical evaluation of methods](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.12826). Glob Change Biol, 21, 2040–2054.

    Polansky, L. & Robbins, M.M. (2013). Generalized additive mixed models for disentangling long‐term trends, local anomalies, and seasonality in fruit tree phenology. Ecology and Evolution, 3, 3141–3151.

    Simpson, G.L. (2018). [Modelling palaeoecological time series using generalised additive models](https://www.frontiersin.org/articles/10.3389/fevo.2018.00149/full). Frontiers in Ecology and the Environment, 6, 1–21.

    Wood, S.N. (2003). Thin plate regression splines. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 65, 95–114.


## Step 3. Scale to the ecosystem by combining statistical model with census data
### What others have done:

- [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906) used tree-rings to charactize growth as a function of DBH, and assumed that the proportional response to climate was consistent across size classes. We want to move beyond that.
- [Teets et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.14120) used stratefied random sample based on community composition. They summed across individuals to get values proportional to ecosystem-level productivity. This approach does not work when sampling is not carefully stratefied relative to community composition.



## Additional Methodological Considerations



### Response metric


- **ring width (radial growth)** - Annual ring width from individual-level chronologies. Current intention is to focus on this. Used in [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906).


#### What others have done:

- **ring width index from species-level chronologies** - Traditional dendro method, used in [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906). Advantage is that it correlates more strongly to climate. Disadvantages are that it greatly reduces sample size and that its very limited when it comes to predicting growth for any given tree.

- **basal area increment (BAI)** - BAI is commonly considered (by dendrochronologists) to be independent of DBH (or at least moresy than ring width) after trees pass a juvenile growth phase, the logic being that if a tree produces a constant amount of stem wood area each year, BAI will hold steady and ring width will decline. However, while this may hold true in some cases such as very old conifers, in most closed-canopy forests around the world, radial/ diameter growth *increases* with DBH (e.g., Muller-Landau et al. 2006, [Anderson-Teixeira et al. 2015](https://www.researchgate.net/publication/277918165_Size-related_scaling_of_tree_form_and_function_in_a_mixed-age_forest)). There are some likely differences between trees represented by tree-rings (survivors) and the community as a whole, but in [Helcoski et al. 2019](https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.15906) we show a positive relationship of radial growth to DBH from tree-rings. The implication--based on basic geometry--is that BAI increases even more strongly with DBH than does radial growth, and hence BAI should *not* be considered independent of DBH. Furthermore, BAI in itself is not very meaningful, as wood is added in three dimensions, and BAI only accounts for two. Thus, KAT sees no *a priori* motivation to use BAI over ring width. BAI was used by and Alexander et al. (in review) but without assumption that BAI is independent of DBH. 

- **biomass increments** - If we go beyond ring width, biomass increment would be preferable to BAI. I currently think it makes more sense to scale to biomass after statistically predicting ring width increment, but biomass increment could work. **It's worth giving further consideration to the issue of whether to scale to biomass before or after running the statistica model.** 

    [Foster et al. (2016)](https://github.com/EcoClimLab/ForestGEO-climate-sensitivity/blob/master/methods/references/Foster_et_al-2016-Global_Change_Biology.pdf) used indivdiual-level biomass increments. 
   
   [Teets et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.14120) summed biomass increments from cored trees at the species and stand level. They state, "Summing annual biomass increments maintains the signals from individual trees while homogenizing the effect of tree-level factors such as tree size and crowding." This logic makes sense. However, it requires that sampling be proportional to community composition. 
   
- **percent biomass change** - potential option that would be less strongly dependent on size than biomass increments



