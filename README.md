# Aggregation Trees
R package to implement aggregation trees, a data-driven methodology to discover heterogeneous subpopulations while avoiding the risk of data snooping and the drawbacks of pre-analysis plans.  
  
Aggregation trees provide a completely nonparametric approach to constructing partitions of the population that can handle covariate spaces of arbitrary dimensions and an arbitrary number of interactions among covariates.  
  
As the data drive the discovery of groups, aggregation trees allow for uncovering unexpected heterogeneity. Moreover, using cost-complexity pruning, a sequence of "optimal" partitions of the population is derived, one for each granularity level. Because the resulting sequence is nested, previous aggregations are never undone when moving to coarser levels. Therefore, consistency is guaranteed across the different granularity levels.

## Installation  
The current development version the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/aggTrees")
library(aggTrees)
```
