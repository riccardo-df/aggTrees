# Aggregation Trees
R package to implement aggregation trees, a data-driven approach to discover heterogeneous subpopulations in a selection-on-observables framework that avoids the risk of data snooping and the drawbacks of pre-analysis plans.

Aggregation trees provide a completely nonparametric approach to constructing partitions of the population that can handle covariate spaces of arbitrary dimensions and an arbitrary number of interactions among covariates. The approach consists of three steps:

1. Estimate the conditional average treatment effects (CATES);
2. Construct a decision tree using the CATEs;
3. Generate a sequence of "optimal” partitions of the covariate space by pruning the tree.

Optimality here refers to the fact that, at each granulairty level, the loss in explained heterogeneity resulting from aggregation is minimized. Notice that the sequence of partitions is nested, as we never undo previous aggregations. This guarantees consistency across the different granularity levels.
  
## Installation  
The current development version the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/aggTrees")
library(aggTrees)
```
