# Aggregation Trees
R package to implement aggregation trees, a data-driven approach to discover heterogeneous subpopulations in a selection-on-observables framework that avoids the risk of data snooping and the drawbacks of pre-analysis plans.

Aggregation trees provide a completely nonparametric approach to constructing partitions of the population that can handle covariate spaces of arbitrary dimensions and an arbitrary number of interactions among covariates. The approach consists of three steps:

1. Estimate the conditional average treatment effects (CATES);
2. Construct a decision tree using the CATEs;
3. Generate a sequence of “optimal” partitions of the covariate space by pruning the tree.

Optimality here refers to the fact that, at each granulairty level, the loss in explained heterogeneity resulting from aggregation is minimized. Notice that the sequence of partitions is nested, as we never undo previous aggregations. This guarantees consistency across the different granularity levels.

We use different samples to estimate the CATEs (step 1) and construct the sequence of partitions (steps 2 and 3). We call these “estimation sample” and “aggregation sample”. To conduct valid inference, we need to estimate honest trees. For this, we need an additional “honest sample”.

Family of \*\_rpart functions.
  
## Installation  
The current development version the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/aggTrees")
library(aggTrees)
```

## Usage Examples
This section demonstrates how to use aggregation trees to discover heterogeneous subpopulations.

First, we need to estimate the CATEs. This must be done outside the `aggTrees` package. The rationale for this is that one can use any estimator (or an ensemble of more estimators). Here we use the `causal_forest` function from the `grf` package:

```
## Generate data.
set.seed(1986)

n <- 3000
k <- 3

X <- matrix(rnorm(n * k), ncol = k)
D <- rbinom(n, size = 1, prob = 0.5)
mu0 <- 0.5 * X[, 1]
mu1 <- 0.5 * X[, 1] + X[, 2] # This implies that tau(x) = X_2.
y <- mu0 + D * mu1 + rnorm(n)

## Split the sample.

library(grf)

```
