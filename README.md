# Aggregation Trees
R package to implement aggregation trees, a nonparametric data-driven approach to discovering heterogeneous subgroups in a selection-on-observables framework. Additionally, it provides useful functions to work with `rpart` objects.

The approach consists of three steps:

1. Estimate conditional average treatment effects (CATEs);
2. Aggregate CATEs into a decision tree;
3. Prune the tree.

This way, we generate a sequence of groupings, one for each granularity level. The resulting sequence is nested in the sense that subgroups formed at a given
level of granularity are never broken at coarser levels. This guarantees consistency of the results across the different granularity levels, generally considered a basic requirement that every classification system should satisfy. Moreover, each grouping features an optimality property in that it ensures that the loss in
explained heterogeneity resulting from aggregation is minimized.

Given the sequence of groupings, we can estimate group average treatment effects (GATEs) as we like. The package supports two estimator, based on differences in mean outcome between treated and control units (unbiased in randomized experiments) and on sample averages of doubly-robust scores (unbiased also in observational studies). The package also allows o get standard errors for the GATEs by estimating via OLS appropriate linear models. An honesty condition is required to construct valid confidence intervals. Thus, different subsamples must be used to construct the tree and estimate the linear models. 

## Installation  
The current development version of the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/aggTrees") # run install.packages("devtools") if needed.
library(aggTrees)
```

## Usage Examples
This section demonstrates how to use the package. Let us generate some data:

```
## Generate data.
set.seed(1986)

n <- 3000
k <- 3

X <- matrix(rnorm(n * k), ncol = k)
colnames(X) <- paste0("x", seq_len(k))
D <- rbinom(n, size = 1, prob = 0.5)
mu0 <- 0.5 * X[, 1]
mu1 <- 0.5 * X[, 1] + X[, 2] # This implies that tau(x) = X_2.
y <- mu0 + D * mu1 + rnorm(n)
```

As a first step, we need to estimate CATEs. We can do this with any estimator we like. Then, in the second step we construct a tree using the CATEs as an outcome. Given the tree, we can compute node predictions (i.e., GATEs) as we like.

I split the data into a training sample and an honest sample. We use the training sample to estimate the CATEs and construct/prune the tree. Then, we use the honest sample to populate each node. This way, we can later construct valid confidence intervals for GATEs. This is achieved in the following chunk of code.

```
## Sample splitting.
splits <- sample_split(n_train, training_frac = 0.5)
training_idx <- splits$training_idx
honest_idx <- splits$honest_idx

y_tr <- y[training_idx]
D_tr <- D[training_idx]
X_tr <- X[training_idx, ]

y_hon <- y[honest_idx]
D_hon <- D[honest_idx]
X_hon <- X[honest_idx, ]

## Estimate CATEs. Use training sample. 
library(grf)
forest <- causal_forest(X_tr, y_tr, D_tr) 
cates <- predict(forest, X)$predictions

## Construct and prune the tree. Also, compute node predictions (GATEs).
groupings <- build_aggtree(y, D, X, method = "aipw", cates = cates, is_honest = 1:length(y) %in% honest_idx)
```

By default, `build_aggtree` estimate CATEs internally via a [https://github.com/grf-labs/grf/blob/master/r-package/grf/R/causal_forest.R](causal forest). Alternatively, we can use the `cates` argument to input estimated CATEs. 



Once the tree is constructed, we can pick the grouping associated with our preferred granularity level. Then, we can estimate group average treatment effects in several ways. This package allows two estimation strategies. The first consists of taking the differences in mean outcomes between treated and control units in each group. This yields unbiased estimates if the assignment to treatment is randomized. Otherwise, in observational studies we can construct and average doubly-robust scores to get unbiased estimates.

To construct confidence intevals, we can regress the observed outcome on leaf dummies and their interactions with the treatment

We use different samples to estimate the CATEs (step 1) and construct the sequence of partitions (steps 2 and 3). We call these “estimation sample” and “aggregation sample.” To conduct valid inference, we need to estimate honest trees. For this, we need an additional “honest sample.”



  



```
## 1.) Estimate the CATEs. Use only estimation sample.
library(grf)

forest <- causal_forest(X[estimation_idx, ], y[estimation_idx], D[estimation_idx])
cates <- predict(forest, X)$predictions
```

To assess treatment effect heterogeneity, one can look at the distribution of the CATEs, e.g., via histograms or kernel density estimates. However, this is not conclusive evidence of heterogeneity: If the histogram is concentrated at one point, it may be that our estimator is not able to detect heterogeneity, and if the histogram is spread, it may be that our estimates are very noisy. 

A more systematic way to assess heterogeneity is to partition the population into groups that differ in the magnitude of their treatment effects. This is where aggregation trees play their role. 

```
## 2.) Construct a decision tree. Use only aggregation sample.
tree <- aggregation_tree(cates[aggregation_idx], X[aggregation_idx, ])

## 3.) Generate the sequence of partitions.
plot(tree, sequence = TRUE) # Get a visualization of the sequence.
```

To conduct valid inference, we can use the honest sample to regress `y` on leaf dummies and interactions of these dummies and `D`. As the treatment is randomized, the coefficients of the interaction terms are unbiased estimates of the average treatment effects in each group (GATEs):

```
## Get standard errors. Use only honest sample.
# First, we need a particular partition. 
subtree <- subtree_aggtree(tree, leaves = 5)

linear_model <- causal_ols_aggtree(subtree, y[honest_idx], X[honest_idx, ], D[honest_idx])
summary(linear_model)
```

Keep in mind that one should not conclude that covariates not used for splitting are not related to heterogeneity. There may exist several ways to form subpopulations
that differ in the magnitude of their treatment effects, and if two covariates are highly correlated, trees generally split on only one of those covariates.

A more systematic way to assess how treatment effects relate to the covariates consists of investigating how the average characteristics of the units vary across the leaves of the tree. For this, I provide a function that directly produces LATEX code for a nice table:

```
## Assess heterogeneity.
avg_characteristics_aggtree(subtree, X[aggregation_idx, ], cates = cates[aggregation_idx], method = "cates")
```

## References

- Athey, S., & Imbens, G. W. (2016).
<b>Recursive Partitioning for Heterogeneous Causal Effects.</b>
<i>Proceedings of the National Academy of Sciences</i>, 113(27).
[<a href="https://www.pnas.org/doi/abs/10.1073/pnas.1510489113">paper</a>]

- Athey, S., Tibshirani, J., & Wager, S. (2019).
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2).
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>]

- Chernozhukov, V., Demirer, M., Duflo, E., & Fernandez-Val, I. (2018).
<b>Generic Machine Learning Inference on Heterogeneous Treatment Effects in Randomized Experiments, with an Application to Immunization in India.</b>
<i>National Bureau of Economic Research</i>.
[<a href="https://www.nber.org/papers/w24678">paper</a>]

- Di Francesco, R. (2022).
<b>Aggregation Trees.</b> <i>CEIS Research Paper, 546.</i>
[<a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4304256">paper</a>]

- Wager, S., & Athey, S. (2018).
<b>Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.</b>
<i>Journal of the American Statistical Association</i>, 113(523).
[<a href="https://www.tandfonline.com/eprint/v7p66PsDhHCYiPafTJwC/full">paper</a>]
