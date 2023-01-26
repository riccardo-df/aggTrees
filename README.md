# Aggregation Trees
R package to implement aggregation trees, a nonparametric data-driven approach to discovering heterogeneous subgroups in a selection-on-observables framework. Additionally, the package provides useful functions to work with `rpart` objects.

The approach consists of three steps:

1. Estimate conditional average treatment effects (CATEs);
2. Aggregate CATEs into a decision tree;
3. Prune the tree.

This way, we generate a sequence of groupings, one for each granularity level. The resulting sequence is nested in the sense that subgroups formed at a given
level of granularity are never broken at coarser levels. This guarantees consistency of the results across the different granularity levels, generally considered a basic requirement that every classification system should satisfy. Moreover, each grouping features an optimality property in that it ensures that the loss in
explained heterogeneity resulting from aggregation is minimized.

Given the sequence of groupings, we can estimate group average treatment effects (GATEs) as we like. The package supports two estimators, based on differences in mean outcomes between treated and control units (unbiased in randomized experiments) and on sample averages of doubly-robust scores (unbiased also in observational studies). The package also allows to get standard errors for the GATEs by estimating via OLS appropriate linear models. An honesty condition is required to construct valid confidence intervals. Thus, different subsamples must be used to construct the tree and estimate the linear models. 

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
mu1 <- 0.5 * X[, 1] + X[, 2]
y <- mu0 + D * mu1 + rnorm(n)
```

As a first step, we need to estimate CATEs. We can do this with any estimator we like. Then, in the second step we construct a tree using the CATEs as an outcome. Given the tree, we can compute node predictions (i.e., GATEs) as we like. In the following chunk of code, I split the data into a training sample and an honest sample (required to conduct valid inference). We use the training sample to estimate CATEs and construct/prune the tree. Then, we use the honest sample to compute node predictions by constructing and averaging doubly-robust scores.

```
## Sample splitting.
splits <- sample_split(length(y), training_frac = 0.5)
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

## Construct and prune the tree using training sample. Also, compute node predictions (GATEs) using honest sample.
groupings <- build_aggtree(y, D, X, method = "aipw", cates = cates, is_honest = 1:length(y) %in% honest_idx)

## We have compatibility with generic S3-methods. 
summary(groupings)
print(groupings)
plot(groupings) # Beware of interpretation here. See below.

## To predict, do the following.
tree <- subtree(groupings$tree, cv = TRUE) # Select by cross-validation.
predict(tree, data.frame(X))
```

By default, `build_aggtree` estimate CATEs internally via a [causal forest](https://github.com/grf-labs/grf/blob/master/r-package/grf/R/causal_forest.R). Alternatively, we can override this by using the `cates` argument to input estimated CATEs, as I did above. When this is the case, we also need to specify `is_honest`, a logical vector denoting which observations we allocated to the honest sample. This way, `build_aggtree` knows which observations must be used to construct the tree and compute node predictions.

Now we have a whole sequence of optimal groupings. We can pick the grouping associated with our preferred granularity level and run some analysis. First, we would like to get standard errors for the GATEs. This is achieved by estimating via OLS appropriate linear models using the honest sample. Then, we would like to assess the driving factors of treatment effects by relating heterogeneity to observed covariates. Keep in mind that one should not conclude that covariates not used for splitting are not related to heterogeneity. There may exist several ways to form groups, and if two covariates are highly correlated, trees generally split on only one of those covariates. A more systematic way to assess how treatment effects relate to the covariates consists of investigating how the average characteristics of the units vary across groups. All of this is done in the following chunk of code:

```
## Analyze grouping with 5 groups.
results <- analyze_aggtree(groupings, n_groups = 5, method = "aipw", scores = groupings$scores)
summary(results$model)
```

`analyze_aggtree` prints LATEX code in the console. To avoid this, set `verbose = FALSE`. The code provides a table with GATEs and confidence intervals, and average characteristics of units in each leaf. This way, we obtain a nice and easy-to-read output that we can plug in papers/reports.

## References

- Athey, S., & Imbens, G. W. (2016).
<b>Recursive Partitioning for Heterogeneous Causal Effects.</b>
<i>Proceedings of the National Academy of Sciences</i>, 113(27).
[<a href="https://www.pnas.org/doi/abs/10.1073/pnas.1510489113">paper</a>]

- Athey, S., Tibshirani, J., & Wager, S. (2019).
<b>Generalized Random Forests.</b> <i>Annals of Statistics</i>, 47(2).
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>]

- Chernozhukov, V., Demirer, M., Duflo, E., & Fernandez-Val, I. (2017).
<b>Generic Machine Learning Inference on Heterogeneous Treatment Effects in Randomized Experiments.</b>
<i>National Bureau of Economic Research</i>.
[<a href="https://www.nber.org/papers/w24678">paper</a>]

- Di Francesco, R. (2022).
<b>Aggregation Trees.</b> <i>CEIS Research Paper, 546.</i>
[<a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4304256">paper</a>]

- Semenova, V., & Chernozhukov, V. (2021).
<b>Debiased Machine Learning of Conditional Average Treatment Effects and Other Causal Functions.</b>
<i>The Econometrics Journal</i>, 24 (2).
[<a href="https://academic.oup.com/ectj/article/24/2/264/5899048">paper</a>]

- Wager, S., & Athey, S. (2018).
<b>Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.</b>
<i>Journal of the American Statistical Association</i>, 113(523).
[<a href="https://www.tandfonline.com/eprint/v7p66PsDhHCYiPafTJwC/full">paper</a>]
