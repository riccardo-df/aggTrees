# Aggregation Trees
R package to implement aggregation trees, a nonparametric approach to discovering heterogeneous subgroups in a selection-on-observables framework. 

The approach consists of three steps:

1. Estimate the conditional average treatment effects (CATEs);
2. Approximate the CATEs by a decision tree;
3. Prune the tree.

This way, we generate a sequence of groupings, one for each granularity level. The resulting sequence is nested in the sense that subgroups formed at a given
level of granularity are never broken at coarser levels. This guarantees consistency of the results across the different granularity levels, generally considered a basic requirement that every classification system should satisfy. Moreover, each grouping features an optimality property in that it ensures that the loss in
explained heterogeneity resulting from aggregation is minimized.

Given the sequence of groupings, we can estimate the group average treatment effects (GATEs) as we like. The package supports two estimators, based on differences in mean outcomes between treated and control units (unbiased in randomized experiments) and on sample averages of doubly-robust scores (unbiased also in observational studies). The package also allows to get standard errors for the GATEs by estimating via OLS appropriate linear models. An honesty condition is required to conduct valid inference. Thus, different subsamples must be used to construct the tree and estimate the linear models. 

## Installation  
The package can be downloaded from CRAN:

```
install.packages("aggTrees")
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```
devtools::install_github("riccardo-df/aggTrees") # run install.packages("devtools") if needed.
```

## Usage 
Please chek the package vignette for a short tutorial.

## References

- Athey, S., & Imbens, G. W. (2016).
<b>Recursive Partitioning for Heterogeneous Causal Effects.</b>
<i>Proceedings of the National Academy of Sciences</i>, 113(27).
[<a href="https://www.pnas.org/doi/abs/10.1073/pnas.1510489113">paper</a>]

- Athey, S., Tibshirani, J., & Wager, S. (2019).
<b>Generalized Random Forests.</b> 
<i>Annals of Statistics</i>, 47(2).
[<a href="https://projecteuclid.org/euclid.aos/1547197251">paper</a>]

- Chernozhukov, V., Demirer, M., Duflo, E., & Fernandez-Val, I. (2017).
<b>Generic Machine Learning Inference on Heterogeneous Treatment Effects in Randomized Experiments.</b>
<i>National Bureau of Economic Research</i>.
[<a href="https://www.nber.org/papers/w24678">paper</a>]

- Cotterman, R., & Peracchi, F. (1992).
<b>Classification and aggregation: An application to industrial classification in cps data.</b> 
<i>Journal of Applied Econometrics</i>, 7(1).
[<a href="https://onlinelibrary.wiley.com/doi/abs/10.1002/jae.3950070105">paper</a>]

- Di Francesco, R. (2022).
<b>Aggregation Trees.</b> 
<i>CEIS Research Paper, 546.</i>
[<a href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4304256">paper</a>]

- Holm, S. (1979).
<b>A Simple Sequentially Rejective Multiple Test Procedure.</b> 
<i>Scandinavian Journal of Statistics</i>, 6(2).
[<a href="https://www.jstor.org/stable/4615733">paper</a>] 

- Semenova, V., & Chernozhukov, V. (2021).
<b>Debiased Machine Learning of Conditional Average Treatment Effects and Other Causal Functions.</b>
<i>The Econometrics Journal</i>, 24(2).
[<a href="https://academic.oup.com/ectj/article/24/2/264/5899048">paper</a>]
