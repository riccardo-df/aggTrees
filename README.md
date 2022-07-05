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

## Example
This section shows how to use aggregation trees to discover heterogeneous subpopulations using a random subsample of the data on maternal smoking and infants' weight at birth (available as a built-in data set in Stata and in the [GitHub repository](https://github.com/mdcattaneo/replication-C_2010_JOE) of Matias Cattaneo).  

Let's load the data and define the relevant variables.

```
dta <- haven::read_dta("http://www.stata-press.com/data/r13/cattaneo2.dta")

Y <- as.matrix(dta[, "bweight"])
D <- as.matrix(dta[, "mbsmoke"])
X_names <- c("bweight", "mbsmoke", "msmoke", "deadkids", "monthslb", "lbweight")
X <- as.matrix(dta[, !(colnames(dta) %in% X_names)])
```

For our analysis, we split the full sample into an **estimation sample** and an **aggregation sample**. This is a standard procedure: treatment effects are estimated in one part of the sample, and heterogeneity is assessed in the held-out sample.

```
set.seed(1986)
n <- dim(dta)[1]
est_idx <- sample(1:n, n / 2, replace = FALSE)

X_est <- X[est_idx, ]
Y_est <- Y[est_idx]
D_est <- D[est_idx]

X_agg <- X[-est_idx, ]
Y_agg <- Y[-est_idx]
D_agg <- D[-est_idx]    
```   

Building an aggregation tree is a three-step methodology.

1. **estimation step** -> CATEs are estimated using only the estimation sample, and the fitted model predicts out-of-sample in the aggregation sample.
2. **tree-growing step** -> The out-of-sample predictions are used to construct a decision tree that provides the set of admissible stratifications of the population.
3. **tree-pruning step** -> The tree from the previous step is pruned, thus producing a sequence of subtrees. Each subtree provides an "optimal" partition of the population, one for each granularity level.[^1]

[^1]: Notice that the tree-growing and the tree-pruning steps are separated because the order in which nodes are formed in the former needs not to correspond to the order in which nodes are collapsed in the latter.
  
### Estimation step
Without further ado, let's jump to the first step: estimation of the CATEs. Notice that aggregation trees are agnostic about the choice of the CATE estimator. Here, we use the [causal forest](https://github.com/grf-labs/grf) estimator, which is simple to understand and efficiently coded. The first line in the following block fits the forest in the estimation sample, and the second line predicts the CATEs on the held-out aggregation sample. The estimated density of the out-of-sample predictions is then displayed using the `ggplot2` package.

```
cates_forest <- grf::causal_forest(X = X_est, Y = Y_est, W = D_est)
cates <- predict(cates_forest, newdata = X_agg)$predictions

library(ggplot2)
ggplot(data.frame(x = cates), aes(x = x)) +
  geom_density(fill = "steelblue", color = "black", alpha = 0.8, adjust = 1) +
  xlab("CATEs") + ylab("Density") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```

![density](/images/cates_density_github.png)

From the plot above, we notice substantial heterogeneity in the treatment effects. Effects range from -382.681 grams to -57.087 grams (you can double-check using `range(cates)`). The mean of the distribution is an estimate of the ATE and corresponds to -257.681 grams. The standard deviation of the distribution is 63.317.  

Fitting a regression forest using the estimated CATEs and computing the importance of each variable allows us to gain some knowledge on which covariates are more related to treatment effect heterogeneity. This is the role of the `var_importance` function. Notice that the forest is grown in the aggregation sample. Results can be displayed by using the `plot_importance` function. To avoid overplotting, it is possible to show only the `k` most important covariates.

```
var_imp <- var_importance(cates, X_agg)
plot_importance(var_imp, k = 10)
```

![importance](/images/importance_github.png)

Results show that the mother's age and the number of prenatal care visits are the main drivers of treatment effect heterogeneity. However, nothing can be inferred about the direction of the effect of these two covariates.

### Tree-growing step
Further insights into treatment effect heterogeneity can be learned from analyzing intermediate aggregation levels. The idea is to discover the subpopulations most and least impacted by the treatment. This is the goal of aggregation trees.  

To build the tree, we use the `aggregation_tree` function, which allows the user to set two stopping criteria:
- `maxdepth` -> The maximum depth of the tree. My suggestion is to keep the tree shallow enough to avoid overplotting, but at the same time grow a tree deep enough to explore finer granularity levels (generally, `maxdepth = 3` or `maxdepth = 4` is a good choice).
- `cp` -> Minimum decrease in the mean squared error of the model that the next split must achieve to be accepted. 

The `plot_aggregation_tree` function asks the user to select a palette for the tree. I suggest using diverging colors so that it is immediate to notice the most and the least affected groups. The plot below uses the blue-red palette from the *Advance: Diverging* type, clickable from the drop-down menu in the upper part of the `colorspace::choose_palette()` window.[^2]

[^2]: If you need to build several trees, it can be annoying selecting the palette each time. The package allows the user to avoid this by first choosing a palette `palette <- colorspace::choose_palette()` and the using it as an input in the `plot_aggregation_tree` function: `plot_aggregation_tree(tree, palette, main = "My example tree"`).

```
tree <- aggregation_tree(cates, X_agg, maxdepth = 3, cp = 0.01)
plot_aggregation_tree(tree, main = "My example tree")
```

![tree](/images/tree_github.png)

The figure above displays the estimated aggregation tree. Associated with each node is a group of units, and split labels provide information about the identity of these groups. The number and percentage of units belonging to the associated group are displayed within each node, together with the point estimate of the corresponding GATE. Blue and red shades denote groups with stronger (i.e., more negative) and lighter (i.e., more positive) GATEs than the ATE.

The aggregation tree shows that treatment effects are increasingly negative with the mother's age and the number of prenatal care visits.[^3] At the finest granularity level we explore, we have six distinct subpopulations. The GATEs for the most and the least impacted units are -311 grams and -130 grams.

[^3]: In the original paper, where I use a more comprehensive data set (435,125 observations and 39 covariates), the number of prenatal care visits does not appear to substantially drive treatment effect heterogeneity. Due to the nature of the data used here, results in the paper are more reliable.

### Tree-pruning step
To explore coarser aggregation levels, the tree-pruning step prunes the tree to derive the whole sequence of "optimal" partitions of the population. To visualize this sequence, it is sufficient to set `sequence = TRUE` in the `plot_aggregation_tree` function.

```
plot_aggregation_tree(tree, main = "My example tree", sequence = TRUE)
```

### Building the tree on a subset of the covariates
One of the main advantages of aggregation trees is that they allow the researcher to study heterogeneity in the treatment effects relative to a subset
of the covariates used in the estimation step. In this way, it is possible to focus the heterogeneity analysis only on the dimensions of interest.  
  
In this subsection, we use only the mother's age (`mage`) and the number of prenatal care visits (`nprenatal`) to build our tree. As we already estimated the CATEs, there is no need to repeat the estimation step. To keep things simple and easy to interpret, I set `maxdepth = 2`.

```
X_subsample <- as.data.frame(X_agg[, c("mage", "nprenatal")])
tree <- aggregation_tree(cates, X_subsample, maxdepth = 2, cp = 0.01)
```

Because the tree is constructed using only two covariates, we can represent it by plotting the axis-aligned splits it performed. This is achieved by the `recursive_partitioning_plot` function, which overlays the splits on a scatter plot of the covariates used in the tree construction. To aid readability, I suggest coloring the points according to the magnitude of the CATEs. By default, `recursive_partitioning_plot` uses yellow for more negative effects and red for more positive effects. Users can override this choice by setting the parameters `low` and `high`. Because in the present exercise effects are always negative, I decide to set `low = "red"` and `high = "yellow"`. 

```
plot <- recursive_partitioning_plot(tree, cates, X_subsample, low = "red", high = "yellow", size = 3.5)
plot
```

![recursive](/images/recursive_github.png)

The tree discovers four distinct subpopulations formed according to the mother's age and the number of prenatal care visits. Thanks to the coloring, it is immediate to notice that the most affected group comprises children born to mothers older than 22 who attended at least eight prenatal care visits. As before, treatment effects are increasingly negative with both variables.

### Comparison with clustering methods
Finally, we can compare the results from the aggregation tree to those obtained by cluster analysis. For this purpose, we form three clusters using the estimated CATEs and compute the average value of our covariates for each cluster. Again, as we are using only two covariates in this example, we can plot the average characteristics of our clusters. Average values are computed under the hood by the `cluster_plot` function.

```
clusters <- kmeans(cates, 3)

plot <- cluster_plot(clusters, cates, X_subsample, low = "red", high = "yellow", size = 3.5)
plot
```

![cluster](/images/cluster_github.png)

The cluster analysis results align with the tree's findings: mothers of the most impacted cluster are, on average, older than those of the least impacted cluster and attended a higher number of prenatal care visits, although the difference in the average values of `nprenatal` is not substantial.

As a last remark, the package allows the user to overlay the recursive partitioning of the aggregation tree and the results from cluster analysis. This is achieved by first generating the tree plot and then passing it as an argument in the `cluster_plot` function. 

```
tree_plot <- recursive_partitioning_plot(tree, cates, X_subsample, low = "red", high = "yellow", size = 3.5)
cluster_plot <- cluster_plot(clusters, cates, X_subsample, plot = tree_plot, low = "red", high = "yellow", size = 3.5)
cluster_plot
```

![comparison](/images/comparison_github.png)
