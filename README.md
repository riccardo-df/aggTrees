# Aggregation Trees <a href="https://riccardo-df.github.io/aggTrees/"><img src="man/figures/logo.svg" align="right" height="200" /></a>

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT) 
[![CRAN status](https://www.r-pkg.org/badges/version/aggTrees)](https://CRAN.R-project.org/package=aggTrees) 
[![Downloads](https://cranlogs.r-pkg.org/badges/aggTrees)](https://CRAN.R-project.org/package=aggTrees)

`aggTrees` is an R package for *aggregation trees*, a nonparametric approach designed to discover heterogeneous subgroups in a selection-on-observables framework. This method helps researchers assess treatment effect heterogeneity by constructing a hierarchical sequence of optimal groupings, enabling estimation and inference at different levels of granularity.

------------------------------------------------------------------------

### Key features
✔ **Data-driven subgroup discovery** for causal inference.
✔ **Optimal partitions** for identifying treatment effect heterogeneity.
✔ **Point estimation and inference** for the average treatment effect of each subgroup.
✔ **Fully nonparametric** — no need to pre-specify effect moderators.
✔ **Seamless integration** with standard econometric workflows.
✔ **Active development & support** - open-source and actively maintained.                                                          |

------------------------------------------------------------------------

## 🚀 Installation

To install the latest stable version from CRAN:

```         
install.packages("aggTrees")
```

Alternatively, the current development version of the package can be installed using the `devtools` package:

```         
devtools::install_github("riccardo-df/aggTrees") # run install.packages("devtools") if needed.
```

------------------------------------------------------------------------

## Contributing

We welcome contributions! If you encounter issues, have feature requests, or want to contribute to the package, please follow the guidelines below.

📌 **Report an issue:** If you encounter a bug or have a suggestion, please open an issue on GitHub: 
[Submit an issue](https://github.com/riccardo-df/aggTrees/issues)

📌 **Contribute code:** We encourage contributions via pull requests. Before submitting, please:
1. Fork the repository and create a new branch.
2. Ensure that your code follows the existing style and documentation conventions.
3. Run tests and check for package integrity.
4. Submit a pull request with a clear description of your changes.

📌 **Feature requests:** If you have ideas for new features or extensions, feel free to discuss them by opening an issue.

------------------------------------------------------------------------

## Citation

If you use `aggTrees` in your research, please cite the corresponding paper:

> Di Francesco, R. (2024). Aggregation trees. arXiv preprint arXiv:2410.11408.

