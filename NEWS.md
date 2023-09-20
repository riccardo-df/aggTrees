## version 2.0.2

Fixed a bug that occurred when using only one covariate.

## version 2.0.1

Now the inference_aggtree function can estimate asymmetric bias-corrected and accelerated 95% confidence intervals for the GATEs based on their bootstrap distributions.

Also, improved automatic check of leaf composition. Functions that estimate the GATEs now:

-   Raise an error if one or more leaves contain only one observation;
-   Raise an error if differences in mean outcomes is the chosen estimator and one or more leaves contain only treated or only control units.

## version 2.0.0

Updated hypothesis testing procedure. Improved LATEX output.

## version 1.0.0

First public release
