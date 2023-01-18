#include <Rcpp.h>
using namespace Rcpp;

//' Honest rpart Estimates
//'
//' Computes honest estimates for a rpart object.
//'
//' @param unique_leaves_honest Vector storing the unique leaf ids of the tree relative to the honest sample.
//' @param y_honest Outcome vector of honest units. This could be any vector. Leaf are replaced with means of this object.
//' @param honest_leaves Vector of size \code{n.samples}. The i-th element stores the id of the leaf where the i-th honest observation falls.
//'
// [[Rcpp::export]]
NumericVector honest_rpart_cpp(List unique_leaves_honest, NumericVector y_honest, NumericVector honest_leaves) { // Taken from https://github.com/okasag/orf/blob/master/orf/src/pred_honest_rcpp.cpp

  // Declaring variables.
  int leaf_ID = 0;

  int n_rows_honest = honest_leaves.size();

  NumericVector obs_same_all_honest(n_rows_honest);

  NumericVector y(n_rows_honest);
  NumericVector y_same(n_rows_honest);

  double y_mean = 0;

  NumericVector honest_pred(n_rows_honest);

  int n_leaves = unique_leaves_honest.size();

  // Looping over the leaves where honest observations fall.
  for(int leaf_idx = 0; leaf_idx < n_leaves; ++leaf_idx) {
    // Select one leaf.
    leaf_ID = unique_leaves_honest[leaf_idx];

    // Looping over honest units. Finds observations that fall in this leaf.
    for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
      obs_same_all_honest[row_idx] = honest_leaves(row_idx) == leaf_ID; // True if honest unit falls into this leaf.
    }

    // If honest unit in this leaf, save associated outcome, otherwise write a zero.
    for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
      if (obs_same_all_honest[row_idx] == 1) {
        y_same[row_idx] = y_honest[row_idx];
      } else {
        y_same[row_idx] = 0;
      }
    }

    double y_count = 0;
    double y_sum = 0;

    // Count and sum honest outcomes within this leaf.
    for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
      if (obs_same_all_honest[row_idx] == 1) {
        y_count = y_count + obs_same_all_honest[row_idx]; // Due to above loop, here we sum 1 if the row_idx-th honest unit is in this leaf, and 0 otherwise.
        y_sum = y_sum + y_same[row_idx]; // Due to the above loop, here we sum outcomes if the row_idx-th honest unit is in this leaf, and 0 otherwise.
      }
    }

    y_mean = y_sum / y_count;

    // If row_idx-th honest unit is in this leaf, assign this prediction.
    for(int row_idx = 0; row_idx < n_rows_honest; ++row_idx) {
      if (obs_same_all_honest[row_idx] == 1) {
        honest_pred(row_idx) = y_mean;
      }
    }
  }

  return honest_pred;
}
