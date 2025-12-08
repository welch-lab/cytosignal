#include <RcppArmadillo.h>
// #include <Rcpp.h>
#include <iostream>
#include <progress.hpp>
#include <cli/progress.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


// gaussian kernel transform for a vector in cpp
// [[Rcpp::export]]
NumericVector gauss_vec_cpp(NumericVector x, double sigma) {
  x = exp(-pow(x, 2) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
  return x;
}

// [[Rcpp::export]]
arma::vec rep_each_cpp(int l, int n) {
  arma::vec x = arma::linspace(1, l, l);
  arma::vec y = arma::repelem(x, n, 1);
  return y;
}


// average a list of matrices in rcpp
// [[Rcpp::export]]
arma::sp_mat meanMat_cpp(const List& l, int nrow, int ncol) {
  arma::sp_mat m;
  m.zeros(nrow, ncol);

  for (int i = 0; i < l.size(); i++) {
    // m += l[i];
    arma::sp_mat x = l[i];
    m = m + x;
  }
  // m /= l.size();
  m = m / l.size();
  return m;
}


// compute standard deviation of each column of a matrix in cpp
// [[Rcpp::export]]
arma::mat stdMat_cpp(const arma::mat& x) {
  arma::mat m = arma::stddev(x, 0);
  return m;
}

// compute pearson correlation between correpsonding columns of two matrices in cpp
// [[Rcpp::export]]
arma::vec pearson_col_cpp(const arma::mat& x, const arma::mat& y) {
  arma::vec res(x.n_cols);
  for (int i = 0; i < x.n_cols; i++) {
    res.subvec(i,i) = arma::cor(x.col(i), y.col(i), 0);
  }
  return res;

}


// [[Rcpp::export]]
arma::sp_mat cbind_list(List& sparse_matrix_list) {

  // combine all sparse matrices in the list by row
  arma::sp_mat combined_sparse_matrix;
  for (int i = 0; i < sparse_matrix_list.size(); i++) {
    arma::sp_mat mat = sparse_matrix_list[i];
    if (i == 0) {
      combined_sparse_matrix = mat;
    } else {
      combined_sparse_matrix = arma::join_rows(combined_sparse_matrix, mat);
    }
  }

  return combined_sparse_matrix;
}


// [[Rcpp::export]]
arma::sp_mat rbind_list(List sparse_matrix_list) {

  // combine all sparse matrices in the list by column
  arma::sp_mat combined_sparse_matrix;
  for (int i = 0; i < sparse_matrix_list.size(); i++) {
    arma::sp_mat mat = sparse_matrix_list[i];
    if (i == 0) {
      combined_sparse_matrix = mat;
    } else {
      combined_sparse_matrix = arma::join_cols(combined_sparse_matrix, mat);
    }
  }

  return combined_sparse_matrix;
}

// Clean up lrscore sparse matrices by forcing <1 values to 0, and round to
// integer numbers. Finally clean up sparsity structure.
// [[Rcpp::export]]
arma::sp_mat cleanLRscore_sparse_cpp(
    arma::uvec i,
    arma::uvec p,
    arma::colvec x,
    const int nrow,
    const int ncol
) {
  x.elem(arma::find(x < 1)).zeros();
  x = arma::round(x);
  arma::sp_mat out(i, p, x, nrow, ncol, true);
  return out;
}

