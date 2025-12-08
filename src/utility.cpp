#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <progress.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

// normalize each column of a sparse matrix in cpp
// [[Rcpp::export]]
arma::sp_mat normalizeSparse_cpp(const arma::sp_mat& x) {
      arma::sp_mat normalized_x = x;
      arma::vec colsum = arma::zeros(x.n_cols);
      // Use sparse matrix pointer to accumulate the colsums
      for (arma::sp_mat::const_iterator it = x.begin(); it != x.end(); ++it) {
         colsum[it.col()] += *it;
      }
      for (arma::sp_mat::iterator it = normalized_x.begin(); it != normalized_x.end(); ++it) {
          if (colsum[it.col()] != 0) {
              *it /= colsum[it.col()];
          }
      }
      return normalized_x;
}


// compute euclidean distances in cpp
// [[Rcpp::export]]
NumericVector euclidean_elementwise_cpp(const arma::mat& x, const arma::mat& y) {
   arma::mat d = arma::sqrt(arma::sum(arma::pow(x - y, 2), 1));
   // copy the only column to a vector
   NumericVector d_vec(d.begin(), d.end());
   return d_vec;
}


void print_sp_mat_header(const arma::sp_mat& mat, arma::uword nrows = 5, arma::uword ncols = 5) {
    arma::uword rows = std::min(nrows, mat.n_rows);
    arma::uword cols = std::min(ncols, mat.n_cols);
    for (arma::uword i = 0; i < rows; ++i) {
        for (arma::uword j = 0; j < cols; ++j) {
            if (mat(i, j) != 0) {
                std::cout << mat(i, j) << "\t";
            } else {
                std::cout << ".\t";
            }
        }
        std::cout << std::endl;
    }
}
