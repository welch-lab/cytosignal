#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cli/progress.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


// x - the location matrix, rows - spots, columns x and y
// eps - radius of the selection
// Rcpp::List select_EB_rcpp(const arma::mat& loc, const double eps) {
// [[Rcpp::export]]
arma::sp_mat select_EB_rcpp2(const arma::mat& loc, const double eps) {
    int n = loc.n_rows;
    arma::colvec x = loc.col(0);
    arma::colvec y = loc.col(1);
    arma::colvec x_box, y_box;
    double xmax, xmin, ymax, ymin;
    arma::uvec idx, subidx;
    arma::sp_mat dist_mat(n, n);
    SEXP bar = PROTECT(cli_progress_bar(n, NULL));
    for (int i = 0; i < n; i++) {
        Rcpp::checkUserInterrupt();
        double x_i = x(i);
        double y_i = y(i);
        xmax = x_i + eps;
        xmin = x_i - eps;
        ymax = y_i + eps;
        ymin = y_i - eps;
        // filter for x < xmax and get the index
        idx = arma::find(x < xmax);
        idx = idx(arma::find(x(idx) > xmin));
        idx = idx(arma::find(y(idx) < ymax));
        idx = idx(arma::find(y(idx) > ymin));
        idx = idx(arma::find(idx != i));
        x_box = x(idx);
        y_box = y(idx);
        arma::vec dists = arma::sqrt(arma::pow(x_box - x_i, 2) + arma::pow(y_box - y_i, 2));
        subidx = arma::find(dists < eps);
        idx = idx(subidx);
        if (idx.n_elem == 0) {
            if (CLI_SHOULD_TICK) cli_progress_set(bar, i);
            continue;
        }
        dists = dists(subidx);
        for (int j = 0; j < idx.n_elem; j++) {
            dist_mat(idx(j), i) = dists(j);
        }
        if (CLI_SHOULD_TICK) cli_progress_set(bar, i);
    }
    cli_progress_done(bar);
    UNPROTECT(1);
    return dist_mat;
}

// Inplace gaussian kernel transform for all non-zero elements in a sparse matrix
// [[Rcpp::export]]
void gauss_vec_inplace_cpp(NumericVector& x, const double sigma) {
    x = exp(-pow(x, 2) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
}
