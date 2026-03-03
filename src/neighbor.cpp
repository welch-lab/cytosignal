#include <RcppArmadillo.h>
#include <cli/progress.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

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


// X - n by 2 matrix of coordinates
// radius - distance threshold
// ncores - number of threads to use for parallel processing
// returns a sparse matrix of size n by n, where non-zero (i,j) entries indicate
// that spot i is within radius threshold of spot j, and the value of the entry
// is the distance between the two spots.
// [[Rcpp::export]]
arma::sp_mat dist_mat_within_r(
        const arma::mat& X,
        const double radius,
        const int ncores = 1
) {
    const arma::uword n = X.n_rows;
    const double R2 = radius * radius;
    const arma::colvec x = X.col(0);
    const arma::colvec y = X.col(1);
    std::vector< std::vector<arma::uword> > col_rowIndices(n);
    std::vector< std::vector<double> > col_values(n);
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores) schedule(dynamic)
#endif
    for (arma::uword i = 0; i < n; i++) {
        // let 'i' be the center spot.
        // Observation for each i-spot recorded in a column of the sparse matrix.
        double x_i = X(i, 0);
        double y_i = X(i, 1);
        
        for (arma::uword j = 0; j < n; j++) {
            // let 'j' be the neighbor spot.
            // Record the distance and the index of in a i-sub-vector.
            double x_j = X(j, 0);
            double y_j = X(j, 1);
            double dx = x_i - x_j;
            // early exits if x/y distance is already greater than radius
            if (dx > radius) continue;
            if (dx < -radius) continue;
            double dy = y_i - y_j;
            if (dy > radius) continue;
            if (dy < -radius) continue;
            double d2 = dx * dx + dy * dy;
            if (d2 <= R2) {
                double dist = sqrt(d2);
                col_rowIndices[i].push_back(j);
                col_values[i].push_back(dist);
            }
        }
    }
    // Now assemble the CSC sparse matrix
    std::vector<arma::uword> row_indices;
    std::vector<arma::uword> col_ptrs(n + 1, 0);
    std::vector<double> values;
    for (arma::uword i = 0; i < n; i++) {
        for (arma::uword j = 0; j < col_rowIndices[i].size(); j++) {
            row_indices.push_back(col_rowIndices[i][j]);
            values.push_back(col_values[i][j]);
        }
        col_ptrs[i + 1] = row_indices.size();
    }
    arma::uvec row_indices_uvec(row_indices);
    arma::uvec col_ptrs_uvec(col_ptrs);
    arma::vec values_vec(values);
    return arma::sp_mat(row_indices_uvec, col_ptrs_uvec, values_vec, n, n);
}
