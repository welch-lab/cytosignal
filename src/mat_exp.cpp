#include <RcppArmadillo.h>
// #include <Rcpp.h>
#include <iostream>
#include <progress.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


// gaussian kernel transform for a vector in cpp
// [[Rcpp::export]]
arma::vec gauss_vec_cpp(const arma::vec& x, double sigma) {
  // arma::mat y = arma::exp(-arma::pow(x, 2) / (2 * sigma * sigma));
  // exp(-((x-mu)^2)/(2*sigma^2)) / (sigma*sqrt(2*pi))
  arma::vec y = arma::exp(-arma::square(x) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
  return y;
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


// compute euclidean distances in cpp
// [[Rcpp::export]]
arma::mat euclidean_cpp(const arma::mat& x, const arma::mat& y) {
  arma::mat d = arma::sqrt(arma::sum(arma::pow(x - y, 2), 1));
  return d;
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

// // compute variance of each column of a matrix in cpp
// // [[Rcpp::export]]
// arma::mat varMat_cpp(const arma::mat& x) {
//   arma::mat m = arma::var(x, 0);
//   return m;
// }


// // generate regular sequence in cpp
// // [[Rcpp::export]]
// arma::vec reg_seq_cpp(int n, int k) {
//   arma::vec x = arma::linspace(0, n*k, n+1);
//   return x;
// }

// // repeat vector by a given number of times
// // [[Rcpp::export]]
// arma::vec rep_vec_cpp(arma::vec x, int n) {
//   arma::vec y = arma::repmat(x, n, 1);
//   return y;
// }

// // repeat each element in a vector by a given number of times
// // [[Rcpp::export]]
// arma::vec rep_each_cpp(arma::vec x, int n) {
//   arma::vec y = arma::repelem(x, n, 1);
//   return y;
// }

// // gaussian kernel transform in cpp
// // [[Rcpp::export]]
// arma::mat gauss_kernel_cpp2(arma::mat x, double sigma) {
//   // arma::mat y = arma::exp(-arma::pow(x, 2) / (2 * sigma * sigma));
//   arma::mat y = arma::exp(-arma::pow(x, 2) / (2 * sigma * sigma));
//   return y;
// }

// // element-wise matrix exponential in cpp
// // [[Rcpp::export]]
// arma::mat expmat_cpp(arma::mat x) {
//   arma::mat y = arma::exp(x);
//   return y;
// }


// // [[Rcpp::export]]
// arma::sp_mat sp_sample_rows(List matrix_list, int n) {

//   // combine all matrices in the list by row
//   arma::sp_mat combined_matrix = arma::sp_join_rows(matrix_list);

//   // sort each row in ascending order
//   for (int i = 0; i < combined_matrix.n_rows; i++) {
//     combined_matrix.row(i) = arma::sort(combined_matrix.row(i));
//   }

//   // sample n elements within each row evenly
//   arma::sp_mat sampled_matrix(combined_matrix.n_rows, combined_matrix.n_cols);
//   for (int i = 0; i < combined_matrix.n_rows; i++) {
//     int nnz = combined_matrix.row(i).n_nonzero;
//     int step = nnz / n;
//     if (step == 0) {
//       step = 1;
//     }
//     int count = 0;
//     for (arma::sp_mat::const_row_iterator it = combined_matrix.begin_row(i); it != combined_matrix.end_row(i); ++it) {
//       if (count % step == 0) {
//         sampled_matrix(i, it.col()) = it.value();
//       }
//       count++;
//     }
//   }

//   return sampled_matrix;
// }


// // [[Rcpp::export]]
// arma::mat sample_rows_cpp(List matrix_list, int n) {

//   // combine all matrices in the list by row
//   arma::mat combined_matrix;
//   for (int i = 0; i < matrix_list.size(); i++) {
//     arma::mat mat = matrix_list[i];
//     if (i == 0) {
//       combined_matrix = mat;
//     } else {
//       combined_matrix = join_rows(combined_matrix, mat);
//     }
//   }

//   // sort each row in ascending order
//   for (int i = 0; i < combined_matrix.n_rows; i++) {
//     combined_matrix.row(i) = sort(combined_matrix.row(i));
//   }

//   // sample n elements within each row evenly
//   arma::mat sampled_matrix(combined_matrix.n_rows, n);
//   for (int i = 0; i < combined_matrix.n_rows; i++) {
//     int step = combined_matrix.n_cols / (n-1);
//     if (step == 0) {
//       step = 1;
//     }
//     int count = 0;
//     for (int j = 0; j < n; j++) {
//       int idx = j * step;
//       sampled_matrix(i, j) = combined_matrix(i, idx);
//     }
//   }

//   return sampled_matrix;
// }


// // gaussian kernel transform in cpp
// // [[Rcpp::export]]
// arma::mat gauss_cpp(const arma::mat& x, double sigma) {
//   // arma::mat y = arma::exp(-arma::pow(x, 2) / (2 * sigma * sigma));
//   // exp(-((x-mu)^2)/(2*sigma^2)) / (sigma*sqrt(2*pi))
//   arma::mat y = arma::exp(-arma::square(x) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
//   return y;
// }

// [[Rcpp::export]]
arma::uvec select_EB_singleSpot_rcpp(
    const arma::mat& loc,
    const double x_i,
    const double y_i,
    const double eps) {
  arma::colvec x = loc.col(0);
  arma::colvec y = loc.col(1);
  double xmax = x_i + eps;
  double xmin = x_i - eps;
  double ymax = y_i + eps;
  double ymin = y_i - eps;
  arma::uvec idx = arma::find((x < xmax) && (x > xmin) && (y < ymax) && (y > ymin));
  return idx;
}

// x - the location matrix, rows - spots, columns x and y
// eps - radius of the selection
// [[Rcpp::export]]
Rcpp::List select_EB_rcpp(const arma::mat& loc, const double eps) {
  int n = loc.n_rows;
  arma::colvec x = loc.col(0);
  arma::colvec y = loc.col(1);
  arma::colvec x_box, y_box;
  double xmax, xmin, ymax, ymin;
  arma::uvec idx, subidx;
  Rcpp::List id, dist;
  Rcpp::Rcerr << "- GauEps: Selecting neighbors within the radius of " << eps <<
    " from " << n << " spots" << std::endl;
  Progress p(n, true);
  for (int i = 0; i < n; i++) {
    if (Progress::check_abort()) {
      return R_NilValue;
    }
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
    dists = dists(subidx);
    id.push_back(Rcpp::IntegerVector(idx.begin(), idx.end()) + 1);
    dist.push_back(Rcpp::NumericVector(dists.begin(), dists.end()));
    p.increment();
  }
  return Rcpp::List::create(Rcpp::Named("id") = id, Rcpp::Named("dist") = dist);
}
