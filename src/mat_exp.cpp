#include <RcppArmadillo.h>
// #include <Rcpp.h>
#include <iostream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

// gaussian kernel transform in cpp
// [[Rcpp::export]]
arma::mat gauss_cpp(const arma::mat& x, double sigma) {
  // arma::mat y = arma::exp(-arma::pow(x, 2) / (2 * sigma * sigma));
  // exp(-((x-mu)^2)/(2*sigma^2)) / (sigma*sqrt(2*pi))
  arma::mat y = arma::exp(-arma::square(x) / (2 * sigma * sigma)) / (sigma * sqrt(2 * M_PI));
  return y;
}

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

// compute euclidean distances in cpp
// [[Rcpp::export]]
arma::mat euclidean_cpp(const arma::mat& x, const arma::mat& y) {
  arma::mat d = arma::sqrt(arma::sum(arma::pow(x - y, 2), 1));
  return d;
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

// compute variance of each column of a matrix in cpp
// [[Rcpp::export]]
arma::mat varMat_cpp(const arma::mat& x) {
  arma::mat m = arma::var(x, 0);
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


// [[Rcpp::export]]
arma::mat sample_rows_cpp(List matrix_list, int n) {
  
  // combine all matrices in the list by row
  arma::mat combined_matrix;
  for (int i = 0; i < matrix_list.size(); i++) {
    arma::mat mat = matrix_list[i];
    if (i == 0) {
      combined_matrix = mat;
    } else {
      combined_matrix = join_rows(combined_matrix, mat);
    }
  }
  
  // sort each row in ascending order
  for (int i = 0; i < combined_matrix.n_rows; i++) {
    combined_matrix.row(i) = sort(combined_matrix.row(i));
  }
  
  // sample n elements within each row evenly
  arma::mat sampled_matrix(combined_matrix.n_rows, n);
  for (int i = 0; i < combined_matrix.n_rows; i++) {
    int step = combined_matrix.n_cols / (n-1);
    if (step == 0) {
      step = 1;
    }
    int count = 0;
    for (int j = 0; j < n; j++) {
      int idx = j * step;
      sampled_matrix(i, j) = combined_matrix(i, idx);
    }
  }
  
  return sampled_matrix;
}


// [[Rcpp::export]]
arma::sp_mat combine_sparse_rows(List sparse_matrix_list) {
  
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
arma::sp_mat combine_sparse_cols(List sparse_matrix_list) {
  
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
