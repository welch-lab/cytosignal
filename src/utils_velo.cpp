#include <RcppArmadillo.h>
// #include <Rcpp.h>
#include <iostream>
// #include <vector>
// #include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat inferScoreLR_cpp(
    const arma::mat& dge_lig,
    const arma::mat& dge_recep,
    const arma::uvec& lig_index,
    const arma::uvec& lig_list,
    const arma::uvec& recep_index,
    const arma::uvec& recep_list
){
    // lig_dge & recep_dge: genes X beads
    // res_mtx: bead X interactions mtx
    arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

    for (uword i = 0; i < res_mtx.n_cols; i++){ // for each intr

        arma::mat lig_dge = dge_lig.rows(lig_list(span(lig_index(i), lig_index(i+1)-1)));
        arma::rowvec lig_vec = arma::sum(lig_dge, 0); // column sum of lig_dge, 1 X n_beads

        arma::mat recep_dge = dge_recep.rows(recep_list(span(recep_index(i), recep_index(i+1)-1)));
        arma::rowvec recep_vec = arma::sum(recep_dge, 0); // column sum of recep_dge, 1 X n_beads

        // arma::rowvec intr_res = lig_vec % recep_vec;
        res_mtx.col(i) = (lig_vec % recep_vec).t();
    }

    return res_mtx;
}



// [[Rcpp::export]]
arma::mat inferVeloLR_cpp(
    const arma::mat& dge_lig,
    const arma::mat& dge_recep,
    const arma::mat& dge_lig_velo,
    const arma::mat& dge_recep_velo,
    const arma::uvec& lig_index,
    const arma::uvec& lig_list,
    const arma::uvec& recep_index,
    const arma::uvec& recep_list
){
    // lig_dge & recep_dge: genes X beads
    // res_mtx: bead X interactions mtx
    arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

    for (uword i = 0; i < res_mtx.n_cols; i++){ // for each intr

        arma::mat lig_dge = dge_lig.rows(lig_list(span(lig_index(i), lig_index(i+1)-1)));
        arma::rowvec lig_vec = arma::sum(lig_dge, 0); // column sum of lig_dge, 1 X n_beads

        arma::mat lig_velo = dge_lig_velo.rows(lig_list(span(lig_index(i), lig_index(i+1)-1)));
        arma::rowvec lig_velo_vec = arma::sum(lig_velo, 0); // column sum of lig_dge, 1 X n_beads

        arma::mat recep_dge = dge_recep.rows(recep_list(span(recep_index(i), recep_index(i+1)-1)));
        arma::rowvec recep_vec = arma::sum(recep_dge, 0); // column sum of recep_dge, 1 X n_beads

        arma::mat recep_velo = dge_recep_velo.rows(recep_list(span(recep_index(i), recep_index(i+1)-1)));
        arma::rowvec recep_velo_vec = arma::sum(recep_velo, 0); // column sum of recep_dge, 1 X n_beads

        // arma::rowvec intr_res = lig_vec % recep_vec;
        // res_mtx.col(i) = (lig_vec % recep_vec).t();
        res_mtx.col(i) = (lig_vec % recep_velo_vec + lig_velo_vec % recep_vec).t();
    }

    return res_mtx;
}



// [[Rcpp::export]]
arma::mat graphNicheLR_cpp(
    const arma::mat& dge_lig,
    const arma::mat& dge_recep,
    const arma::uvec& lig_index,
    const arma::uvec& lig_list,
    const arma::uvec& recep_index,
    const arma::uvec& recep_list,
    const arma::uvec& nb_index,
    const arma::uvec& nb_list
){
    // lig_dge & recep_dge: genes X beads
    // res_mtx: bead X interactions mtx
    arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

    // for (i in ncol){ for (j in nrow) {...}}
    for (uword i = 0; i < res_mtx.n_cols; i++){
        for (uword j = 0; j < res_mtx.n_rows; j++){
            // All indices in Cpp starts with 0 --> all values -1
            // mean(x, dim): each column (dim = 0), or each row (dim = 1)
            arma::mat lig_dge = dge_lig.submat( lig_list(span(lig_index(i), lig_index(i+1)-1)), 
                                            nb_list(span(nb_index(j), nb_index(j+1)-1)) );
            // arma::vec lig_vec = mean(lig_dge, 1); // n_ligs X 1
            double lig_val = arma::accu(mean(lig_dge, 1)); // n_ligs X 1

            arma::mat recep_dge = dge_recep.submat( recep_list(span(recep_index(i), recep_index(i+1)-1)), 
                                            nb_list(span(nb_index(j), nb_index(j+1)-1)) );
            // arma::rowvec recep_vec = mean(recep_dge, 1).t(); // 1 X n_receps
            double recep_val = arma::accu(mean(recep_dge, 1)); // 1 X n_receps

            // res_mtx(j, i) = arma::accu(lig_vec * recep_vec);
            res_mtx(j, i) = lig_val*recep_val;
        }
    }

    return res_mtx;
}



// [[Rcpp::export]]
arma::mat VelographNicheLR_cpp(
    const arma::mat& dge_lig,
    const arma::mat& dge_recep,
    const arma::mat& dge_velo,
    const arma::uvec& lig_index,
    const arma::uvec& lig_list,
    const arma::uvec& recep_index,
    const arma::uvec& recep_list,
    const arma::uvec& nb_index,
    const arma::uvec& nb_list
){
    // lig_dge & recep_dge: genes X beads
    // res_mtx: bead X interactions mtx
    arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

    // for (i in ncol){ for (j in nrow) {...}}
    for (uword i = 0; i < res_mtx.n_cols; i++){
        for (uword j = 0; j < res_mtx.n_rows; j++){
            // All indices in Cpp starts with 0 --> all values -1
            // mean(x, dim): each column (dim = 0), or each row (dim = 1)
            // arma::mat lig_dge = dge_lig.submat(nb_list(span(nb_index(j), nb_index(j+1)-1)), 
            //                                 lig_list(span(lig_index(i), lig_index(i+1)-1)));
            arma::mat lig_dge = dge_lig.submat( lig_list(span(lig_index(i), lig_index(i+1)-1)), 
                                            nb_list(span(nb_index(j), nb_index(j+1)-1)) );
            // arma::vec lig_vec = mean(lig_dge, 1); // n_ligs X 1
            double lig_val = arma::accu(mean(lig_dge, 1)); // n_ligs X 1

            arma::mat lig_velo = dge_velo.submat( lig_list(span(lig_index(i), lig_index(i+1)-1)), 
                                            nb_list(span(nb_index(j), nb_index(j+1)-1)) );
            // arma::vec lig_velo_vec = mean(lig_velo, 1); // n_ligs X 1
            double lig_velo_val = arma::accu(mean(lig_velo, 1)); // n_ligs X 1

            // arma::mat recep_dge = dge_recep.submat(nb_list(span(nb_index(j), nb_index(j+1)-1)), 
            //                                 recep_list(span(recep_index(i), recep_index(i+1)-1)));
            arma::mat recep_dge = dge_recep.submat( recep_list(span(recep_index(i), recep_index(i+1)-1)), 
                                            nb_list(span(nb_index(j), nb_index(j+1)-1)) );
            // arma::rowvec recep_vec = mean(recep_dge, 1).t(); // 1 X n_receps
            double recep_val = arma::accu(mean(recep_dge, 1)); // 1 X n_receps


            arma::mat recep_velo = dge_velo.submat( recep_list(span(recep_index(i), recep_index(i+1)-1)), 
                                            nb_list(span(nb_index(j), nb_index(j+1)-1)) );
            // arma::rowvec recep_velo_vec = mean(recep_velo, 1).t(); // 1 X n_receps
            double recep_velo_val = arma::accu(mean(recep_velo, 1)); // 1 X n_receps

            res_mtx(j, i) = lig_val * recep_velo_val + recep_val * lig_velo_val;
            // res_mtx(j, i) = arma::accu(lig_vec * recep_vec);
            // lr.res = recep.exp * lig.velo + lig.exp * recep.velo
        }
    }

    return res_mtx;
}



// // [[Rcpp::export]]
// arma::mat inferScoreLR_cpp_old(
//     const arma::mat& dge_lig,
//     const arma::mat& dge_recep,
//     const arma::uvec& lig_index,
//     const arma::uvec& lig_list,
//     const arma::uvec& recep_index,
//     const arma::uvec& recep_list
// ){
//     // lig_dge & recep_dge: genes X beads
//     // res_mtx: bead X interactions mtx
//     arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

//     for (uword i = 0; i < res_mtx.n_cols; i++){ // for each ligand

//         for (uword j = 0; j < res_mtx.n_rows; j++){ // for each bead
//             // All indices in Cpp starts with 0 --> all values - 1

//             arma::vec lig_vec = dge_lig.col(j);
//             lig_vec = lig_vec.elem(lig_list(span(lig_index(i), lig_index(i+1)-1))); // n_ligs X 1
//             double lig_val = arma::accu(lig_vec); // n_ligs X 1

//             arma::vec recep_vec = dge_recep.col(j);
//             recep_vec = recep_vec.elem(recep_list(span(recep_index(i), recep_index(i+1)-1))); // n_ligs X 1
//             double recep_val = arma::accu(recep_vec); // 1 X n_receps

//             // res_mtx(j, i) = arma::accu(lig_vec * recep_vec.t());
//             res_mtx(j, i) = lig_val*recep_val;
//         }
        
//         // cout << i << endl;
//     }

//     return res_mtx;
// }



// // [[Rcpp::export]]
// arma::mat inferVeloLR_cpp_old(
//     const arma::mat& dge_lig,
//     const arma::mat& dge_recep,
//     const arma::mat& dge_lig_velo,
//     const arma::mat& dge_recep_velo,
//     const arma::uvec& lig_index,
//     const arma::uvec& lig_list,
//     const arma::uvec& recep_index,
//     const arma::uvec& recep_list
// ){
//     // lig_dge & recep_dge: genes X beads
//     // res_mtx: bead X interactions mtx
//     arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

//     for (uword i = 0; i < res_mtx.n_cols; i++){ // for each ligand

//         for (uword j = 0; j < res_mtx.n_rows; j++){ // for each bead
//             // All indices in Cpp starts with 0 --> all values - 1

//             arma::vec lig_vec = dge_lig.col(j);
//             lig_vec = lig_vec.elem(lig_list(span(lig_index(i), lig_index(i+1)-1))); // n_ligs X 1
//             double lig_val = arma::accu(lig_vec); // n_ligs X 1

//             // same for dge_velo
//             arma::vec lig_velo = dge_lig_velo.col(j);
//             lig_velo = lig_velo.elem(lig_list(span(lig_index(i), lig_index(i+1)-1))); // n_ligs X 1
//             double lig_velo_val = arma::accu(lig_velo); // n_ligs X 1

//             arma::vec recep_vec = dge_recep.col(j);
//             recep_vec = recep_vec.elem(recep_list(span(recep_index(i), recep_index(i+1)-1))); // n_ligs X 1
//             double recep_val = arma::accu(recep_vec); // 1 X n_receps

//             arma::vec recep_velo = dge_recep_velo.col(j);
//             recep_velo = recep_velo.elem(recep_list(span(recep_index(i), recep_index(i+1)-1))); // n_ligs X 1
//             double recep_velo_val = arma::accu(recep_velo); // 1 X n_receps

//             // res_mtx(j, i) = arma::accu(lig_vec * recep_vec.t());
//             res_mtx(j, i) = lig_val * recep_velo_val + recep_val * lig_velo_val;
//         }
        
//         // cout << i << endl;
//     }

//     return res_mtx;
// }




// // [[Rcpp::export]]
// arma::mat inferScoreLR_sp_cpp(
//     const arma::sp_mat& dge_lig,
//     const arma::sp_mat& dge_recep,
//     const arma::uvec& lig_index,
//     const arma::uvec& lig_list,
//     const arma::uvec& recep_index,
//     const arma::uvec& recep_list
// ){
//     // lig_dge & recep_dge: genes X beads
//     // res_mtx: bead X interactions mtx
//     arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

//     // !!!! this function is not runnable since sp_mat only supports .cols() but not .rows()
//     // !!!! need to modify the function if using sp_mat in the future

//     for (uword i = 0; i < res_mtx.n_cols; i++){ // for each intr
//         arma::sp_mat lig_dge = dge_lig.rows(lig_list(span(lig_index(i), lig_index(i+1)-1)));
//         // column sum of lig_dge
//         arma::vec lig_vec = arma::sum(lig_dge, 0); // 1 X n_beads

//         arma::sp_mat recep_dge = dge_recep.rows(recep_list(span(recep_index(i), recep_index(i+1)-1)));
//         // column sum of recep_dge
//         arma::vec recep_vec = arma::sum(recep_dge, 0); // 1 X n_beads

//         res_mtx.col(i) = (lig_vec % recep_vec).t();
    
//         // cout << i << endl;
//     }

//     return res_mtx;
// }


// // [[Rcpp::export]]
// arma::mat graphMeanLR_cpp(
//     const arma::vec& alpha_list,
//     const arma::mat& dge_lig,
//     const arma::uvec& lig_index,
//     const arma::uvec& lig_list,
//     const arma::uvec& nb_index,
//     const arma::uvec& nb_list
// ){
//     // lig_dge & recep_dge: genes X beads
//     // res_mtx: bead X interactions mtx
//     arma::mat res_mtx = arma::zeros<arma::mat>(dge_lig.n_cols, lig_index.size()-1);

//     // for (i in ncol){ for (j in nrow) {...}}
//     for (uword i = 0; i < res_mtx.n_cols; i++){
//         for (uword j = 0; j < res_mtx.n_rows; j++){
//             // All indices in Cpp starts with 0 --> all values -1
//             // mean(x, dim): each column (dim = 0), or each row (dim = 1)
//             arma::mat lig_dge = dge_lig.submat( lig_list(span(lig_index(i), lig_index(i+1)-1)), 
//                                             nb_list(span(nb_index(j), nb_index(j+1)-1)) );
//             arma::vec lig_vec = mean(lig_dge, 1); // n_ligs X 1

//             // the length of a vec
//             double lig_list_len = lig_list(span(lig_index(i), lig_index(i+1)-1)).n_elem;

//             // count the number of elements in the column which are greater than 1
//             double lig_val = arma::accu(lig_vec > alpha_list(j))/lig_list_len;
            
//             // cout << lig_vec << endl;
//             // cout << arma::accu(lig_vec > alpha_list(j)) <<endl;
//             // cout << lig_list_len << endl;
//             // cout << lig_val << endl;
            
//             res_mtx(j, i) = lig_val;
//         }
//     }

//     return res_mtx;
// }



