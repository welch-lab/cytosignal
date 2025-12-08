#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <cli/progress.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
using namespace arma;

// Simply found that doing the matrix multiplication in C++ is faster than in R
// Use L/Rmap to collapse complex expression into single component
// impL - Imputed Ligand value, gene x spot
// Lmap - Binary ligand gene mapping per interaction, gene x interaction
// impR - Imputed Receptor value, gene x spot
// Rmap - Binary receptor gene mapping per interaction, gene x interaction
// returns spot x interaction matrix
// [[Rcpp::export]]
arma::sp_mat multiply_lr_cpp(
        const arma::sp_mat& impL,
        const arma::sp_mat& Lmap,
        const arma::sp_mat& impR,
        const arma::sp_mat& Rmap
) {
    arma::sp_mat L = impL.t() * Lmap;
    arma::sp_mat R = impR.t() * Rmap;
    arma::sp_mat result = L % R;
    return result;
}

// Replicate R findInterval(left_open = TRUE)
// [[Rcpp::export]]
arma::uword findInterval_leftOpen_cpp(double x, const arma::vec& breaks) {
    auto it = std::lower_bound(breaks.begin(), breaks.end(), x);
    // # of elements < x  (index in [0, n])
    return static_cast<int>(it - breaks.begin());
}

// Inplace update the sampled averaged null L/R score per permutation per interaction
// raw - Transposed raw count matrix, [spot x valid-gene]
// map - Binary ligand/receptor gene mapping per interaction, [valid-gene x interaction]
// nullgraph - Transposed null neighbor graph, [spot x spot]
// nullLibSize - Library size per spot under null graph [spot]
// dtAvgNullGraph - Transposed from to_mean(contGraph)[, permIdx], [sampled-spot x spot]
// nullSub - Pre-allocated spot x 1 matrix to hold the collapsed null ligand
//            or receptor value per permutation per interaction
// nullSampled - Pre-allocated sampled-averaged sampleSize x 1 matrix to hold the final null
//                ligand or receptor value per permutation per interaction
// i - Interaction index
void collapse_null_value_perm_intr(
    const arma::sp_mat& raw,
    const arma::sp_mat& map,
    const arma::sp_mat& nullgraph,
    const arma::vec& nullLibSize,
    const arma::sp_mat& dtAvgNullGraph,
    arma::vec& nullSub,
    arma::vec& nullSampled,
    arma::uword i
) {
    nullSub.zeros();
    for (arma::sp_mat::const_col_iterator it = map.begin_col(i);
         it != map.end_col(i); ++it) {
        arma::uword gene_idx = it.row();
        arma::sp_mat gene_expr = raw.col(gene_idx);
        // Rcpp::Rcout << "nullgraph density: " << double(nullgraph.n_nonzero)/double(nullgraph.n_rows*nullgraph.n_cols) << "\n";
        // Rcpp::Rcout << "gene_expr density: " << double(gene_expr.n_nonzero)/double(gene_expr.n_rows*gene_expr.n_cols) << "\n";
        gene_expr = nullgraph * gene_expr;
        for (arma::sp_mat::iterator gene_it = gene_expr.begin();
             gene_it != gene_expr.end(); ++gene_it) {
            (*gene_it) /= nullLibSize[gene_it.row()];
            (*gene_it) *= 1.0e4;
            (*gene_it) = std::log1p(*gene_it);
        }
        nullSub += gene_expr;
    }
    nullSampled = dtAvgNullGraph * nullSub;
}

// Permutation test for LR score significance
// raw - Transposed raw count matrix, [spot x valid-gene]
//       It's important to subset columns instead of subset rows for the sake of
//       sparse matrix efficiency
// libSize - Library size per spot [spot]
// intrType - Interaction type per interaction, 0-contact, 1-diffusion [intr]
// "NGV" -- Null graph vector; "diff"/"cont" for diffusion/contact-dependent interactions
// lig_NGV - null graphs for imputing ligand gene expression [spot x spot x nPerm]
// recep_NGV - null graphs for imputing receptor gene expression [spot x spot x nPerm]
// dtAvg_NGV - null graphs for sampling averaged spot expression [sampled-spot x spot x nPerm]
// smooth_NGV - null graphs for smoothing null LR score [spot x spot x nPerm]
// Lmap - Binary ligand gene mapping per interaction, [valid-gene x interaction]
// Rmap - Binary receptor gene mapping per interaction, [valid-gene x interaction]
// lrscore - Observed LR score matrix, [spot x interaction]
arma::sp_mat permute_test_cpp(
    const arma::sp_mat& raw,
    const arma::vec& libSize,
    const arma::uvec& intrType,
    const std::vector<arma::sp_mat>& diff_lig_NGV,
    const std::vector<arma::sp_mat>& diff_recep_NGV,
    const std::vector<arma::sp_mat>& diff_dtAvg_NGV,
    const std::vector<arma::sp_mat>& diff_smooth_NGV,
    const std::vector<arma::sp_mat>& cont_lig_NGV,
    const std::vector<arma::sp_mat>& cont_recep_NGV,
    const std::vector<arma::sp_mat>& cont_dtAvg_NGV,
    const std::vector<arma::sp_mat>& cont_smooth_NGV,
    const arma::sp_mat& Lmap,
    const arma::sp_mat& Rmap,
    const arma::sp_mat& lrscore
) {
    arma::uword nSpot = raw.n_rows;
    arma::uword nIntr = Lmap.n_cols;
    arma::uword nPerm = diff_lig_NGV.size();
    arma::uword sampleSize = diff_dtAvg_NGV[0].n_rows;
    // Copy lrscore to ecdf. P-value will be obtained via $1 - ecdf$.
    // Zero-LRscore will obviously have the empirical proportion of zero, so no
    // need to calculate for those entries, meaning that we can directly use the
    // non-zero entries of LRscore.
    arma::sp_mat ecdf = lrscore;

    // NLSV - Null library size vector
    // Precalculate these Null-imputed library sizes since they are reused
    // by each interaction according to the original implementation.
    std::vector<arma::mat> diff_lig_NLSV;
    diff_lig_NLSV.reserve(nPerm);
    std::vector<arma::mat> diff_recep_NLSV;
    diff_recep_NLSV.reserve(nPerm);
    std::vector<arma::mat> cont_lig_NLSV;
    cont_lig_NLSV.reserve(nPerm);
    std::vector<arma::mat> cont_recep_NLSV;
    cont_recep_NLSV.reserve(nPerm);

    for (arma::uword i = 0; i < nPerm; i++) {
        Rcpp::checkUserInterrupt();
        const arma::sp_mat& diff_lig_NG = diff_lig_NGV[i];
        arma::mat diff_lig_NLS = diff_lig_NG * libSize;
        diff_lig_NLSV.push_back(std::move(diff_lig_NLS));

        const arma::sp_mat& diff_recep_NG = diff_recep_NGV[i];
        arma::mat diff_recep_NLS = diff_recep_NG * libSize;
        diff_recep_NLSV.push_back(std::move(diff_recep_NLS));

        const arma::sp_mat& cont_lig_NG = cont_lig_NGV[i];
        arma::mat cont_lig_NLS = cont_lig_NG * libSize;
        cont_lig_NLSV.push_back(std::move(cont_lig_NLS));

        const arma::sp_mat& cont_recep_NG = cont_recep_NGV[i];
        arma::mat cont_recep_NLS = cont_recep_NG * libSize;
        cont_recep_NLSV.push_back(std::move(cont_recep_NLS));
    }

    // Now start to calculate null score for each interaction. Doing so because
    // holding the null score (assuming 100,000 sample size) of all interactions
    // will be memory intensive.

    // *nullSub is pre-allocated for holding collapsed normalized null-imputed
    // L/R value per permutation per interaction
    // *nullSampled is pre-allocated for holding sampled-dt-averaged *nullSub
    // values per permutation per interaction
    arma::vec LnullSub(nSpot);
    arma::vec LnullSampled(sampleSize);

    arma::vec RnullSub(nSpot);
    arma::vec RnullSampled(sampleSize);

    // nullSub is the final null LR score per permutation per interaction
    arma::vec nullSub(sampleSize);
    // nullScore holds all null LR scores for all permutations per interaction
    arma::vec nullScore(sampleSize*nPerm);

    Rcpp::List cli_config = Rcpp::List::create(
        Rcpp::Named("name") = "Calculating ECDF",
        Rcpp::Named("clear") = true
    );
    SEXP bar = PROTECT(cli_progress_bar(nIntr, cli_config));

    for (arma::uword i = 0; i < nIntr; i++) {
        const std::vector<arma::sp_mat>* lig_NGV = (intrType[i] == 0) ? &cont_lig_NGV : &diff_lig_NGV;
        const std::vector<arma::mat>* lig_NLSV = (intrType[i] == 0) ? &cont_lig_NLSV : &diff_lig_NLSV;
        const std::vector<arma::sp_mat>* recep_NGV = (intrType[i] == 0) ? &cont_recep_NGV : &diff_recep_NGV;
        const std::vector<arma::mat>* recep_NLSV = (intrType[i] == 0) ? &cont_recep_NLSV : &diff_recep_NLSV;
        const std::vector<arma::sp_mat>* dtAvg_NGV = (intrType[i] == 0) ? &cont_dtAvg_NGV : &diff_dtAvg_NGV;
        const std::vector<arma::sp_mat>* smooth_NGV = (intrType[i] == 0) ? &cont_smooth_NGV : &diff_smooth_NGV;

        for (arma::uword j = 0; j < nPerm; j++) {
            const arma::sp_mat& lig_NG = (*lig_NGV)[j];
            const arma::mat& lig_NLS = (*lig_NLSV)[j];
            const arma::sp_mat& recep_NG = (*recep_NGV)[j];
            const arma::mat& recep_NLS = (*recep_NLSV)[j];
            const arma::sp_mat& dtAvg_NG = (*dtAvg_NGV)[j];
            const arma::sp_mat& smooth_NG = (*smooth_NGV)[j];

            collapse_null_value_perm_intr(raw, Lmap, lig_NG, lig_NLS, dtAvg_NG,
                                          LnullSub, LnullSampled, i);
            collapse_null_value_perm_intr(raw, Rmap, recep_NG, recep_NLS, dtAvg_NG,
                                          RnullSub, RnullSampled, i);
            nullSub = LnullSampled % RnullSampled;
            nullSub =  smooth_NG * nullSub;
            nullScore.subvec(j*sampleSize, (j+1)*sampleSize-1) = nullSub;
        }

        // then calculate ecdf
        arma::vec nullScoreSorted = arma::sort(nullScore, "ascend", 0);
        for (arma::sp_mat::col_iterator it = ecdf.begin_col(i); it != ecdf.end_col(i); ++it) {
            arma::uword ecdf = findInterval_leftOpen_cpp(*it, nullScoreSorted);
            *it = double(ecdf)/double(nPerm*sampleSize);
        }

        if (CLI_SHOULD_TICK) cli_progress_set(bar, i);
    }
    cli_progress_done(bar);
    UNPROTECT(1);
    return ecdf;
}


// Rcpp wrapper for perm_test_cpp, just to convert Rcpp::List to std::vector
// Mainly for allowing portability of the code.
// [[Rcpp::export]]
arma::sp_mat perm_test_Rcpp(
    const arma::sp_mat& raw,
    const arma::vec& libSize,
    const arma::uvec& intrType,
    Rcpp::List diff_lig_NGL,
    Rcpp::List diff_recep_NGL,
    Rcpp::List diff_dtAvg_NGL,
    Rcpp::List diff_smooth_NGL,
    Rcpp::List cont_lig_NGL,
    Rcpp::List cont_recep_NGL,
    Rcpp::List cont_dtAvg_NGL,
    Rcpp::List cont_smooth_NGL,
    const arma::sp_mat& Lmap,
    const arma::sp_mat& Rmap,
    const arma::sp_mat& lrscore
) {
    arma::uword nPerm = diff_lig_NGL.size();
    std::vector<arma::sp_mat> diff_lig_NGV;
    diff_lig_NGV.reserve(nPerm);
    std::vector<arma::sp_mat> diff_recep_NGV;
    diff_recep_NGV.reserve(nPerm);
    std::vector<arma::sp_mat> diff_dtAvg_NGV;
    diff_dtAvg_NGV.reserve(nPerm);
    std::vector<arma::sp_mat> diff_smooth_NGV;
    diff_smooth_NGV.reserve(nPerm);
    std::vector<arma::sp_mat> cont_lig_NGV;
    cont_lig_NGV.reserve(nPerm);
    std::vector<arma::sp_mat> cont_recep_NGV;
    cont_recep_NGV.reserve(nPerm);
    std::vector<arma::sp_mat> cont_dtAvg_NGV;
    cont_dtAvg_NGV.reserve(nPerm);
    std::vector<arma::sp_mat> cont_smooth_NGV;
    cont_smooth_NGV.reserve(nPerm);

    for (arma::uword i = 0; i < nPerm; i++) {
        diff_lig_NGV.push_back(std::move(diff_lig_NGL[i]));
        diff_recep_NGV.push_back(std::move(diff_recep_NGL[i]));
        diff_dtAvg_NGV.push_back(std::move(diff_dtAvg_NGL[i]));
        diff_smooth_NGV.push_back(std::move(diff_smooth_NGL[i]));
        cont_lig_NGV.push_back(std::move(cont_lig_NGL[i]));
        cont_recep_NGV.push_back(std::move(cont_recep_NGL[i]));
        cont_dtAvg_NGV.push_back(std::move(cont_dtAvg_NGL[i]));
        cont_smooth_NGV.push_back(std::move(cont_smooth_NGL[i]));
    }

    arma::sp_mat ecdf = permute_test_cpp(
        raw,
        libSize,
        intrType,
        diff_lig_NGV,
        diff_recep_NGV,
        diff_dtAvg_NGV,
        diff_smooth_NGV,
        cont_lig_NGV,
        cont_recep_NGV,
        cont_dtAvg_NGV,
        cont_smooth_NGV,
        Lmap,
        Rmap,
        lrscore
    );
    return ecdf;
}
