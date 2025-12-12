#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <cli/progress.h>
#include "progress.hpp"
#include "eta_progress_bar.hpp"

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

    ETAProgressBar pb;
    Progress p(nIntr, true, pb);
    for (arma::uword i = 0; i < nIntr; i++) {
        if (Progress::check_abort()) {
            return ecdf;
        }
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

        p.increment();
    }
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


// ecdf - Empirical CDF matrix, output of previous functions [spot x interaction]
// intrType - Interaction type per interaction, 0-contact, 1-diffusion [intr]
// diffGraph - Diffusion-based neighbor graph [spot x spot]
// contGraph - Contact-based neighbor graph [spot x spot]
// Returns adjusted p-value matrix [spot x interaction]
// [[Rcpp::export]]
arma::mat spatialGraphFDR_cpp(
        const arma::mat& pval,
        const arma::uvec& intrType,
        const arma::sp_mat& contGraph
) {
    arma::uword nIntr = pval.n_cols;
    arma::uword nSpot = pval.n_rows;
    arma::mat fdr(pval.n_rows, nIntr, arma::fill::ones);

    // Initialize the weighting for the weighted BH adjustment from Cydar
    // Basically divide 1 by the maximum neighbor weight per spot
    arma::vec weights(nSpot, arma::fill::ones);
    double max_weight;
    for (arma::uword i = 0; i < nSpot; i++) {
        max_weight = 0.0;
        for (arma::sp_mat::const_col_iterator it = contGraph.begin_col(i);
             it != contGraph.end_col(i); ++it) {
            if (it.row() == i) continue;
            if (*it > max_weight) {
                max_weight = *it;
            }
            if (max_weight > 0.0) weights[i] = 1.0 / max_weight;
        }
    }
    double weight_sum = arma::sum(weights);

    arma::vec pval_vec(nSpot);
    arma::vec psorted(nSpot);
    arma::vec wsub(nSpot);

    for (arma::uword i = 0; i < nIntr; i++) {
        Rcpp::checkUserInterrupt();
        pval_vec = pval.col(i);

        arma::uvec naIdx = arma::find_nonfinite(pval_vec);
        arma::uvec o = arma::sort_index(pval_vec, "ascend");
        psorted = pval_vec(o);
        wsub = weights(o);
        arma::vec padj(nSpot);
        padj.fill(-1.0);
        arma::vec weight_cumsum = arma::cumsum(wsub);
        arma::vec padj_sorted = weight_sum * psorted / weight_cumsum;
        padj_sorted = arma::reverse(padj_sorted);
        for (arma::uword j = 1; j < nSpot; j++) {
            double last = padj_sorted[j-1];
            double now = padj_sorted[j];
            if (now > last) {
                padj_sorted[j] = last;
            }
        }
        padj_sorted = arma::reverse(padj_sorted);
        for (arma::uword j = 0; j < nSpot; j++) {
            if (padj_sorted[j] > 1.0) {
                padj_sorted[j] = 1.0;
            }
        }
        padj(o) = padj_sorted;
        padj(naIdx).fill(NA_REAL);
        fdr.col(i) = padj;
    }
    return fdr;
}




// raw - Gene expression matrix [(valid gene) x spot]
// libSize - Library size per spot [spot]
// clusterInt - 0-based integer cluster assignment per spot [spot]
// nCluster - Number of clusters
// intrType - Interaction type per interaction, 0-contact, 1-diffusion [intr]
// diff_lig_graph - Diffusion-based ligand neighbor graph [(sender)spot x (receiver)spot]
// cont_lig_graph - Contact-based ligand neighbor graph [(sender)spot x (receiver)spot]
// recep_graph - Receptor neighbor graph [(sender)spot x (receiver)spot]
// dtAvg_graph - Delauney triangulation neighbor-based averaging graph [(sender)spot x (receiver)spot]
// Lmap - Binary ligand gene mapping per interaction, [valid-gene x interaction]
// Rmap - Binary receptor gene mapping per interaction, [valid-gene x interaction]
// Returns cluster-wise LR score cube [(sender)cluster x (receiver)cluster x interaction]
// [[Rcpp::export]]
arma::cube clusterWiseLRscore_cpp(
    const arma::sp_mat& raw,
    const arma::uvec& libSize,
    const arma::uvec& clusterInt,
    const arma::uword& nCluster,
    const arma::uvec& intrType,
    const arma::sp_mat& diff_lig_graph,
    const arma::sp_mat& cont_lig_graph,
    const arma::sp_mat& recep_graph,
    const arma::sp_mat& dtAvg_graph,
    const arma::sp_mat& Lmap,
    const arma::sp_mat& Rmap
) {
    arma::cube clusterWiseLRScore(nCluster, nCluster, Lmap.n_cols, arma::fill::zeros);
    arma::uword nSpot = raw.n_cols;
    arma::uword nIntr = Lmap.n_cols;
    arma::vec clusterSizes(nCluster, arma::fill::zeros);
    for (arma::uword i = 0; i < nSpot; i++) {
        clusterSizes[clusterInt[i]] += 1.0;
    }

    // The basic idea is to do the LR imputation for each receiver spot, per
    // cluster. So below Lgene(c,r) holds the weighted-sum of ligand genes from
    // only neighboring c-type cells received by r.
    arma::mat Lgene(nCluster, nSpot);
    arma::mat L(nCluster, nSpot, arma::fill::zeros);
    // Whereas receptor amount is on the receiver spot itself, which already
    // belongs to a certain cluster. We later aggregate columns by assignment.
    arma::mat Rgene(1, nSpot);
    arma::mat R(1, nSpot, arma::fill::zeros);
    // r for receiver spot index
    // s for sender spot index
    // g for gene index in raw rows
    // i for interaction index
    // c for cluster index
    arma::uword r, s, g, i, c;
    arma::mat diff_libSizes(nCluster, nSpot);
    arma::mat cont_libSizes(nCluster, nSpot);
    for (r = 0; r < nSpot; r++) {
        for (arma::sp_mat::const_col_iterator graph_it = diff_lig_graph.begin_col(r);
             graph_it != diff_lig_graph.end_col(r); ++graph_it) {
            s = graph_it.row();
            c = clusterInt[s];
            diff_libSizes(c, r) += libSize[s] * (*graph_it);
        }
        for (arma::sp_mat::const_col_iterator graph_it = cont_lig_graph.begin_col(r);
             graph_it != cont_lig_graph.end_col(r); ++graph_it) {
            s = graph_it.row();
            c = clusterInt[s];
            cont_libSizes(c, r) += libSize[s] * (*graph_it);
        }
    }
    // Progress bar
    ETAProgressBar pb;
    Progress p(nIntr, true, pb);
    for (i = 0; i < nIntr; i++) {
        if (Progress::check_abort()) {
            return clusterWiseLRScore;
        }
        const arma::sp_mat& lig_graph = (intrType[i] == 0) ? cont_lig_graph : diff_lig_graph;
        arma::mat& libSizes = (intrType[i] == 0) ? cont_libSizes : diff_libSizes;
        L.zeros();
        for (arma::sp_mat::const_col_iterator intr_it = Lmap.begin_col(i);
             intr_it != Lmap.end_col(i); ++intr_it) {
            g = intr_it.row();
            Lgene.zeros();
            for (r = 0; r < nSpot; r++) {
                for (arma::sp_mat::const_col_iterator graph_it = lig_graph.begin_col(r);
                     graph_it != lig_graph.end_col(r); ++graph_it) {
                    s = graph_it.row();
                    c = clusterInt[s];
                    Lgene(c, r) += raw(g, s) * (*graph_it) * 1.0e4;
                }
                for (c = 0; c < nCluster; c++) {
                    if (Lgene(c, r) > 0) {
                        Lgene(c, r) /= libSizes(c, r);
                    };
                }
            }
            Lgene = arma::log1p(Lgene);
            // Apparently, this is not the right place to take this DT-average, if we want to do so
            // Lgene = Lgene * dtAvg_graph;
            L += Lgene;
        }

        R.zeros();
        for (arma::sp_mat::const_col_iterator intr_it = Rmap.begin_col(i);
             intr_it != Rmap.end_col(i); ++intr_it) {
            g = intr_it.row();
            Rgene.zeros();
            for (r = 0; r < nSpot; r++) {
                for (arma::sp_mat::const_col_iterator graph_it = recep_graph.begin_col(r);
                     graph_it != recep_graph.end_col(r); ++graph_it) {
                    s = graph_it.row();
                    Rgene(0, r) += raw(g, s) * (*graph_it) * 1.0e4;
                }
                if (Rgene(0, r) > 0) {
                    Rgene(0, r) /= libSize[r];
                };
            }
            Rgene = arma::log1p(Rgene);
            // Apparently, this is not the right place to take this DT-average, if we want to do so
            // Rgene = Rgene * dtAvg_graph;
            R += Rgene;
        }
        for (r = 0; r < nSpot; r++) {
            c = clusterInt[r];
            clusterWiseLRScore.slice(i).col(c) += (L.col(r) * R(0, r))/clusterSizes[c];
        }
        p.increment();
    }
    return(clusterWiseLRScore);
}

