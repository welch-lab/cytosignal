#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// [[Rcpp::export]]
std::string check_omp_threads() {
#ifdef _OPENMP
    int ncores = omp_get_max_threads();
    return "OpenMP is available. Max threads: " + to_string(ncores);
#else
    return "OpenMP is not available.";
#endif
}
