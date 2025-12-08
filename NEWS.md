## CytoSignal 0.99.9001 - 12/08/2025

- Refactored object class
    - Reduced slots for storing unnecessary yet large-sized information, which was kept mainly for debugging purpose in the past.
    - Simplified structure to allow intuitive user navigation
    - Added useful S3/S4 methods for R-style object manipulation
- Refactored and reproducible permutation test fully in Rcpp
    - Slightly improved speed
    - Largely reduced memory usage
    - Allow efficient access of p-values instead of storing only names of significant interactions
- Added new visualization methods
    - `plotScale()` and `plotNeighbor*()` for visually checking if the scale factor and the neighborhood inference look appropriate
- Improved documentation
    - Algorithmic details intuitively explained
    - Data structure manipulation explained
