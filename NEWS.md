## CytoSignal 0.99.9002 - 12/12/2025

- Refactored object class
    - Reduced slots for storing unnecessary yet large-sized information, which was kept mainly for debugging purpose in the past.
    - Simplified structure to allow intuitive user navigation
    - Added useful S3/S4 methods for R-style object manipulation
- Refactored and reproducible permutation test fully in Rcpp
    - Improved speed
    - Largely reduced memory usage
    - Allowed direct access of p-values and spatialFDR values instead of storing only binary filtering results
- Added new visualization methods
    - `plotScale()` and `plotNeighbor*()` for visually checking if the scale factor and the neighborhood inference look appropriate
    - `plotSpatial*()` functions for intuitive and efficient visualization with more information
- Added new functionality
    - `inferClusterLRScore()` allows partial ligand imputation from spots belonging to each cluster in a neighborhood, and hence returns
      mean LRscore for each sender-cluster receiver-cluster pair
- Improved documentation
    - Algorithmic details intuitively explained
    - Data structure manipulation explained
