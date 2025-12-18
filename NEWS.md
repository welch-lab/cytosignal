## CytoSignal 0.99.9003 - 12/18/2025

- Simplified user-side workflow
    - Merged unnecessary operations into single functions
    - Reduced number of parameters to set
    - Renamed functions for intuitive understanding of their purposes
    - See [new vignettes for basic CytoSignal analysis](https://htmlpreview.github.io/?https://github.com/welch-lab/cytosignal/blob/refactor/doc/cytosignal2.html)
- Refactored object class
    - Reduced slots for storing unnecessary yet large-sized information, which was kept mainly for debugging purpose in the past.
    - Simplified structure to allow intuitive user navigation
    - Added useful S3/S4 methods for R-style object manipulation
- Refactored and reproducible permutation test fully in Rcpp
    - Improved speed
    - Largely reduced memory usage
    - Allowed direct access of spatialFDR values instead of storing only binary filtering results
- Allowed flexible result access with `getTopIntr()`
    - Returning database with significance metrics attached, that can be freely manipulated
    - Results are now more readable
- Added new visualization methods
    - `plotScale()` and `plotNeighbor*()` for visually checking if the scale factor and the neighborhood inference look appropriate
    - `plotSpatial*()` functions for intuitive and efficient visualization with more information
- Added new functionality
    - `inferClusterLRScore()` allows partial ligand imputation from spots belonging to each cluster in a neighborhood, and hence returns
      mean LRscore for each sender-cluster receiver-cluster pair
- Improved documentation
    - Algorithmic details intuitively explained
    - Data structure manipulation explained
