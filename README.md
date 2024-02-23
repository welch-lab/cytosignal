# CytoSignal & VelcoCytoSignal

CytoSignal & VeloCytoSignal is a tool for detecting static and dynamic cell-cell signaling interactions at single-cell resolution from spatial transcriptomic data.

## Installation

The package is developed and tested under R>=4.2.0. Users can install R following the [instruction provided on CRAN](https://cran.r-project.org/). [RStudio](https://posit.co/downloads/) is a recommended IDE for working with R projects. 

To install CytoSimplex in R, run the following command in an R console:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("welch-lab/CytoSimplex")
```

We are currently working on releasing a CRAN version. Stay tuned!

## Tutorial

For usage examples and guided walkthroughs, check the vignettes directory of the repo.

* [Infer spatially resolved cell-cell communication signaling at cellular resolution]()
* [Infer spatially resolved temporal dynamics of cell-cell communication at cellular resolution]()

## Citation

If you used CytoSignal in your work, please cite the following work:

>TBA
