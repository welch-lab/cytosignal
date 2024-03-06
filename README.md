# CytoSignal & VeloCytoSignal
<img alt="ViewCount" src="https://views.whatilearened.today/views/github/welch-lab/cytosignal.svg">

`CytoSignal & VeloCytoSignal` is a tool for detecting static and dynamic cell-cell signalings at single-cell resolution from spatial transcriptomic data. Our tools are applicable to most sequencing-based and probe-based spatial transcriptomic techniques including Slide-seq, Stereo-seq, MERFISH, requiring only a cell-by-gene matrix and a cell-by-spatial-position matrix.

Check out our paper for a more complete description of the methods and analyses:

>[TBD]()

## Installation

The package is developed and tested under R>=4.2.0. Users can install R following the [instruction provided on CRAN](https://cran.r-project.org/). To install `CytoSimplex & VelcoCytoSignal` in R, run the following command in an R console:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("welch-lab/cytosignal")
```

We are currently working on releasing a CRAN version. Stay tuned!

## Tutorial

For usage examples and guided walkthroughs, check the vignettes directory of the repo.

* [Infer spatially resolved cell-cell communication signaling at cellular resolution](https://htmlpreview.github.io/?https://github.com/welch-lab/cytosignal/blob/master/doc/cytosignal.html)
* [Infer spatially resolved temporal dynamics of cell-cell communication at cellular resolution](https://htmlpreview.github.io/?https://github.com/welch-lab/cytosignal/blob/master/doc/velocytosignal.html)

## License
This project is covered under the **GNU General Public License 3.0**.
