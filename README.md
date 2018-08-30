# CNEr
Conserved Noncoding Elements (CNEs) Identification and Visualisation

## Installation of the stable version of `CNEr` from Bioconductor

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("CNEr")
```

## Installation of the development version of `CNEr` from github
**Prerequsite**:

  * Mac: Install "Command Line Tools" via `gcc` on terminal
  * Linux: Install a compiler and various development libraries (details vary across different flavors of Linux).
  * Windows: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

```R
devtools::install_github("ge11232002/CNEr")
```

## Vignette
Latest vignette is available at http://rpubs.com/yang2/CNEr3
