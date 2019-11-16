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

<b>Citation</b>:

```
G. Tan, D. Polychronopoulos, B.Lenhard: 
CNEr: A toolkit for exploring extreme noncoding conservation.
PLoS Comput Biol. 2019 Aug 26;15(8):e1006940. doi: 10.1371/journal.pcbi.1006940. eCollection 2019 Aug. 
```
