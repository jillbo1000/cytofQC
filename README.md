# cytofQC

## A package for cleaning CyTOF data

CyTOF data includes several event types that should be removed from the data prior to analysis. The package cytofQC uses learning algorithms to identify and label each event in a CyTOF dataset as a bead, doublet, debris, dead cell, GDP zero (meaning at least one of the QC variables was 0), or a viable cell. 

Authors: Jill Lundell and Kelly Street

Maintainer: Jill Lundell <jflundell@gmail.com>

### Installation

To install this package, start R and enter: 

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("cytofQC")

```

### Documenation

To view the documentation for the version of this package installed on your system, start R and enter: 

```
browseVignettes("cytofQC")
```
