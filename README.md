[![Travis-CI Build Status](https://travis-ci.org/kevinrue/SeuratConverter.svg?branch=master)](https://travis-ci.org/kevinrue/SeuratConverter)
[![codecov Status](https://codecov.io/gh/kevinrue/SeuratConverter/branch/master/graph/badge.svg?token=By8PFPNXX7)](https://codecov.io/gh/kevinrue/SeuratConverter)

# `SeuratConverter`

## Overview

This package aims at making
[`seurat`](https://CRAN.R-project.org/package=Seurat) objects
accessible for downstream
processing using common
[_Bioconductor_](https://www.bioconductor.org) methods and classes.
In particular, the most readily equivalent _Bioconductor_ class apt
to store the various components of a `seurat` object is the
`SingleCellExperiment` class
([`SingleCellExperiment`](http://bioconductor.org/packages/SingleCellExperiment/)
package).
Indeed, the `SingleCellExperiment` class contains slots to store:
* any number of _assays_, that encompasses the `@raw.data`, `@data`, and
    `@scaled.data` slots of a `seurat` object.
* sample meta-data, that corresponds to the `@meta.data` slot of a `seurat`
    object.
* feature meta-data, that encompasses the presence of individual genes in the
    `@var.genes` slot of a `seurat` object
* any number of reduced dimensions sets of coordinates, that are usually
    stored under the `@dr` slot of a `seurat` object

## Installation

The `r Githubpkg("kevinrue/SeuratConverter")` package is currently hosted
on GitHub, and may be installed as follows:

```{r install, eval=FALSE}
devtools::install_github("kevinrue/SeuratConverter", build_vignettes = TRUE)
```
