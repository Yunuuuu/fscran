
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fast single-cell RNA-seq data analysis (fscran)

⚠️ ⚠️ ⚠️ ⚠️ **This repository is deprecated, use the
[scend](https://github.com/Yunuuuu/scend) package instead.** ⚠️ ⚠️ ⚠️ ⚠️

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of fscran from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("pak")) {
    install.packages("pak",
        repos = sprintf(
            "https://r-lib.github.io/p/pak/devel/%s/%s/%s",
            .Platform$pkgType, R.Version()$os, R.Version()$arch
        )
    )
}
pak::pkg_install("Yunuuuu/fscran")
```

## Introduction

[scran.chan](https://github.com/LTLA/scran.chan) is a package provides
methods for an end-to-end analysis of a single-cell RNA-sequencing
(scRNA-seq) analysis, starting from the count matrix and finishes with
clusters, markers, and various embeddings (i.e., t-SNE and UMAP). It’s
pretty fast and memory-efficient. The `fscran` package serves as a
bridging tool connecting the `scran.chan` package with either `Seurat`
or `SingleCellExperiment` objects.

## Available functions

| functions         | Description                                   |
|-------------------|-----------------------------------------------|
| `logNormCounts`   | Log-transformed normalized expression         |
| `runPCA`          | Principal component analysis                  |
| `logNormAndPCA`   | `logNormCounts` + `runPCA`                    |
| `runMNN`          | Fast mutual nearest neighbors correction      |
| `quickMNN`        | `logNormAndPCA` + `runMNN`                    |
| `downsample`      | Downsample cells based on their neighbors     |
| `runTSNE`         | t-stochastic neighbor embedding               |
| `runUMAP`         | uniform manifold approximation and projection |
| `clusterSNNGraph` | Graph-based clustering                        |
| `subCluster`      | Find subclusters                              |
| `scoreMarkers`    | Score marker genes                            |
| `scoreFeatureSet` | Score feature set activity for each cell      |
