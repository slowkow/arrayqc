# arrayqc

An R script for quality control of microarray data, including:

-   Pseudo-images of all arrays:

    -   Raw values and cube-root of values.

    -   Probe-level model fit estimates: residuals, signed residuals, and
        weights.

-   Boxplots of raw data and RMA normalized data.

-   Boxplots of NUSE and RLE:

    -   Normalized Unscaled Standard Errors (should be centered around 1).

    -   Relative Log-Expression (should be centered around 0).

-   MA plots for all arrays:

    -   M = difference in log expression between the array and the median
        across arrays.

    -   A = average log expression between the array and the median across
        arrays.

## Installation

Install the dependencies:

```r
source("http://bioconductor.org/biocLite.R")

biocLite(
  pkgs = c("getopt", "oligo", "stringr", "reshape2", "data.table", "biomaRt"),
  suppressUpdates = TRUE,
  suppressAutoUpdate = TRUE
)
```

You also need [Imagemagick]. On OSX, use [Homebrew] to install it:

```bash
# Mac OSX
brew install imagemagick
```

## Usage

I use data from [GSE58203] to illustrate the usage and results. Please see the
figures in [figures/][figures].

```bash
opt=(
  --cel_dir data-raw
  --fig_dir figures
  --array_name affy_hg_u133_plus_2
  --mart_file data-raw/mart.csv.gz
  --out_file data/GSE58203.RData 
)

Rscript R/arrayqc.R ${opt[*]} &> arrayqc.log
```

## Metrics

Below, I summarize the included metrics and briefly explain how to interpret
them.

For more detailed information, please see:

[Bioinformatics and Computational Biology Solutions Using R and
Bioconductor][book1]

### Array Images

Large hybridization artifacts affecting thousands of probes on an array are
immediately visible upon visual inspection. It is normal to see small
blotches, lines, or gradients on the array. Taking the cube root of the raw
intensity eases identification of these artifacts. However, artifacts are most
visible when viewing the residuals and weights of a probe-level model fit.

### Probe-level Model Fit

The Robust Multi-array Average (RMA) procedure assumes a linear model for
background-adjusted normalized probe-level data:

![Probe-level Model Fit](http://mathurl.com/qznc3u6.png)

It is useful to examine the residuals and weights assigned to probes after
fitting this probe-level model (PLM).

### NUSE and RLE

The Normalized Unscaled Standard Error (NUSE) plot shows standard errors from
the PLM fit, scaled so that the median standard error across arrays is 1 for
each gene. Low quality arrays have median standard error much greater than 1,
or standard error that is more variable than the other arrays.

The Relative Log-Expression plot shows the difference between a gene's
estimated expression and the median for that gene. It assumes that most genes
are not differentially expressed, so the boxes should be centered at 0 and
have low variability.

![NUSE and RLE](http://mathurl.com/pvpduyp.png)

### MA

The MA-plot is used to compare two arrays, or to compare one array against the
summary of a group of arrays. With the logarithmic expression values, this
plot shows each gene's log fold change on the y-axis and average log intensity
on the x-axis. If most genes are not differentially expressed, then the [Loess
curve][loess] fit to this plot should follow the horizontal line at M = 0.

![MA Plot](http://mathurl.com/pafd9th.png)

## Related Work

-   [arrayQualityMetrics]

    -   This package failed to install on my system due to a failed
        installation of the [Cairo] package.

-   [arrayQuality]

    -   I have not tried this package.

[GSE58203]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58203
[book1]: http://bioconductor.org/help/publications/books/bioinformatics-and-computational-biology-solutions/
[figures]: https://github.com/slowkow/arrayqc/tree/master/figures
[loess]: https://en.wikipedia.org/wiki/Local_regression
[arrayQualityMetrics]: http://master.bioconductor.org/packages/release/bioc/html/arrayQualityMetrics.html
[arrayQuality]: http://master.bioconductor.org/packages/release/bioc/html/arrayQuality.html
[Cairo]: http://cran.r-project.org/web/packages/Cairo/index.html
[Imagemagick]: http://www.imagemagick.org/
[Homebrew]: http://brew.sh/
