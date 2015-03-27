# Pseudo-Images of all Arrays

## Raw Values

![raw values](https://github.com/slowkow/arrayqc/blob/master/figures/rawData-images.jpg)

## Cube-root of Raw Values

Taking the cube-root enhances ability to notice linear trends in color
intensity. It is already clear that the first batch of 30 arrays shows
lower intensity (deeper blue) than the second batch of 30 arrays.

![cuberoot of raw values](https://github.com/slowkow/arrayqc/blob/master/figures/rawData-images-cuberoot.jpg)

## PLM Fit Residuals

![residuals](https://github.com/slowkow/arrayqc/blob/master/figures/plmFit-residuals.jpg)

## PLM Fit Signed Residuals

![signed residuals](https://github.com/slowkow/arrayqc/blob/master/figures/plmFit-sign_residuals.jpg)

## PLM Fit Weights

The weights often show the clearest picture of hybridization artifacts.
I always expect to see some blotches, snakes, or gradients. None of the
artifacts below are bad enough to discard an array from further analysis.

![weights](https://github.com/slowkow/arrayqc/blob/master/figures/plmFit-weights.jpg)

# Boxplots of Raw and RMA Normalized Data

## Raw Data

The samples were clearly done in two batches.

![raw data](https://github.com/slowkow/arrayqc/blob/master/figures/rawData-boxplot.png)

## RMA Normalized Data

RMA uses [quantile normalization][1], which forces all distributions to be
identical after normalization.

![RMA data](https://github.com/slowkow/arrayqc/blob/master/figures/rmaData-boxplot.png)

# Boxplots of NUSE and RLE

## NUSE

After RMA normalization, it is still abundantly clear that the samples were
done in two batches. The batch-effect is not fully corrected by RMA. All
distributions should be centered around 1. Instead, the first batch is
elevated, indicating greater variance in expression values, and the second
batch is depressed.

![NUSE](https://github.com/slowkow/arrayqc/blob/master/figures/plmFit-NUSE.png)

## RLE

![RLE](https://github.com/slowkow/arrayqc/blob/master/figures/plmFit-RLE.png)

# MA Plots

It is typical to assume that most genes in an experiment are not
differentially expressed. If this assumption holds, then the red Loess curves
fit to these scatter plots should be centered around 0. Ideally, each of these
plots should share the same x and y axes to facilitate comparison.

![MA](https://github.com/slowkow/arrayqc/blob/master/figures/rmaData-ma.jpg)

[1]: https://en.wikipedia.org/wiki/Quantile_normalization
