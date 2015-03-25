#!/usr/bin/env Rscript
# arrayqc.R
# Kamil Slowikowski
# 
# Perform quality control analysis on CEL files.

main <- function() {
  library(getopt)
  args <- parse_args()

  suppressMessages(library(oligo))
  suppressMessages(library(stringr))
  suppressMessages(library(reshape2))
  suppressMessages(library(data.table))
  suppressMessages(library(biomaRt))
  suppressMessages(library(R.utils))

  # If you have multiple CPUs, You can speed up RMA and PLM with more threads.
  Sys.setenv(R_THREADS = 2)

  cel_files <- list.celfiles(
    args[['cel_dir']],
    listGzipped = TRUE,
    full.names = TRUE
  )
  fig_dir <- args[['fig_dir']]
  array_name <- args[['array_name']]
  out_file <- args[['out_file']]
  mart_file <- args[['mart_file']]

  dir.create(fig_dir, showWarnings = FALSE)

  n_cel_files <- length(cel_files)

  # Read microarray data -----------------------------------------------------
  rawData <- read.celfiles(cel_files)
  colnames(rawData) <- 1:ncol(rawData)

  png(file.path(fig_dir, "rawData-boxplot.png"),
      width = 14 * n_cel_files, height = 600)
  boxplot(rawData, main = "Raw Data", cex.axis = 1.2, cex.lab = 1.2)
  dev.off()

  # RMA ----------------------------------------------------------------------
  rmaData <- rma(rawData)
  rmaExp <- exprs(rmaData)

  png(file.path(fig_dir, "rmaData-boxplot.png"),
      width = 14 * n_cel_files, height = 600)
  boxplot(rmaData, main = "RMA Data", cex.axis = 1.2, cex.lab = 1.2)
  dev.off()

  # Probe-level models -------------------------------------------------------
  plmFit <- fitProbeLevelModel(rawData)

  # NUSE and RLE -------------------------------------------------------------

  # Normalized Unscaled Standard-Errors (NUSE), should be distributed around 1
  png(file.path(fig_dir, "plmFit-NUSE.png"),
      width = 14 * n_cel_files, height = 600)
  NUSE(plmFit, main = "NUSE", cex.axis = 1.2, cex.lab = 1.2)
  dev.off()

  # Relative Log-Expression (RLE), should be distributed around 0
  png(file.path(fig_dir, "plmFit-RLE.png"),
      width = 14 * n_cel_files, height = 600)
  RLE(plmFit, main = "RLE", cex.axis = 1.2, cex.lab = 1.2)
  dev.off()

  # Array images -------------------------------------------------------------
  make_array_images(fig_dir, rawData, plmFit)

  # MA plot ------------------------------------------------------------------
  # M = x1 - x2
  # A = (x1 + x2) / 2
  ma_plots(fig_dir, rmaData)

  # Entrez IDs and Gene Symbols ----------------------------------------------
  mart_dat <- get_mart(rownames(rmaData), array_name, mart_file)

  # Select rows that have an Entrez Gene ID.
  mart_dat <- unique(mart_dat[!is.na(mart_dat$entrezgene), ])

  # Dictionaries to convert probeset IDs to useful identifiers.
  # Exclude probesets that map to multiple genes, or genes that map to
  # multiple gene symbols.
  probeset_to_entrez <- one_to_one(
    mart_dat, array_name, "entrezgene", exclude = TRUE)
  probeset_to_gene <- one_to_one(
    mart_dat, array_name, "external_gene_name", exclude = TRUE)
  entrez_to_gene <- one_to_one(
    mart_dat, "entrezgene", "external_gene_name", exclude = TRUE)

  # Average of multi-probe genes ---------------------------------------------

  # Save a copy of the RMA expression values for all probes.
  rmaExp_p <- rmaExp

  # Proportion of probesets that are assigned an Entrez Gene ID.
  prop_probesets_assigned <-
    100 * sum(rownames(rmaExp) %in% names(probeset_to_entrez)) / nrow(rmaExp)
  cat(sprintf(
    "%.2g%% of probesets are assigned exactly one Entrez Gene ID\n",
    prop_probesets_assigned
  ))

  # Exclude probes that are not assigned exactly one Entrez Gene ID.
  rmaExp <- rmaExp[rownames(rmaExp) %in% names(probeset_to_entrez), ]

  # Mean of multiple probesets assigned to the same Entrez Gene ID.
  rmaExp <- mean_by(
    rmaExp,
    probeset_to_entrez[as.character(rownames(rmaExp))]
  )
  entrez_ids <- rownames(rmaExp)
  gene_symbols <- entrez_to_gene[entrez_ids]

  # Save the workspace -------------------------------------------------------
  save(
    list = c(
        # Probeset and gene expression matrices.
        "rmaExp_p", "rmaExp",
        # Corresponding vectors of Entrez Gene IDs and gene symbols.
        "entrez_ids", "gene_symbols",
        # Dictionaries for converting ids.
        "probeset_to_entrez", "probeset_to_gene", "entrez_to_gene"
    ),
    file = out_file
  )

  cat("Finished successfully\n")
}

parse_args <- function() {
  # Parse these command-line options.
  spec <- matrix(c(
      'help'        , 'h', 0, "logical",
        "Print this help message",
      'cel_dir'     , 'c', 1, "character",
        "Directory with CEL files",
      'fig_dir'     , 'f', 1, "character",
        "Directory to write figures",
      'array_name'  , 'a', 1, "character",
        "Ensembl Biomart name of your array",
      'mart_file'   , 'm', 1, "character",
        "Output Biomart file with Entrez Gene IDs",
      'out_file'    , 'o', 1, "character",
        "Output RData file with expression and dictionaries"
  ), byrow = TRUE, ncol = 5)
  opt <- getopt(spec = spec)

  # Check for missing options.
  opts <- c("cel_dir", "fig_dir", "array_name", "mart_file", "out_file")
  all_args <- sapply(opts, function(x) !is.null(opt[[x]]))
  cat(sprintf("ERROR: missing option --%s\n", opts[!all_args]))

  # Print a help message and quit.
  if (!all(all_args) || !is.null(opt$help) || !length(commandArgs(TRUE))) {
    usage <- '
Usage: ./arrayqc.R [-h] -c <DIR> -f <DIR> -a <character> -m <FILE> -o <FILE>
    -h|--help          Print this help message.
    -c|--cel_dir       Directory with CEL files.
    -f|--fig_dir       Directory to write figures.
    -a|--array_name    Ensembl Biomart name of your array.
    -m|--mart_file     Output Biomart file with Entrez Gene IDs.
    -o|--out_file      Output RData file with expression and dictionaries.
'
    cat(usage)
    quit()
  }

  return(opt)
}

#' Take the mean of all columns of a matrix or dataframe, where rows are
#' aggregated by a vector of values. 100 times faster than stats::aggregate.
#' Depends on data.table and reshape2.
#' @param dat A numeric matrix or data.frame.
#' @param xs A vector of groups (e.g. gene names).
#' @return A data.table with the aggregated mean for each group.
#' @seealso stats::aggregate
mean_by <- function(dat, xs) {
  dat <- data.table(dat)
  dat$agg_var <- xs
  dat <- melt(dat, id.vars = "agg_var")
  dat <- dcast.data.table(
      dat, agg_var ~ variable, value.var = "value",
      fun.aggregate = mean, na.rm = TRUE
  )
  rowns <- dat$agg_var
  dat[ , agg_var := NULL]
  dat <- as.data.frame(dat)
  rownames(dat) <- rowns
  dat
}

#' Make images of all the arrays in a SummarizedExperiment object.
#' @param fig_dir The directory to create a folder called 'array-images'.
#' @param rawData A SummarizedExperiment object with raw microarray data.
#' @param plmFit The object returned by oligo::fitProbeLevelModel.
make_array_images <- function(fig_dir, rawData, plmFit) {
  array_image_files <- file.path(fig_dir, c(
    "rawData-images.jpg", "rawData-images-cuberoot.jpg",
    "plmFit-weights.jpg", "plmFit-residuals.jpg",
    "plmFit-sign_residuals.jpg"
  ))
  if (!all(file.exists(array_image_files))) {
    # Make the folder if it does not exist.
    dir.create(file.path(fig_dir, "array-images"), showWarnings = FALSE)
    
    height <- 960
    
    # Linear trends are easier to visualize after cuberoot transformation.
    cuberoot <- function(x) { y <- log2(x); sign(y) * abs(y) ^ (1/3) }
    
    # This will take several minutes to complete.
    for (i in 1:ncol(rawData)) {
      cat(sprintf("%s # Writing %i/%i\n", date(), i, ncol(rawData)))
      
      jpg.file <- sprintf("%s/array-images/rawData-image-%02d.jpg", fig_dir, i)
      jpeg(jpg.file, height = height, width = height)
      image(rawData[, i])
      dev.off()
      
      jpg.file <- sprintf(
        "%s/array-images/rawData-image-cuberoot-%02d.jpg", fig_dir, i)
      jpeg(jpg.file, height = height, width = height)
      image(rawData[, i], transfo = cuberoot)
      dev.off()
      
      jpg.file <- sprintf("%s/array-images/plmFit-weights-%02d.jpg", fig_dir, i)
      jpeg(jpg.file, height = height, width = height)
      image(plmFit, which = i, type = "weights")
      dev.off()
      
      jpg.file <- sprintf("%s/array-images/plmFit-residuals-%02d.jpg", fig_dir, i)
      jpeg(jpg.file, height = height, width = height)
      image(plmFit, which = i, type = "residuals")
      dev.off()
      
      jpg.file <- sprintf(
        "%s/array-images/plmFit-sign_residuals-%02d.jpg", fig_dir, i)
      jpeg(jpg.file, height = height, width = height)
      image(plmFit, which = i, type = "sign.residuals")
      dev.off()
    }
    # Depends on imagemagick.
    system(paste(
      "montage", file.path(fig_dir, "array-images/rawData-image-??.jpg"),
      "-geometry 180x180+0+0", file.path(fig_dir, "rawData-images.jpg")
    ))
    system(paste(
      "montage", file.path(fig_dir, "array-images/rawData-image-cuberoot-??.jpg"),
      "-geometry 180x180+0+0", file.path(fig_dir, "rawData-images-cuberoot.jpg")
    ))
    system(paste(
      "montage", file.path(fig_dir, "array-images/plmFit-weights-??.jpg"),
      "-geometry 180x180+0+0", file.path(fig_dir, "plmFit-weights.jpg")
    ))
    system(paste(
      "montage", file.path(fig_dir, "array-images/plmFit-residuals-??.jpg"),
      "-geometry 180x180+0+0", file.path(fig_dir, "plmFit-residuals.jpg")
    ))
    system(paste(
      "montage", file.path(fig_dir, "array-images/plmFit-sign_residuals-??.jpg"),
      "-geometry 180x180+0+0", file.path(fig_dir, "plmFit-sign_residuals.jpg")
    ))
  }
}

#' Create an MA plot for each array in a SummarizedExperiment object.
#' @param fig_dir The directory to create a folder called 'array-images'.
#' @param rmaData A SummarizedExperiment object with normalized data.
ma_plots <- function(fig_dir, rmaData) {
  ma_image_files <- file.path(fig_dir, "rmaData-ma.jpg")
  if (!all(file.exists(ma_image_files))) {
    # Make the folder if it does not exist.
    dir.create(file.path(fig_dir, "ma-plots"), showWarnings = FALSE)

    for (i in 1:ncol(rmaData)) {
      cat(sprintf("%s # Writing %i/%i\n", date(), i, ncol(rmaData)))

      maplot_file <- file.path(
        fig_dir, sprintf("ma-plots/rmaData-ma-%02d.jpg", i)
      )

      if (!file.exists(maplot_file)) {
        png(maplot_file, width = 800, height = 800)
        MAplot(rmaData, which = i)
        dev.off()
      }
    }
    system(paste(
      "montage", file.path(fig_dir, "ma-plots/rmaData-ma-??.jpg"),
      "-geometry 180x180+0+0", file.path(fig_dir, "rmaData-ma.jpg")
    ))
  }
}

#' Download Entrez Gene IDs and gene symbols from Ensembl Biomart, or load the
#' previously downloaded file.
#' @param probesets A character vector of probeset IDs.
#' @param array_name The name of the array in Ensembl Biomart.
#' @param mart_file Write output to this file. It will be gzipped.
#' @return A data.frame with Entrez Gene ID, symbol, and probeset ID.
get_mart <- function(probesets, array_name, mart_file) {
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  if (!file.exists(mart_file)) {
    mart_dat <- getBM(
        values = probesets,
        filters = c(array_name),
        attributes = c("entrezgene", "external_gene_name", array_name),
        mart = mart
    )
    # Write a compressed CSV file.
    mart_nogz <- str_replace(mart_file, ".gz", "")
    write.csv(mart_dat, mart_nogz)
    gzip(mart_nogz)
  } else {
    mart_dat <- read.csv(mart_file, row.names = 1)
  }
  return(mart_dat)
}

#' Create a one-to-one dictionary from one column to another.
#' @param dat A dataframe.
#' @param x The name or index of a column in the dataframe to use as names.
#' @param y The name or index of a column in the dataframe to use as values.
#' @param exclude If true, exclude x values that map to more than one y.
#' @return A named vector.
one_to_one <- function(dat, x, y, exclude = FALSE) {
  d <- unique(dat[ , c(x, y)])
  n_names <- length(unique(d[ , x]))
  n_values <- length(unique(d[ , y]))
  if (!exclude && n_names != n_values) {
    stop(sprintf("One %s maps to many %s", x, y))
  }
  if (exclude && n_names != n_values) {
    counts <- table(d[ , x])
    n_excluded <- sum(counts[counts > 1])
    warning(sprintf(
      "Excluding %i '%s' that map to many '%s'", n_excluded, x, y))
    singles <- names(counts[counts == 1])
    d <- d[d[ , x] %in% singles, ]
  }
  # Names.
  xs <- d[ , x]
  # Values.
  ys <- d[ , y]
  names(ys) <- xs
  ys
}

main()
