#!/usr/bin/env bash
# example.sh
# Kamil Slowikowski

opt=(
  --cel_dir data-raw
  --fig_dir figures
  --array_name affy_hg_u133_plus_2
  --mart_file data-raw/mart.csv.gz
  --out_file data/GSE58203.RData 
)
Rscript R/arrayqc.R ${opt[*]} &> arrayqc.log

