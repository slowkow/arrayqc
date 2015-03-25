#!/usr/bin/env bash
# download.sh
# Kamil Slowikowski

# Download the CEL files.
wget -N ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58203/suppl/GSE58203_RAW.tar

# Extract the CEL files.
tar xf GSE58203_RAW.tar
