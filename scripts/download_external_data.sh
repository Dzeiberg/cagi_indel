#!/bin/bash

# This script downloads external data files required for the project.
# It is assumed that the script is run from the root directory of the project.
# Ensure the script is executable
# chmod +x scripts/download_external_data.sh

# download MANE-select data
base_url="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4"
mkdir -p data/MANE_human/release_1.4/

for file in MANE.GRCh38.v1.4.refseq_protein.faa.gz MANE.GRCh38.v1.4.summary.txt.gz; do
  curl -o "data/MANE_human/release_1.4/${file}" "${base_url}/${file}"
done

# Download ClinVar variant_summary
mkdir -p data/clinvar/
curl -o "data/clinvar/variant_summary_2025-05.txt.gz" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2025-05.txt.gz"
