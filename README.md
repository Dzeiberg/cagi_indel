# CAGI Indel Project

This repository contains code for the **CAGI Annotate Many Indels Project**, which focuses on predicting the effects of insertion and deletion (indel) variants.

## Repository Structure

- `data/`
    Contains external-input files and stores the resulting challenge varaint sets.

- `scripts/`
    Includes scripts for generating variant sets

## Requirements

- Python 3.13+
- Required Python packages are listed in `requirements.txt`.

## Installation

1. Clone the repository:
     ```bash
     git clone https://github.com/Dzeiberg/cagi_indel.git
     cd cagi_indel
     ```

2. Install dependencies:
     ```bash
     pip install -r requirements.txt
     ```

3. Download External Data
    ```bash
    bash scripts/download_external_data.sh
    ```

## Example Usage

Generate all single amino acid deletions:
```bash
python scripts/generate_inframe_files.py generate_inframe_file \
--mane_summary_filepath data/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz \
--mane_fasta_filepath data/MANE_human/release_1.4/MANE.GRCh38.v1.4.refseq_protein.faa.gz \
--insertion_or_deletion deletion \
--length 1 \
--output_file data/inframe_indel_files/MANE.GRCh38.v1.4.refseq_protein.singleDeletion.txt
```

## Contact

For questions or feedback, please contact [d.zeiberg@northeastern.edu].
