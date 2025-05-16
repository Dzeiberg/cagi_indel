import pandas as pd
from pathlib import Path
from fire import Fire


def get_disease_genes(
    clinvar_variant_summary_filepath: Path,
    mane_summary_filepath: Path,
    output_filepath: Path,
) -> None:
    """
    Get the subset of the MANE Select summary file that contains genes with at least one >= 1-star P/LP annotation in ClinVar.
    """
    # Read the ClinVar variant summary file
    clinvar_df = pd.read_csv(
        clinvar_variant_summary_filepath, sep="\t", compression="gzip"
    )
    pathogenic_clinical_significance_labels = {
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic",
    }
    zero_star_review_statuses = {
        "no assertion criteria provided",
        "no classification for the single variant",
        "no classifications from unflagged records",
        "no classification provided",
    }
    # Filter for pathogenic/likely pathogenic variants with >= 1-star review status
    clinvar_pathogenic = clinvar_df[
        (
            clinvar_df["ClinicalSignificance"].isin(
                pathogenic_clinical_significance_labels
            )
        )
        & (~clinvar_df["ReviewStatus"].isin(zero_star_review_statuses))
    ]
    # Get the unique HGNC_ID values from the filtered DataFrame
    hgnc_ids = clinvar_pathogenic["HGNC_ID"].unique()
    # Read the MANE Select summary file
    mane_df = pd.read_csv(mane_summary_filepath, sep="\t", compression="gzip")
    # Filter the MANE Select summary file for the disease gene HGNC_IDs
    mane_disease_genes = mane_df[mane_df.HGNC_ID.isin(hgnc_ids)]
    # Save the filtered MANE Select summary file to the output path
    mane_disease_genes.to_csv(
        output_filepath, sep="\t", index=False, compression="gzip"
    )
    print(f"Filtered MANE Select summary file saved to {output_filepath}")


if __name__ == "__main__":
    Fire(get_disease_genes)
    # Example usage:
    # get_disease_genes(
    #     clinvar_variant_summary_filepath=Path("data/clinvar/variant_summary_2025-05.txt.gz"),
    #     mane_summary_filepath=Path("data/MANE_human_release_1.4/MANE/MANE.GRCh38.v1.4.summary.txt.gz"),
    #     output_filepath=Path("data/MANE_human_release_1.4/MANE/MANE.GRCh38.v1.4.summary.disease_genes.txt.gz")
    # )
