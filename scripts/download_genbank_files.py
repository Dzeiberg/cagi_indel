from Bio import Entrez, SeqIO
import time
from pathlib import Path
import fire
from tqdm import tqdm
import pandas as pd

Entrez.email = "dzeiberg@me.com"
Entrez.api_key = "d5cf9c7eb077239979c01fb25c5994bb7a08"


def download_genbank_files(
    mane_summary_file: Path, output_dir: Path, delay: float = 0.1
):
    """
    Download GenBank files for all transcripts listed in the MANE summary file.

    Args:
        mane_summary_file (Path): Path to the MANE summary file.
        output_dir (Path): Directory to save the downloaded GenBank files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    mane_summary_file = Path(mane_summary_file)
    if not mane_summary_file.exists():
        raise FileNotFoundError(f"File {mane_summary_file} does not exist.")
    summary_df = pd.read_csv(mane_summary_file, sep="\t")
    reseq_ids = summary_df.RefSeq_nuc.values
    for reseq_id in tqdm(reseq_ids, desc="Downloading GenBank files", unit="file"):
        try:
            with Entrez.efetch(
                db="nucleotide", id=reseq_id, rettype="gb", retmode="text"
            ) as handle:
                record = handle.read()
            output_file = output_dir / f"{reseq_id}.gb"
            with open(output_file, "w") as f:
                f.write(record)
            time.sleep(delay)
        except Exception as e:
            print(f"Error downloading {reseq_id}: {e}")
            continue


if __name__ == "__main__":
    # fire.Fire()
    download_genbank_files(
        mane_summary_file=Path(
            "/home/dzeiberg/cagi_indel/data/MANE_human_release_1.4/MANE/MANE.GRCh38.v1.4.summary.disease_genes.txt.gz"
        ),
        output_dir=Path(
            "/home/dzeiberg/cagi_indel/data/MANE_human_release_1.4/MANE/GenBank_files"
        ),
        delay=0.5,
    )
