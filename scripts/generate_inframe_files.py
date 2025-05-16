from pathlib import Path
import sys
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from typing import List
import gzip
from tqdm import tqdm
from fire import Fire
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
from cagi_indel.generate_indels import generate_deletions, generate_insertions


def parse_mane_select(summary_filename: Path, fasta_filepath: Path) -> List[SeqRecord]:
    """
    Parse the MANE Select file and return a list of SeqRecord objects.

    Args:
        summary_filepath (Path): Path to the MANE Select summary file.
        fasta_filepath (Path): Path to the MANE FASTA file.
    Returns:
        List[SeqRecord]: List of SeqRecord objects.
    """
    summary_filename = Path(summary_filename)
    if not summary_filename.exists():
        raise FileNotFoundError(f"File {summary_filename} does not exist.")
    summary_df = pd.read_csv(
        summary_filename,
        sep="\t",
        compression="gzip" if summary_filename.suffix == ".gz" else None,
    )
    disease_gene_transcripts = summary_df.RefSeq_prot.values
    # Check if the MANE FASTA file exists
    fasta_filepath = Path(fasta_filepath)
    if not fasta_filepath.exists():
        raise FileNotFoundError(f"File {fasta_filepath} does not exist.")
    # Check if the MANE FASTA file is gzipped
    if fasta_filepath.suffix == ".gz":
        # If the file is gzipped, use gzip to open it
        # and parse it with Bio.SeqIO
        with gzip.open(fasta_filepath, "rt") as handle:
            records_index = {record.id: record for record in parse(handle, "fasta")}
    else:
        # If the file is not gzipped, parse it directly
        records_index = {record.id: record for record in parse(fasta_filepath, "fasta")}
    missing_transcripts = set(disease_gene_transcripts) - set(records_index.keys())
    if missing_transcripts:
        print(
            f"Warning: The following transcripts are missing in the FASTA file: {', '.join(list(map(str, missing_transcripts)))}"
        )
    return [
        records_index[transcript_id]
        for transcript_id in disease_gene_transcripts
        if transcript_id in records_index
    ]


def generate_inframe_file(
    mane_summary_filepath: Path,
    mane_fasta_filepath: Path,
    insertion_or_deletion: str,
    length: int,
    output_file: Path,
):
    """
    Generate in-frame insertions of a given length for each sequence in the MANE summary file.

    Args:
        mane_summary_filepath (Path): Path to the MANE summary file.
        mane_fasta_filepath (Path): Path to the MANE FASTA file.
        insertion_or_deletion (str): Type of in-frame mutation to generate ("insertion" or "deletion").
        length (int): The maximum length of insertions/deletions.
        output_file (Path): Path to the output file.
    """
    records = parse_mane_select(mane_summary_filepath, mane_fasta_filepath)
    counter = 0
    if insertion_or_deletion == "insertion":
        fn = lambda seq, seq_id: generate_insertions(seq, seq_id, length)
    elif insertion_or_deletion == "deletion":
        fn = lambda seq, seq_id: generate_deletions(seq, seq_id, length)
    else:
        raise ValueError(
            "insertion_or_deletion must be either 'insertion' or 'deletion'"
        )
    with open(output_file, "w") as f:
        for record in tqdm(records):
            seq = str(record.seq)
            seq_identifier = str(record.id)
            insertions = fn(seq, seq_identifier)
            counter += len(insertions)
            for insertion in insertions:
                f.write(f"{insertion}\n")
    print(
        f"{counter:,d} in-frame {insertion_or_deletion}s generated and saved to {output_file}"
    )


if __name__ == "__main__":
    Fire()
