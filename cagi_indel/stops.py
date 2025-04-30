from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path


def parse_transcript_genbank_file(file_path: Path):
    """
    Parse the GenBank file of a transcript

    Args:
        file_path (Path): Path to the GenBank file.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File {file_path} does not exist.")
    if not file_path.is_file():
        raise ValueError(f"Path {file_path} is not a file.")
    if file_path.suffix != ".gb":
        raise ValueError(f"File {file_path} is not a GenBank file.")
    # Parse the GenBank file
    with open(file_path, "r") as file:
        seq_record = next(SeqIO.parse(file, "genbank"))
    return seq_record


def get_cds(seq_record: SeqRecord) -> SeqRecord:
    """
    Extract the coding sequence (CDS) from a SeqRecord object.

    Args:
        seq_record (SeqRecord): The SeqRecord object containing the sequence data.

    Returns:
        SeqRecord: A new SeqRecord object containing the CDS.
    """
    # Extract the CDS from the SeqRecord object
    cd_features = [feature for feature in seq_record.features if feature.type == "CDS"]
    if len(cd_features) == 0:
        raise ValueError("No CDS feature found in the GenBank file.")
    if len(cd_features) > 1:
        raise ValueError("Multiple CDS features found in the GenBank file.")
    cds_feature = cd_features[0]
    cds_sequence = cds_feature.extract(seq_record.seq)
    if not isinstance(cds_sequence, Seq):
        raise ValueError("Extracted CDS is not a Seq object.")
    cds_record = SeqRecord(cds_sequence, id=seq_record.id, description="CDS")
    if not isinstance(cds_record, SeqRecord):
        raise ValueError("Extracted CDS is not a SeqRecord object.")
    return cds_record


def generate_stop_gains(transcript_genbank_file: Path, seq_identifier: str) -> list:
    """
    Generate all possible stop gains from a SeqRecord object.

    Args:
        transcript_seq_record (SeqRecord): The SeqRecord object for the gene transcript containing the sequence data.
        seq_identifier (str): The identifier for the sequence.

    Returns:
        list: A list of all possible stop gain variants
    """
    stop_gains = []
    # Parse the GenBank file
    seq_record = parse_transcript_genbank_file(transcript_genbank_file)
    # Extract the CDS from the SeqRecord object
    cds_record = get_cds(seq_record)
    # Iterate over each codon in the CDS
    template = (
        f"{seq_identifier}:c.{{position}}{{reference_sequence}}>{{alternate_sequence}}"
    )
    for codon_num, codon_start_idx in enumerate(range(0, len(cds_record.seq), 3)):  # type: ignore
        # Get the current codon
        codon_start_pos = codon_start_idx + 1
        codon = cds_record.seq[codon_start_idx : codon_start_idx + 3]  # type: ignore
        is_snv_nonsense, (
            position,
            reference_sequence,
            alternate_sequence,
        ) = get_snv_nonsense(codon, codon_start_pos)
        if is_snv_nonsense:
            # Generate the HGVS notation for the stop gain
            hgvs = template.format(
                position=position,
                reference_sequence=reference_sequence,
                alternate_sequence=alternate_sequence,
            )
            # Append the HGVS notation to the list
            stop_gains.append(hgvs)
    return stop_gains


def get_snv_nonsense(codon, codon_start_pos):
    """
    Check if the codon is 1 edit away from a stop codon. If so, return the position of the edit and the reference and alternate bases.

    Args:
        codon (str): The codon to check.
        codon_start_pos (int): The start position of the codon in the sequence.


    Returns:
        bool: True if the codon is 1 edit away from a stop codon, False otherwise.
        tuple: A tuple containing the position of the edit, the reference base, and the alternate base.
    """
    stop_codons = ["TAA", "TAG", "TGA"]
    # check if the codon is 1 edit away from a stop codon
    edit_distances = [
        sum([ai != bi for ai, bi in zip(codon, stop_codon)])
        for stop_codon in stop_codons
    ]
    # check if the codon is 1 edit away from a stop codon
    if min(edit_distances) == 1:
        # get the position of the edit
        stop_codon = stop_codons[edit_distances.index(min(edit_distances))]
        # get the reference and alternate
        position, reference_sequence, alternate_sequence = get_substitution(
            codon, stop_codon, codon_start_pos
        )
        return True, (position, reference_sequence, alternate_sequence)
    return False, (None, None, None)


def get_substitution(codon, stop_codon, codon_start_pos):
    """
    Get the position of the edit, the reference base, and the alternate base.
    Args:
        codon (str): The codon to check.
        stop_codon (str): The stop codon to check against.
        codon_start_pos (int): The start position of the codon in the sequence.
    Returns:
        tuple: A tuple containing the position of the edit, the reference base, and the alternate base.
    """
    # get the position at which codon and stop_codon differ
    edit_position = [i for i, (a, b) in enumerate(zip(codon, stop_codon)) if a != b]
    if len(edit_position) != 1:
        raise ValueError("Codon and stop codon differ at more than one position.")
    # get the position of the edit
    edit_position = edit_position[0]
    # get the reference and alternate
    reference_sequence = codon[edit_position]
    alternate_sequence = stop_codon[edit_position]
    position = codon_start_pos + edit_position
    return position, reference_sequence, alternate_sequence
