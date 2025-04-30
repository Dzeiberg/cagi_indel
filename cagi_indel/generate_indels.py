from itertools import product
from Bio.Data.PDBData import protein_letters_3to1


def validate_input(seq: str):
    """
    Validate the input sequence to ensure it contains only valid amino acids.

    Args:
        seq (str): The sequence to validate.

    Raises:
        ValueError: If the sequence contains invalid amino acids.
    """
    valid_aa_single = set(protein_letters_3to1.values())
    for ref_index, ref_aa in enumerate(seq, start=1):
        if ref_aa not in valid_aa_single:
            raise ValueError(
                f"Invalid amino acid '{ref_aa}' at position {ref_index} in sequence '{seq}'"
            )


def generate_deletions(seq: str, seq_identifier: str, del_length: int):
    """
    Generate all possible amino acid deletions of a given length in a sequence.

    Args:
        seq (str): The original sequence.
        seq_identifier (str): The identifier for the sequence.
        del_length (int): The maximum length of deletions.

    Returns:
        list: A list of all possible deletion sequences.
    """
    protein_letters_1to3 = {
        one: three.title() for three, one in protein_letters_3to1.items()
    }
    deletions = list()
    single_position_template = f"{seq_identifier}:p.{{reference_aa}}{{aa_position}}del"
    position_range_template = f"{seq_identifier}:p.{{start_aa}}{{start_position}}_{{end_aa}}{{end_position}}del"
    # validate the sequence
    validate_input(seq)
    # generate deletions
    for length in range(1, del_length + 1):
        for pos in range(len(seq) - length + 1):
            start_aa = protein_letters_1to3[seq[pos]]
            start_position = pos + 1
            if length == 1:
                deletion = single_position_template.format(
                    reference_aa=start_aa, aa_position=start_position
                )
            else:
                end_aa = protein_letters_1to3[seq[pos + length - 1]]
                end_position = pos + length
                deletion = position_range_template.format(
                    start_aa=start_aa,
                    start_position=start_position,
                    end_aa=end_aa,
                    end_position=end_position,
                )
            deletions.append(deletion)
    return deletions


def generate_insertions(seq: str, seq_identifier: str, insertion_length: int):
    """
    Generate all possible amino acid insertions of a given maximum length in a sequence

    Args:
        seq (str): The original sequence.
        seq_identifier (str): The identifier for the sequence.
        insertion_length (int): The maximum length of insertions.

    Returns:
        list: A list of all possible insertion sequences.
    """
    # set the HGVS template for insertions
    insertion_template = f"{seq_identifier}:p.{{aa_range_start}}{{start_position}}_{{aa_range_end}}{{end_position}}ins{{insertion_sequence}}"
    valid_aa_single = set(protein_letters_3to1.values())
    insertions = list()
    protein_letters_1to3 = {
        one: three.title() for three, one in protein_letters_3to1.items()
    }
    # validate the sequence
    validate_input(seq)
    # generate insertions
    for length in range(1, insertion_length + 1):
        for pos in range(len(seq) - 1):
            start_position = pos + 1
            end_position = start_position + 1
            aa_range_start = protein_letters_1to3[seq[pos]]
            aa_range_end = protein_letters_1to3[seq[end_position - 1]]
            for insertion_sequence in product(valid_aa_single, repeat=length):
                insertion_sequence_str = "".join(
                    [protein_letters_1to3[single] for single in insertion_sequence]
                )
                insertion = insertion_template.format(
                    aa_range_start=aa_range_start,
                    start_position=start_position,
                    aa_range_end=aa_range_end,
                    end_position=end_position,
                    insertion_sequence=insertion_sequence_str,
                )
                insertions.append(insertion)
    return insertions
