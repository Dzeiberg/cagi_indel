from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent))
from cagi_indel.generate_indels import generate_deletions, generate_insertions
import re


def test_generate_deletions():
    # Example usage
    seq = "TWYV"
    seq_identifier = "TEST"
    del_length = 2
    deletions = generate_deletions(seq, seq_identifier, del_length)
    # Expected deletions
    expected_deletions = {
        "TEST:p.Thr1del",
        "TEST:p.Trp2del",
        "TEST:p.Tyr3del",
        "TEST:p.Val4del",
        "TEST:p.Thr1_Trp2del",
        "TEST:p.Trp2_Tyr3del",
        "TEST:p.Tyr3_Val4del",
    }
    # Check if the generated deletions match the expected deletions
    assert (
        set(deletions) == expected_deletions
    ), f"Expected {expected_deletions}, but got {set(deletions)}"
    for deletion in deletions:
        if deletion not in expected_deletions:
            raise ValueError(f"Unexpected deletion: {deletion}")


def test_generate_insertions():
    # Example usage
    seq = "TWYV"
    seq_identifier = "Test"
    insertion_length = 2
    insertions = generate_insertions(seq, seq_identifier, insertion_length)
    assert (
        len(insertions) == 1260
    ), f"Expected 1260 insertions, but got {len(insertions)}"
    # Expected insertions
    assert validate_hgvs_insertions(
        insertions
    ), "Not all insertions match the HGVS format"


def validate_hgvs_insertions(hgvs_list):
    """
    Validates each string in the list against the HGVS protein insertion format
    with protein ID prefix.

    Args:
        hgvs_list (list of str): List of HGVS protein insertion strings.

    Returns:
        Bool: True if all strings match the format, False otherwise.
    """
    pattern = re.compile(
        r"^[a-zA-Z0-9_\.]+:p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)ins(([A-Z][a-z]{2})+)$"
    )

    return all([bool(pattern.match(s)) for s in hgvs_list])


if __name__ == "__main__":
    test_generate_insertions()
    test_generate_deletions()
    print("All tests passed!")
