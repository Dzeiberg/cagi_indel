import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
from cagi_indel.stops_and_frameshifts import generate_stop_gains


def test_generate_stop_gains():
    """
    Test the generate_stop_gains function.
    """
    # Define the path to the test GenBank file
    test_genbank_file = Path(__file__).parent.parent / "examples" / "NM_000518_5.gb"
    stop_gains = generate_stop_gains(test_genbank_file, "NM_000518.5")
    assert len(stop_gains) > 0, "No stop gains generated."
    assert "NM_000518.5:c.438T>A" in stop_gains


if __name__ == "__main__":
    test_generate_stop_gains()
    print("All tests passed.")
