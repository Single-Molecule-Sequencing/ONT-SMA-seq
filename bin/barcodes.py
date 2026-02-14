"""ONT native barcode sequences and classification functions.

Contains the 96 ONT native barcode sequences (24bp each) and functions for
barcode classification using edlib semi-global (HW) alignment.
"""

from __future__ import annotations

import edlib

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BARCODE_LENGTH = 24

BARCODES: dict[str, str] = {
    "nb01": "CACAAAGACACCGACAACTTTCTT",
    "nb02": "ACAGACGACTACAAACGGAATCGA",
    "nb03": "CCTGGTAACTGGGACACAAGACTC",
    "nb04": "TAGGGAAACACGATAGAATCCGAA",
    "nb05": "AAGGTTACACAAACCCTGGACAAG",
    "nb06": "GACTACTTTCTGCCTTTGCGAGAA",
    "nb07": "AAGGATTCATTCCCACGGTAACAC",
    "nb08": "ACGTAACTTGGTTTGTTCCCTGAA",
    "nb09": "AACCAAGACTCGCTGTGCCTAGTT",
    "nb10": "GAGAGGACAAAGGTTTCAACGCTT",
    "nb11": "TCCATTCCCTCCGATAGATGAAAC",
    "nb12": "TCCGATTCTGCTTCTTTCTACCTG",
    "nb13": "AGAACGACTTCCATACTCGTGTGA",
    "nb14": "AACGAGTCTCTTGGGACCCATAGA",
    "nb15": "AGGTCTACCTCGCTAACACCACTG",
    "nb16": "CGTCAACTGACAGTGGTTCGTACT",
    "nb17": "ACCCTCCAGGAAAGTACCTCTGAT",
    "nb18": "CCAAACCCAACAACCTAGATAGGC",
    "nb19": "GTTCCTCGTGCAGTGTCAAGAGAT",
    "nb20": "TTGCGTCCTGTTACGAGAACTCAT",
    "nb21": "GAGCCTCTCATTGTCCGTTCTCTA",
    "nb22": "ACCACTGCCATGTATCAAAGTACG",
    "nb23": "CTTACTACCCAGTGAACCTCCTCG",
    "nb24": "GCATAGTTCTGCATGATGGGTTAG",
    "nb25": "GTAAGTTGGGTATGCAACGCAATG",
    "nb26": "CATACAGCGACTACGCATTCTCAT",
    "nb27": "CGACGGTTAGATTCACCTCTTACA",
    "nb28": "TGAAACCTAAGAAGGCACCGTATC",
    "nb29": "CTAGACACCTTGGGTTGACAGACC",
    "nb30": "TCAGTGAGGATCTACTTCGACCCA",
    "nb31": "TGCGTACAGCAATCAGTTACATTG",
    "nb32": "CCAGTAGAAGTCCGACAACGTCAT",
    "nb33": "CAGACTTGGTACGGTTGGGTAACT",
    "nb34": "GGACGAAGAACTCAAGTCAAAGGC",
    "nb35": "CTACTTACGAAGCTGAGGGACTGC",
    "nb36": "ATGTCCCAGTTAGAGGAGGAAACA",
    "nb37": "GCTTGCGATTGATGCTTAGTATCA",
    "nb38": "ACCACAGGAGGACGATACAGAGAA",
    "nb39": "CCACAGTGTCAACTAGAGCCTCTC",
    "nb40": "TAGTTTGGATGACCAAGGATAGCC",
    "nb41": "GGAGTTCGTCCAGAGAAGTACACG",
    "nb42": "CTACGTGTAAGGCATACCTGCCAG",
    "nb43": "CTTTCGTTGTTGACTCGACGGTAG",
    "nb44": "AGTAGAAAGGGTTCCTTCCCACTC",
    "nb45": "GATCCAACAGAGATGCCTTCAGTG",
    "nb46": "GCTGTGTTCCACTTCATTCTCCTG",
    "nb47": "GTGCAACTTTCCCACAGGTAGTTC",
    "nb48": "CATCTGGAACGTGGTACACCTGTA",
    "nb49": "ACTGGTGCAGCTTTGAACATCTAG",
    "nb50": "ATGGACTTTGGTAACTTCCTGCGT",
    "nb51": "GTTGAATGAGCCTACTGGGTCCTC",
    "nb52": "TGAGAGACAAGATTGTTCGTGGAC",
    "nb53": "AGATTCAGACCGTCTCATGCAAAG",
    "nb54": "CAAGAGCTTTGACTAAGGAGCATG",
    "nb55": "TGGAAGATGAGACCCTGATCTACG",
    "nb56": "TCACTACTCAACAGGTGGCATGAA",
    "nb57": "GCTAGGTCAATCTCCTTCGGAAGT",
    "nb58": "CAGGTTACTCCTCCGTGAGTCTGA",
    "nb59": "TCAATCAAGAAGGGAAAGCAAGGT",
    "nb60": "CATGTTCAACCAAGGCTTCTATGG",
    "nb61": "AGAGGGTACTATGTGCCTCAGCAC",
    "nb62": "CACCCACACTTACTTCAGGACGTA",
    "nb63": "TTCTGAAGTTCCTGGGTCTTGAAC",
    "nb64": "GACAGACACCGTTCATCGACTTTC",
    "nb65": "TTCTCAGTCTTCCTCCAGACAAGG",
    "nb66": "CCGATCCTTGTGGCTTCTAACTTC",
    "nb67": "GTTTGTCATACTCGTGTGCTCACC",
    "nb68": "GAATCTAAGCAAACACGAAGGTGG",
    "nb69": "TACAGTCCGAGCCTCATGTGATCT",
    "nb70": "ACCGAGATCCTACGAATGGAGTGT",
    "nb71": "CCTGGGAGCATCAGGTAGTAACAG",
    "nb72": "TAGCTGACTGTCTTCCATACCGAC",
    "nb73": "AAGAAACAGGATGACAGAACCCTC",
    "nb74": "TACAAGCATCCCAACACTTCCACT",
    "nb75": "GACCATTGTGATGAACCCTGTTGT",
    "nb76": "ATGCTTGTTACATCAACCCTGGAC",
    "nb77": "CGACCTGTTTCTCAGGGATACAAC",
    "nb78": "AACAACCGAACCTTTGAATCAGAA",
    "nb79": "TCTCGGAGATAGTTCTCACTGCTG",
    "nb80": "CGGATGAACATAGGATAGCGATTC",
    "nb81": "CCTCATCTTGTGAAGTTGTTTCGG",
    "nb82": "ACGGTATGTCGAGTTCCAGGACTA",
    "nb83": "TGGCTTGATCTAGGTAAGGTCGAA",
    "nb84": "GTAGTGGACCTAGAACCTGTGCCA",
    "nb85": "AACGGAGGAGTTAGTTGGATGATC",
    "nb86": "AGGTGATCCCAACAAGCGTAAGTA",
    "nb87": "TACATGCTCCTGTTGTTAGGGAGG",
    "nb88": "TCTTCTACTACCGATCCGAAGCAG",
    "nb89": "ACAGCATCAATGTTTGGCTAGTTG",
    "nb90": "GATGTAGAGGGTACGGTTTGAGGC",
    "nb91": "GGCTCCATAGGAACTCACGCTACT",
    "nb92": "TTGTGAGTGGAAAGATACAGGACC",
    "nb93": "AGTTTCCATCACTTCAGACTTGGG",
    "nb94": "GATTGTCCTCAAACTGCCACCTAC",
    "nb95": "CCTGTCTGGAAGAAGAATGGACTT",
    "nb96": "CTGAACGGTCATAGAGTCCACCAT",
}

# ---------------------------------------------------------------------------
# Complement table
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence containing only uppercase A, C, G, T characters.

    Returns
    -------
    str
        The reverse complement of *seq*.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def classify_barcode(
    segment: str,
    expected_barcodes: dict[str, str],
) -> dict[str, str | int | float]:
    """Classify a read segment against a set of expected barcodes.

    Uses edlib semi-global (HW) alignment to find the best-matching barcode.
    In HW mode the query (barcode) is fully aligned while the target (segment)
    can be partially consumed, which handles flanking DNA gracefully.

    Parameters
    ----------
    segment : str
        The read segment to classify (may be longer than 24bp).
    expected_barcodes : dict[str, str]
        Mapping of barcode_id -> sequence to search against.

    Returns
    -------
    dict
        ``barcode_id`` : str -- best-matching barcode name
        ``edit_distance`` : int -- edit distance of best alignment
        ``confidence`` : float -- 1.0 - (edit_distance / BARCODE_LENGTH)
    """
    best_id: str | None = None
    best_ed: int = BARCODE_LENGTH + 1  # worse than any real alignment

    for bc_id, bc_seq in expected_barcodes.items():
        result = edlib.align(bc_seq, segment, mode="HW", task="distance")
        ed: int = result["editDistance"]
        if ed < best_ed:
            best_ed = ed
            best_id = bc_id
            if ed == 0:
                break  # can't do better than perfect

    confidence = max(0.0, 1.0 - best_ed / BARCODE_LENGTH)

    return {
        "barcode_id": best_id,
        "edit_distance": best_ed,
        "confidence": confidence,
    }
