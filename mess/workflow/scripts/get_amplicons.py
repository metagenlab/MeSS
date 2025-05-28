import argparse
import sys
from Bio.Emboss import PrimerSearch
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


"""
functions
"""


# argparse
def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract amplicon squences from primersearch output file"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=argparse.FileType("r"),
        required=True,
        help="Primersearch results file",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=argparse.FileType("r"),
        required=True,
        help="Input fasta file",
    )
    parser.add_argument(
        "-p",
        "--primers",
        type=argparse.FileType("r"),
        required=True,
        help="Primers input table",
    )
    parser.add_argument(
        "-c",
        "--cut",
        action="store_true",
        help="Cut primers from amplicons sequences",
    )
    parser.add_argument(
        "-r",
        "--orient",
        action="store_true",
        help="Orient reverse sequences to forward",
    )
    parser.add_argument(
        "-m",
        "--minlen",
        type=int,
        required=False,
        help="Minimum amplicon sequence length",
        default=0,
    )
    parser.add_argument(
        "-M",
        "--maxlen",
        type=int,
        required=False,
        help="Maximum amplicon sequence length",
        default=3000,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        required=True,
        help="Output fasta file",
    )
    parser.add_argument(
        "-l",
        "--log",
        action="store_true",
        help="Print summary table to stderr",
    )
    parser.set_defaults(orient=False, cut=False, log=False)
    return parser.parse_args()


# IUPAC code mapping: code -> set of nucleotides
IUPAC_CODES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

# Reverse lookup for sets
SET_TO_IUPAC = {frozenset(v): k for k, v in IUPAC_CODES.items()}


def get_full_base_set(chars):
    """Expand each IUPAC/degenerate character into a full set of bases."""
    full_set = set()
    for char in chars:
        bases = IUPAC_CODES.get(char.upper(), {char})
        full_set.update(bases)
    return frozenset(full_set)


def convert_iupac(seq):
    # Replace brackets
    def bracket_replacer(match):
        content = match.group(1)
        base_set = get_full_base_set(content)
        return SET_TO_IUPAC.get(base_set, "N")  # Fallback to N

    seq = re.sub(r"\[([A-Z]+)\]", bracket_replacer, seq, flags=re.IGNORECASE)

    # Replace ? with N
    seq = seq.replace("?", "N")

    return seq


def parse_primersearch(results, seq2strand, minlen, maxlen):
    d = []
    for name, hits in results:
        for n, hit in enumerate(hits):
            amplimer = f"amplimer_{n + 1}"
            if hit.length < minlen or hit.length > maxlen:
                continue
            info = hit.hit_info.split("\n\t")
            mismatches = 0
            for match in info[1:]:
                mismatches += int(match.split("with ")[1].split(" ")[0])
                if "forward" in match:
                    fw_primer = convert_iupac(match.split()[0])
                    if seq2strand[fw_primer] == "forward":
                        forward_on_reverse = False
                    else:
                        forward_on_reverse = True
                    start = int(match.split("strand at ")[1].split(" ")[0])

                elif "reverse" in match:
                    rv_primer = convert_iupac(match.split()[0])
                    end = int(
                        (
                            match.split("strand at ")[1]
                            .split(" ")[0]
                            .replace("]", "")
                            .replace("[", "")
                        )
                    )

            d.append(
                {
                    "primer": name,
                    "seq_id": info[0].rstrip(),
                    "amplimer": amplimer,
                    "fw_primer": fw_primer,
                    "rv_primer": rv_primer,
                    "start_forward": start,
                    "end_reverse": end,
                    "total_mismatch": mismatches,
                    "length": hit.length,
                    "forward_match_on_reverse": forward_on_reverse,
                }
            )
        return pd.DataFrame.from_records(d)


def get_amplicon_record(row, seqid2record, cut, orient):
    amp_id = f"{row.seq_id}_{row.amplimer}_{row.primer}"
    amp_seq = seqid2record[row.seq_id].seq[row.start_forward - 1 : -row.end_reverse + 1]

    if cut:
        amp_seq = amp_seq[len(row.fw_primer) : len(amp_seq) - len(row.rv_primer)]

    if orient and row.forward_match_on_reverse:
        amp_seq = amp_seq.reverse_complement()
    return SeqRecord(amp_seq, id=amp_id, description=amp_id)


"""
main
"""


def main():
    args = parse_args()
    results = PrimerSearch.read(args.input).amplifiers.items()
    seq2strand = (
        pd.read_csv(
            args.primers,
            sep="\t",
            names=["primer", "forward", "reverse"],
        )
        .melt(
            id_vars=["primer"],
            value_vars=["forward", "reverse"],
            var_name="strand",
            value_name="seq",
        )
        .set_index("seq")
    ).to_dict()["strand"]

    seqid2record = SeqIO.to_dict(
        SeqIO.parse(
            args.fasta,
            "fasta",
        )
    )
    df = parse_primersearch(results, seq2strand, args.minlen, args.maxlen)
    if df.empty:
        amplicons = SeqRecord("", id="no_amps", description="no_amps")
    else:
        amplicons = df.apply(
            lambda row: get_amplicon_record(row, seqid2record, args.cut, args.orient),
            axis=1,
        ).tolist()
    SeqIO.write(
        amplicons,
        args.output,
        "fasta",
    )
    if args.log:
        print(df.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
