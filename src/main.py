import sys
import argparse
from fasta import read_fasta
from fastq import scan_reads
from sam import ssam_record
from st import mccreight_st_construction as st


def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching using a suffix tree")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()

    genome = {
        name: st(seq) for name, seq in read_fasta(args.genome).items()
    }
    for read_name, read_seq in scan_reads(args.reads):
        for chr_name, chr_st in genome.items():
            for i in chr_st.search(read_seq):
                ssam_record(sys.stdout,
                            read_name, chr_name,
                            i, f"{len(read_seq)}M",
                            read_seq)


if __name__ == '__main__':
    main()
