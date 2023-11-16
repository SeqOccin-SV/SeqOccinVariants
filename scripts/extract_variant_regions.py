#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import sys

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pysam import VariantFile, FastaFile
import pandas as pd

def eprint(*args, **kwargs):
    "Simple utility function to print on stderr"
    print(*args,  file=sys.stderr, **kwargs)


def svlen(record):
    """
    Returns the length of a given variant based on the size of the alleles
    """
    return max(len(record.alts[0]),len(record.ref)) - 1


def valid_record(record): 
    """
    Check if the record is valid :
       A valid record must have a reference and an alternate allele entirely
       descdibed as a sequence of nucleotide (eg. not <INS>).
       The first bases of ref an alternate alleles must match
    Returns true if valid and False otherwise
    """
    if ("<" not in record.alts[0] and
        record.alts[0][0] == record.ref[0]):
        return True
    return False


def get_region_seq(record, fasta, flanking=5000, chrom_prefix=""):
    """
    Retrieve the region surrounding a given variant, ie the genome segment
    defined on the variant chromosome by the interval [start - flanking, stop + flanking]
    """
    chrom = record.chrom
    pos = record.pos
    stop = record.stop
    start = max(1, pos - flanking)
    end = stop + flanking
    sequence = fasta.fetch(reference=chrom, start=start-1, end=end)

    # the position of the variant in the sequence
    pos_seq = min(pos, flanking + 1)
    name = "%s" % record.id
    region = "%s:%d-%d" % (chrom, start, end)
    seq_record = SeqRecord(Seq(sequence).upper(), id=name,
                           description="region=%s pos_var=%d" % (region,pos_seq))
    chrom = chrom_prefix + chrom
    coordinates = [name, chrom, start, end]
    return seq_record, coordinates


def main(vcffile, genome, outfasta, outtsv, minsize, maxsize, flankingsize, chrom_prefix):
    fasta = FastaFile(genome)
    vcf_in = VariantFile(vcffile)
    seq_records = []
    regions = []
    for record in vcf_in:
        if (valid_record(record) and svlen(record) > minsize
                and svlen(record) < maxsize):
            seq, coord = get_region_seq(record, fasta, flanking=flankingsize,
                                        chrom_prefix=chrom_prefix)
            seq_records.append(seq)
            regions.append(coord)
    SeqIO.write(seq_records, outfasta, "fasta")
    df = pd.DataFrame(regions, columns=['variant', 'chrom', 'start', 'end'])
    df.to_csv(outtsv, sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract SV surrounding sequence")
    parser.add_argument('-v', '--vcf', required=True,
                        help='the vcf file file')
    parser.add_argument('-g', '--genome', required=True,
                        help='the genome fasta file')
    parser.add_argument('--outfasta', required=True,
                        help='the output fasta file')
    parser.add_argument('--outtsv', required=True,
                        help='the output tsv file')
    parser.add_argument('--min-size', type=int, default=0,
                        help='the variant minimum size')
    parser.add_argument('--max-size', type=int, default=1000000,
                        help='the variant maximum size')
    parser.add_argument('--flanking-size', type=int, default=100,
                        help='the flanking')
    parser.add_argument('--chrom-prefix', type=str, default="",
                        help='the chromosome prefix')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.vcf, args.genome, args.outfasta, args.outtsv,
         args.min_size, args.max_size, args.flanking_size, args.chrom_prefix)
