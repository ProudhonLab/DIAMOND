#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Concatenate fasta files by sequence ids.

@author: Marc Michel
@email: marc.michel@curie.fr
@project: https://github.com/michel-m
"""
import argparse

from Bio.SeqIO.FastaIO import SimpleFastaParser


def trim_sequences(fasta_file, n_bases):
    for title, sequence in SimpleFastaParser(fasta_file):
        trimmed_sequence = sequence[:-n_bases]
        umi = sequence[-n_bases:]
        print(f">{title}:{umi}\n{trimmed_sequence}")
        # print(f">{title}\n{umi}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file', type=argparse.FileType('r'), nargs=1,
                        help="fasta file to trim")
    parser.add_argument("-n", "--nbases",
                        type=int,
                        default=16,
                        help="number of bases to trim")
    args = parser.parse_args()
    trim_sequences(args.fasta_file[0], args.nbases)


if __name__ == '__main__':
    main()
