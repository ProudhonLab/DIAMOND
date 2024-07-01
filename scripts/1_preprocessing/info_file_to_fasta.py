#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert cutadapt/atropos info file to fasta.

@author: Marc Michel
@email: marc.michel@curie.fr
@project: https://github.com/michel-m
"""
import argparse


ENDS = {
    "forward_umi": 4,
    "forward_primer": 5,
    "reverse_primer": 5,
    "reverse_umi": 6
}


def info_file_to_fasta(info_file, primer, end):
    for line in info_file:
        split_line = line.split('\t')
        read_id = split_line[0]
        read_sequence = split_line[ENDS[end]]
        if read_sequence:
            (print(f">{read_id}:{primer}\n{read_sequence}") if primer else
                print(f">{read_id}\n{read_sequence}"))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('info_file', type=argparse.FileType('r'),
                        help="input cutadapt info file")
    parser.add_argument('-p', '--primer',
                        type=str,
                        default=None,
                        help="primer name")
    parser.add_argument('-e', '--end',
                        type=str,
                        default="reverse_umi",
                        help="which end to recover, e.g. before 'forward' \
primer or after 'reverse' primer")
    args = parser.parse_args()
    info_file_to_fasta(args.info_file, args.primer, args.end)
