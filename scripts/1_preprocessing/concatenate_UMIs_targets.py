#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Concatenate fasta files by sequence ids.

@author: Marc Michel
@email: marc.michel@curie.fr
@project: https://github.com/michel-m
"""
import argparse
import os
import re
import sys
import threading

from Bio.SeqIO.FastaIO import SimpleFastaParser
from itertools import islice
from math import ceil


def concatenate_fastas(fasta_list, merged_sequences):
    for key in fasta_list[0].keys():
        merged_seq = ""
        for fasta in fasta_list:
            merged_seq += fasta[key]
        merged_sequences.append(f">{key}\n{merged_seq}")


def concatenate_fastas2(fasta0_slice, fasta1, end):
    for key in fasta0_slice.keys():
        try:
            merged_seq = (fasta0_slice[key] + fasta1[key] if end == "3prime" else fasta1[key] + fasta0_slice[key])
            sys.stdout.write(f">{key}\n{merged_seq}\n")
        except KeyError:
            continue
    sys.stdout.flush()


def load_sequences(fasta_files):
    seqid_pattern = re.compile("[a-zA-Z0-9_-]+:[0-9]+:[a-zA-Z0-9_-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+") # illumina see https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
    fasta_list = list()
    for fasta_file in fasta_files:
        sequences_dict = dict()
        for title, sequence in SimpleFastaParser(fasta_file):
            title_formated = title.split(' ')[0].replace('_', ':')
            title_formated = seqid_pattern.findall(title_formated)[0]
            sequences_dict[title_formated] = sequence
        fasta_list.append(sequences_dict)
    return fasta_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_files', type=argparse.FileType('r'), nargs='+', help="fasta files to concatenate")
    parser.add_argument("-t", "--threads", type=int, default=None, help="number of threads")
    parser.add_argument("-e", "--end", type=str, default="5prime", help="add fasta1 sequences to which end of fasta0 sequences? (choices: '5prime' or '3prime')")
    args = parser.parse_args()
    fasta_list = load_sequences(args.fasta_files)
    fasta0 = fasta_list[0]
    fasta1 = fasta_list[1]
    seq_count = len(fasta_list[0])
    cpu_count = args.threads if args.threads else os.cpu_count()
    slice_len = ceil(seq_count / cpu_count)
    for c in range(0, cpu_count):
        t = threading.Thread(target=concatenate_fastas2,
                             args=(dict(islice(fasta0.items(),
                                               c * slice_len,
                                               c * slice_len + slice_len)),
                                   fasta1, args.end))
        t.start()
    main_thread = threading.main_thread()
    for t in threading.enumerate():
        if t is not main_thread:
            t.join()

    exit()

    merged_sequences = list()
    concatenate_fastas(fasta_list, merged_sequences)
    print(*merged_sequences, sep="\n")
    exit()

    cpu_count = os.cpu_count()

    for c in cpu_count:
        t = threading.Thread(target=concatenate_fastas,
                             args=(fasta_list, merged_sequences,
                                   c * seq_count / cpu_count))
        t.start()
    main_thread = threading.main_thread()
    for t in threading.enumerate():
        if t is not main_thread:
            t.join()


if __name__ == '__main__':
    main()
