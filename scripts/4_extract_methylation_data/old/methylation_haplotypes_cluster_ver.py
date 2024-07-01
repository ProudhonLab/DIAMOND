#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Extract methylation haplotypes from a .fasta alignment file.

@author: Marc Michel
@email: marc.michel@curie.fr
@project: https://github.com/michel-m
"""
import argparse
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path

mpl.use("Agg")
if plt.get_backend() == 'Qt5Agg':
    from matplotlib.backends.qt_compat import QtWidgets
    qApp = QtWidgets.QApplication(sys.argv)
    plt.matplotlib.rcParams['figure.dpi'] = qApp.desktop().physicalDpiX()

ACCEPTABLE_PERC_OF_EXTRACTED_CG = 20
SAMPLES_INFO_FILE = "/home/genouest/cnrs_umr6074/kdasilva/proudhon_lab/psl1_meth/scripts/marcLegacy/sample_biological-class_annotation.csv"
FIG_DPI = 300


def load_samples_info(samples_info_file_path):
    return pd.read_csv(samples_info_file_path, index_col='sample')


def plot_dinucleotides(fasta_alignment_file, primer_name, bioclass_filter):
    samples_dict = dict()
    sample_id_pattern = re.compile("^[A-Za-z0-9]+")
    cg_tg_pattern = re.compile("CG|TG|TA")
    col_candidates = list()
    samples_info = load_samples_info(SAMPLES_INFO_FILE)
    for title, sequence in SimpleFastaParser(fasta_alignment_file):
        sample_id = sample_id_pattern.findall(title)[0]
        try:
            biological_class = samples_info.loc[sample_id]['biological_class']
        except KeyError:
            biological_class = "unknown_biological_class"
        if bioclass_filter is not None:
            if not biological_class.endswith(bioclass_filter):
                continue
        if sample_id not in samples_dict:
            samples_dict[sample_id] = {'cg_count': dict(), 'seq_count': 0}
        samples_dict[sample_id]['seq_count'] += 1
        it = cg_tg_pattern.finditer(sequence)
        while True:
            try:
                cg = next(it)
            except StopIteration:
                break
            else:
                if cg.start() not in samples_dict[sample_id]['cg_count']:
                    samples_dict[sample_id]['cg_count'][cg.start()] = {
                        'CG': 0, 'TG': 0, 'TA': 0}
                samples_dict[sample_id]['cg_count'][cg.start()
                                                    ][cg.group()] += 1
    global_cg_count = dict()
    for sample_id in samples_dict:
        for key in sorted(samples_dict[sample_id]['cg_count']):
            cg_per = samples_dict[sample_id]['cg_count'][key]['CG'] * \
                100 / samples_dict[sample_id]['seq_count']
            tg_per = samples_dict[sample_id]['cg_count'][key]['TG'] * \
                100 / samples_dict[sample_id]['seq_count']
            ta_per = samples_dict[sample_id]['cg_count'][key]['TA'] * \
                100 / samples_dict[sample_id]['seq_count']
            try:
                global_cg_count[key]['cg_per'] += cg_per
                global_cg_count[key]['tg_per'] += tg_per
                global_cg_count[key]['ta_per'] += ta_per
            except KeyError:
                global_cg_count[key] = dict()
                global_cg_count[key]['cg_per'] = cg_per
                global_cg_count[key]['tg_per'] = tg_per
                global_cg_count[key]['ta_per'] = ta_per

    samples_count = len(samples_dict.keys())
    sorted_cg_count = sorted(global_cg_count)
    for key in sorted_cg_count:
        col_per = (global_cg_count[key]['cg_per'] +
                   global_cg_count[key]['tg_per']) / samples_count
        if col_per >= ACCEPTABLE_PERC_OF_EXTRACTED_CG:
            if (global_cg_count[key]['tg_per'] / samples_count * 100 / col_per
                    <= 95):  # TODO: use global variable for 95% threshold for real TG identification
                col_candidates.append(key)

    # Plot CG/TG proportion by column
    cg_prop = [global_cg_count[k]['cg_per'] /
               samples_count for k in sorted_cg_count]
    tg_prop = [global_cg_count[k]['tg_per'] /
               samples_count for k in sorted_cg_count]
    ta_prop = [global_cg_count[k]['ta_per'] /
               samples_count for k in sorted_cg_count]
    threshold_line = plt.axhline(ACCEPTABLE_PERC_OF_EXTRACTED_CG,
                                 color='crimson',
                                 linestyle='dashed',
                                 linewidth=0.7)
    tg_bar = plt.bar(sorted_cg_count, tg_prop, color='cornflowerblue')
    cg_bar = plt.bar(sorted_cg_count, cg_prop, color='coral', bottom=tg_prop)
    ta_bar = plt.bar(sorted_cg_count, ta_prop, color='turquoise',
                     bottom=np.add(cg_prop, tg_prop))
    plt.xlabel("Position in sequence")
    plt.ylabel("Dinucleotides proportions (%)")
    plt.legend((tg_bar, cg_bar, ta_bar, threshold_line),
               ("TG", "CG", "TA", "Threshold"))
    plt.savefig(f"{fasta_alignment_file.name[:-6]}.{bioclass_filter}.png",
                dpi=FIG_DPI)


def methylation_haplotypes(fasta_alignment_file, primer_name, annotate,
                           include_ta, include_all):
    samples_list = list()
    sample_id_pattern = re.compile("^[A-Za-z0-9]+")
    cg_tg_pattern = re.compile("CG|TG")
    seq_count = 0
    cg_count = dict()
    col_candidates = list()
    samples_info = load_samples_info(SAMPLES_INFO_FILE)
    for title, sequence in SimpleFastaParser(fasta_alignment_file):
        sample_id = sample_id_pattern.findall(title)[0]
        if sample_id not in samples_list:
            samples_list.append(sample_id)
        try:
            if samples_info.loc[sample_id]['biological_class'] != "healthy_plasma":
                continue
        except KeyError:
            continue
        seq_count += 1
        it = cg_tg_pattern.finditer(sequence)
        while True:
            try:
                cg = next(it)
            except StopIteration:
                break
            else:
                if cg.start() not in cg_count:
                    cg_count[cg.start()] = {'CG': 0, 'TG': 0}
                cg_count[cg.start()][cg.group()] += 1
    for key in sorted(cg_count):
        cg_per = cg_count[key]['CG'] * 100 / seq_count
        tg_per = cg_count[key]['TG'] * 100 / seq_count
        col_per = cg_per + tg_per
        if col_per >= ACCEPTABLE_PERC_OF_EXTRACTED_CG:
            if tg_per * 100 / col_per <= 95:
                col_candidates.append(key)

    samples_dict = dict()
    for sample in samples_list:
        samples_dict[sample] = dict()
        for x in map(''.join, itertools.product(
                '01', repeat=len(col_candidates))):
            samples_dict[sample][x] = 0

    fasta_alignment_file.seek(0)
    for title, sequence in SimpleFastaParser(fasta_alignment_file):
        sample_id = sample_id_pattern.findall(title)[0]
        haplotype = list()
        for col in col_candidates:
            if sequence[col:col + 2] == "CG":
                haplotype.append("1")
            elif sequence[col:col + 2] == "TG":
                haplotype.append("0")
            elif include_ta is True and sequence[col:col + 2] == "TA":
                haplotype.append("0")
            elif include_all is True:
                haplotype.append("0")
            else:
                haplotype.append("2")
        if "2" not in haplotype:
            samples_dict[sample_id][''.join(haplotype)] += 1

    samples_info = load_samples_info(SAMPLES_INFO_FILE)
    annotation = f"{primer_name}_"
    if annotate:
        print("sample,biological_class,", end='')
    print(','.join([annotation + k for k in samples_dict[sample_id].keys()]))
    for sample in samples_list:
        try:
            biological_class = samples_info.loc[sample]['biological_class']
        except KeyError:
            biological_class = "unknown_biological_class"
        haplotypes_dict = samples_dict[sample]
        hap_sum = sum(haplotypes_dict.values())
        # for haplotype, count in haplotypes_dict.items():
        #     print(f"{haplotype},{count / haplotypes_sum}")
        if annotate:
            print(f"{sample},{biological_class},", end='')
        # TODO: div by zero
        try:
            print(','.join(
                [str(v / hap_sum) for v in haplotypes_dict.values()]))
        except ZeroDivisionError:
            print("ZeroDivisionError")


def methylation_cg_methyl(fasta_alignment_file, primer_name, annotate,
                          include_ta, include_all):
    samples_list = list()
    sample_id_pattern = re.compile("^[A-Za-z0-9]+")
    cg_tg_pattern = re.compile("CG|TG")
    seq_count = 0
    cg_count = dict()
    col_candidates = list()
    samples_info = load_samples_info(SAMPLES_INFO_FILE)
    for title, sequence in SimpleFastaParser(fasta_alignment_file):
        sample_id = sample_id_pattern.findall(title)[0]
        if sample_id not in samples_list:
            samples_list.append(sample_id)
        try:
            if samples_info.loc[sample_id][
                    'biological_class'] != "healthy_plasma":
                continue
        except KeyError:
            continue
        seq_count += 1
        it = cg_tg_pattern.finditer(sequence)
        while True:
            try:
                cg = next(it)
            except StopIteration:
                break
            else:
                if cg.start() not in cg_count:
                    cg_count[cg.start()] = {'CG': 0, 'TG': 0}
                cg_count[cg.start()][cg.group()] += 1
    for key in sorted(cg_count):
        cg_per = cg_count[key]['CG'] * 100 / seq_count
        tg_per = cg_count[key]['TG'] * 100 / seq_count
        col_per = cg_per + tg_per
        if col_per >= ACCEPTABLE_PERC_OF_EXTRACTED_CG:
            if tg_per * 100 / col_per <= 95:
                col_candidates.append(key)

    samples_dict = dict()
    for sample in samples_list:
        samples_dict[sample] = dict()
        for col in col_candidates:
            samples_dict[sample][col] = {'seq_count': 0, 'cg_count': 0}

    fasta_alignment_file.seek(0)
    for title, sequence in SimpleFastaParser(fasta_alignment_file):
        sample_id = sample_id_pattern.findall(title)[0]
        for col in col_candidates:
            if sequence[col:col + 2] == "CG":
                samples_dict[sample_id][col]['cg_count'] += 1
                samples_dict[sample_id][col]['seq_count'] += 1
            elif include_all is True and sequence[col:col + 2] != "CG":
                samples_dict[sample_id][col]['seq_count'] += 1
            elif include_ta is True and sequence[col:col + 2] == "TA":
                samples_dict[sample_id][col]['seq_count'] += 1
            elif sequence[col:col + 2] == "TG":
                samples_dict[sample_id][col]['seq_count'] += 1

    samples_info = load_samples_info(SAMPLES_INFO_FILE)
    if annotate:
        print("sample,biological_class,", end='')
    print(','.join([f"{primer_name}_cg{i + 1}_pos{k}"
                    for i, k in enumerate(samples_dict[sample_id].keys())]))
    for sample in samples_list:
        try:
            biological_class = samples_info.loc[sample]['biological_class']
        except KeyError:
            biological_class = "unknown_biological_class"
        methyl_level_dict = samples_dict[sample]
        # for haplotype, count in haplotypes_dict.items():
        #     print(f"{haplotype},{count / haplotypes_sum}")
        if annotate:
            print(f"{sample},{biological_class},", end='')
        # TODO: div by zero
        try:
            print(','.join([str(v['cg_count'] / v['seq_count'])
                            for v in methyl_level_dict.values()]))
        except ZeroDivisionError:
            print("ZeroDivisionError")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=argparse.FileType('r'), nargs=1,
                        help="input .fasta alignment file")
    parser.add_argument("-p", "--primer",
                        type=str,
                        default=None,
                        help="primer name")
    parser.add_argument("-a", "--annotate",
                        help="annotate with sample name & biological class",
                        default=False,
                        action="store_true")
    parser.add_argument("--include_ta",
                        help="include TA towards unmethylated score (along with TG)",
                        default=False,
                        action="store_true")
    parser.add_argument("--include_all",
                        help="include all dinucleotides that isn't CG towards unmethylated score (along with TG)",
                        default=False,
                        action="store_true")
    parser.add_argument("--cg_methyl",
                        help="format methylation data as average per CG",
                        default=False,
                        action="store_true")
    parser.add_argument("--plot",
                        help="plot dinucleotides CG, TG, TA proportions",
                        default=False,
                        action="store_true")
    parser.add_argument("--filter",
                        type=str,
                        default=None,
                        help="filter by biological class (endswith)")
    args = parser.parse_args()

    if args.cg_methyl:
        methylation_cg_methyl(args.file[0], args.primer, args.annotate,
                              args.include_ta, args.include_all)
    elif args.plot is True:
        plot_dinucleotides(args.file[0], args.primer, args.filter)
    else:
        methylation_haplotypes(args.file[0], args.primer, args.annotate,
                               args.include_ta, args.include_all)
