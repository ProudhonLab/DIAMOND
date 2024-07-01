#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse # for arguments management
import re # for regular expressions
import itertools # construct all possible haplotypes
import pandas as pd # samples_info csv into dataframe
import pickle
import os # check if inputfile exists

def load_samples_info(samples_info_file_path):
    return pd.read_csv(samples_info_file_path, index_col='sample')

def fasta_iter(fasta_name):
    # https://www.biostars.org/p/710/#1412
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in itertools.groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq

def get_methyl_and_haplo(input_file, selected_pos_pickle):

    with open(selected_pos_pickle, 'rb') as f:
        selected_pos = pickle.load(f)
    haplotypes = [f"h{x}" for x in map(''.join, itertools.product('01', repeat=len(selected_pos)))]

    d_methyl = {}
    d_haplo = {}
    for pos in selected_pos:
        d_methyl[pos] = {'CG' : 0, 'TG' : 0, 'TA' : 0, 'XX' : 0}
    for haplo in haplotypes:
        d_haplo[haplo] = 0

    align = fasta_iter(input_file)

    # for each read
    for header, sequence in align:

        sample_id = re.search('^([A-Za-z0-9]+)_.+',header).group(1)

        haplotype = list()
        for pos in selected_pos:
            dinucl = sequence[pos:(pos+2)]
            if dinucl == 'CG':
                d_methyl[pos]['CG'] += 1
                haplotype.append('1')
            elif dinucl == 'TG':
                d_methyl[pos]['TG'] += 1
                haplotype.append('0')
            elif dinucl == 'TA':
                d_methyl[pos]['TA'] += 1
                haplotype.append('0')
            else:
                d_methyl[pos]['XX'] += 1
                haplotype.append('0')
        haplotype = ''.join(haplotype)
        d_haplo[f"h{haplotype}"] += 1

    return d_methyl,d_haplo,sample_id

def write_methyl_csv(d_methyl, sample_id, outputDir, primer, mode):
    # write methylation results in csv

    selected_pos = sorted(d_methyl)

    ## headers
    headers = "sample,biological_class,"
    headers += ','.join([f"{primer}_cg{i + 1}_pos{k}" for i, k in enumerate(selected_pos)])
    vs_all = open(f"{outputDir}/cg_methyl.{sample_id}.{primer}.csv", "w")
    vs_all.write(f"{headers}\n")

    line = f"{sample_id},"
    try:
        biological_class = samples_info.loc[sample_id]['biological_class']
    except KeyError:
        biological_class = "unknown_biological_class"
    line += f"{biological_class},"

    for pos in selected_pos:
        pos_counts = d_methyl[pos]
        if mode=="cgtg":
            ref_sum = pos_counts['CG']+pos_counts['TG']
        else: # all
            ref_sum = pos_counts['CG']+pos_counts['TG']+pos_counts['TA']+pos_counts['XX']
        cg_perc = pos_counts['CG']/ref_sum
        line += f"{cg_perc},"

    vs_all.write(f"{line[:-1]}\n")
    vs_all.close()

def write_haplo_csv(d_haplo, sample_id, outputDir, primer):
    # write haplotype results in csv

    haplotypes = sorted(d_haplo)

    ## headers
    headers = "sample,biological_class,"
    headers += ','.join([f"{primer}_{x}" for x in haplotypes])
    vs_all = open(f"{outputDir}/haplotypes.{sample_id}.{primer}.csv", "w")
    vs_all.write(f"{headers}\n")

    line = f"{sample_id},"
    try:
        biological_class = samples_info.loc[sample_id]['biological_class']
    except KeyError:
        biological_class = "unknown_biological_class"
    line += f"{biological_class},"

    sum_haplo = sum(d_haplo.values())
    for haplo in haplotypes:
        haplo_vs_all = d_haplo[haplo]/sum_haplo
        line += f"{haplo_vs_all},"

    vs_all.write(f"{line[:-1]}\n")
    vs_all.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default=None, help="Input file")
    parser.add_argument("-o", "--outputDir", type=str, default=None, help="Output directory")
    parser.add_argument("-s", "--samplesInfo", type=str, default=None, help="Samples info file")
    parser.add_argument("-c", "--selectedPos", type=str, default=None, help="Pickle of selected positions for the primer")
    parser.add_argument("-m", "--mode", type=str, default="all", help="Method for computing methylation levels: all or cgtg [default: all]")
    args = parser.parse_args()

    print(f"ARGUMENTS:")
    print(f"--input: {args.input}")
    print(f"--outputDir: {args.outputDir}")
    print(f"--samplesInfo: {args.samplesInfo}")
    print(f"--selectedPos: {args.selectedPos}\n")
    print(f"--mode: {args.mode}\n")

    inputfile_split = args.input.split(".")
    sample_id = inputfile_split[-3]
    primer = inputfile_split[-2]

    print(f"DERIVED FROM ARGUMENTS:")
    print(f"- sample_id: {sample_id}")
    print(f"- primer: {primer}\n")

    print(f"PROGRESS - Loading {args.samplesInfo}\n")
    samples_info = load_samples_info(args.samplesInfo)

    if not os.path.exists(args.input) or os.stat(args.input).st_size == 0:

        print(f"STOP - File does not exist or is empty\n")

        with open(args.selectedPos, 'rb') as f:
            selected_pos = pickle.load(f)
        haplotypes = [f"h{x}" for x in map(''.join, itertools.product('01', repeat=len(selected_pos)))]

        ## methyl
        headers = "sample,biological_class,"
        headers += ','.join([f"{primer}_cg{i + 1}_pos{k}" for i, k in enumerate(selected_pos)])
        vs_all = open(f"{args.outputDir}/cg_methyl.{sample_id}.{primer}.csv", "w")
        vs_all.write(f"{headers}\n")

        line = f"{sample_id},"
        try:
            biological_class = samples_info.loc[sample_id]['biological_class']
        except KeyError:
            biological_class = "unknown_biological_class"
        line += f"{biological_class},"

        for pos in selected_pos:
            line += f"NA,"

        vs_all.write(f"{line[:-1]}\n")
        vs_all.close()

        ## haplotypes
        headers = "sample,biological_class,"
        headers += ','.join([f"{primer}_{x}" for x in haplotypes])
        vs_all = open(f"{args.outputDir}/haplotypes.{sample_id}.{primer}.csv", "w")
        vs_all.write(f"{headers}\n")

        line = f"{sample_id},"
        try:
            biological_class = samples_info.loc[sample_id]['biological_class']
        except KeyError:
            biological_class = "unknown_biological_class"
        line += f"{biological_class},"

        for haplo in haplotypes:
            line += f"NA,"

        vs_all.write(f"{line[:-1]}\n")
        vs_all.close()

    else:

        print(f"PROGRESS - Get methylation levels and haplotypes\n")
        d_methyl,d_haplo,sample_id = get_methyl_and_haplo(args.input,args.selectedPos)

        print(f"PROGRESS - Writing results\n")
        write_methyl_csv(d_methyl,sample_id,args.outputDir,primer, args.mode)
        write_haplo_csv(d_haplo,sample_id,args.outputDir,primer)