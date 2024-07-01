#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse # for arguments management
import re # for regular expressions
import pandas as pd # samples_info csv into dataframe
import itertools # construct all possible haplotypes
import matplotlib.pyplot as plt # to draw alignment plots
import numpy as np
import pickle

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

def get_cg_pos_count_per_sample(list_of_align_files, samples_info_file):

    samples_info = load_samples_info(samples_info_file)

    cg_pos_count_per_sample = {}

    with open(list_of_align_files) as input_list:

        # For each alignment file
        for align_file in input_list:
            
            align = fasta_iter(align_file.strip())

            # for each read
            for header, sequence in align:

                # only healthy plasma are used for position selection
                sample_id = re.search('^([A-Za-z0-9]+)_.+',header).group(1)
                try:
                    if samples_info.loc[sample_id]['biological_class'] != "healthy_plasma":
                        break
                except KeyError:
                    break
                
                if sample_id not in cg_pos_count_per_sample:
                    cg_pos_count_per_sample[sample_id] = {}
                    cg_pos_count_per_sample[sample_id]["seq_count"] = 0                
                cg_pos_count_per_sample[sample_id]["seq_count"] += 1

                all_cg_pos = [(m.start(),m.group()) for m in re.finditer('CG|TG|TA', sequence)]
                for cg_pos in all_cg_pos:
                    if cg_pos[0] not in cg_pos_count_per_sample[sample_id]:
                        cg_pos_count_per_sample[sample_id][cg_pos[0]] = {}
                        cg_pos_count_per_sample[sample_id][cg_pos[0]]['CG'] = 0
                        cg_pos_count_per_sample[sample_id][cg_pos[0]]['TG'] = 0
                        cg_pos_count_per_sample[sample_id][cg_pos[0]]['TA'] = 0
                    cg_pos_count_per_sample[sample_id][cg_pos[0]][cg_pos[1]] += 1
    
    return cg_pos_count_per_sample
    
def per_sample_to_global(cg_pos_count_per_sample):
    # convert into percentage and create global count (mean across samples)
    global_cg_pos_count = {}
    for sample_id in cg_pos_count_per_sample.keys():
        for cg_pos in cg_pos_count_per_sample[sample_id].keys():
            if cg_pos != "seq_count":
                if cg_pos not in global_cg_pos_count:
                    global_cg_pos_count[cg_pos] = {}
                    global_cg_pos_count[cg_pos]['CG'] = 0
                    global_cg_pos_count[cg_pos]['TG'] = 0
                    global_cg_pos_count[cg_pos]['TA'] = 0
                    global_cg_pos_count[cg_pos]['XX'] = 0
                global_cg_pos_count[cg_pos]['CG'] += (cg_pos_count_per_sample[sample_id][cg_pos]['CG'] * 100 / cg_pos_count_per_sample[sample_id]["seq_count"]) / len(cg_pos_count_per_sample.keys())
                global_cg_pos_count[cg_pos]['TG'] += (cg_pos_count_per_sample[sample_id][cg_pos]['TG'] * 100 / cg_pos_count_per_sample[sample_id]["seq_count"]) / len(cg_pos_count_per_sample.keys())
                global_cg_pos_count[cg_pos]['TA'] += (cg_pos_count_per_sample[sample_id][cg_pos]['TA'] * 100 / cg_pos_count_per_sample[sample_id]["seq_count"]) / len(cg_pos_count_per_sample.keys())
                global_cg_pos_count[cg_pos]['XX'] += ( (cg_pos_count_per_sample[sample_id]["seq_count"] - cg_pos_count_per_sample[sample_id][cg_pos]['CG'] - cg_pos_count_per_sample[sample_id][cg_pos]['TG'] - cg_pos_count_per_sample[sample_id][cg_pos]['TA']) * 100 / cg_pos_count_per_sample[sample_id]["seq_count"] ) / len(cg_pos_count_per_sample.keys())
    return global_cg_pos_count

def plot_alignments(global_cg_pos_count, threshold, outputDir, primer):
    # Plot CG/TG proportion by column
    sorted_pos = sorted(global_cg_pos_count)
    tg_prop = [global_cg_pos_count[cg_pos]['TG'] for cg_pos in sorted_pos]
    print(tg_prop)
    cg_prop = [global_cg_pos_count[cg_pos]['CG'] for cg_pos in sorted_pos]
    #ta_prop = [global_cg_pos_count[cg_pos]['TA'] for cg_pos in sorted_pos]
    #other_prop = [global_cg_pos_count[cg_pos]['XX'] for cg_pos in sorted_pos]
    threshold_line = plt.axhline(threshold, color='crimson', linestyle='dashed', linewidth=0.7)
    tg_bar = plt.bar(sorted_pos, tg_prop, color='cornflowerblue')
    cg_bar = plt.bar(sorted_pos, cg_prop, color='coral', bottom=tg_prop)
    #ta_bar = plt.bar(sorted_pos, ta_prop, color='turquoise', bottom=np.add(cg_prop, tg_prop))
    #other_bar = plt.bar(sorted_pos, other_prop, color='crimson', bottom=np.add(np.add(cg_prop, tg_prop),ta_prop))
    plt.xlabel("Position in sequence")
    plt.ylabel("Dinucleotides proportions (%)")
    #plt.legend((tg_bar, cg_bar, ta_bar, other_bar, threshold_line), ("TG", "CG", "TA", "Other Dinucleotide","Threshold"))
    plt.savefig(f"{outputDir}/alignment_plot_{primer}.png", dpi=300)

def select_positions(global_cg_pos_count, threshold):
    # select positions
    selected_pos = list()
    for cg_pos in global_cg_pos_count:
        cgtg_perc = global_cg_pos_count[cg_pos]['CG'] + global_cg_pos_count[cg_pos]['TG']
        if cgtg_perc >= threshold and global_cg_pos_count[cg_pos]['TG']/cgtg_perc <= 0.95:
            selected_pos.append(cg_pos)
    selected_pos.sort()
    return selected_pos

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, default=None, help="A list of align files")
    parser.add_argument("-o", "--outputDir", type=str, default=None, help="Output directory")
    parser.add_argument("-s", "--samples_info", type=str, default=None, help="A csv with samples info")
    parser.add_argument("-p", "--primer", type=str, default=None, help="Primer name")
    parser.add_argument("-t", "--threshold", type=str, default=None, help="Threshold for CG selection")
    args = parser.parse_args()

    cg_pos_count_per_sample = get_cg_pos_count_per_sample(args.input, args.samples_info)
    with open(f'{args.outputDir}/cg_pos_count_per_sample_{args.primer}.pickle', 'wb') as f:
        pickle.dump(cg_pos_count_per_sample, f)

    global_cg_pos_count = per_sample_to_global(cg_pos_count_per_sample)
    with open(f'{args.outputDir}/global_cg_pos_count_{args.primer}.pickle', 'wb') as f:
        pickle.dump(global_cg_pos_count, f)
    
    plot_alignments(global_cg_pos_count, float(args.threshold), args.outputDir, args.primer)

    selected_pos = select_positions(global_cg_pos_count, float(args.threshold))
    print(selected_pos)
    with open(f'{args.outputDir}/selected_pos_{args.primer}.pickle', 'wb') as f:
        pickle.dump(selected_pos, f)