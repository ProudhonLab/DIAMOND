#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse # for arguments management
import os
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputDir", type=str, default=None, help="Input directory")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output file (and directory)")
    args = parser.parse_args()

    with open(args.output,"w") as output_file:
        output_file.write("classes,features,primers,subsampling,accuracy,accuracy_lower_ci,accuracy_upper_ci,f1,f1_lower_ci,f1_upper_ci,roc_auc,roc_auc_lower_ci,roc_auc_upper_ci,pr_auc,pr_auc_lower_ci,pr_auc_upper_ci\n")
        for dir in os.listdir(args.inputDir):

            folder_name_splitted = dir.split("/")[-1].split(".")
            if "multicancer" not in dir:
                path_to_csv = f'{args.inputDir}{dir}/scores.csv'
            else:
                path_to_csv = "" #TODO

            if os.path.isfile(path_to_csv):
                output_file.write(f"{'_'.join(folder_name_splitted[2].split('_')[2:])},{folder_name_splitted[3]},{folder_name_splitted[4]},{folder_name_splitted[5]},")
                data = pd.read_csv(path_to_csv)
                accuracy = float(data.loc[0,"score"])
                accuracy_lower_ci = float(data.loc[0,"lower_ci"])
                accuracy_upper_ci = float(data.loc[0,"upper_ci"])
                f1 = float(data.loc[1,"score"])
                f1_lower_ci = float(data.loc[0,"lower_ci"])
                f1_upper_ci = float(data.loc[1,"upper_ci"])
                roc_auc = float(data.loc[2,"score"])
                roc_lower_ci = float(data.loc[2,"lower_ci"])
                roc_upper_ci = float(data.loc[2,"upper_ci"])
                pr_auc = float(data.loc[3,"score"])
                pr_lower_ci = float(data.loc[3,"lower_ci"])
                pr_upper_ci = float(data.loc[3,"upper_ci"])
                output_file.write(f"{accuracy},{accuracy_lower_ci},{accuracy_upper_ci},{f1},{f1_lower_ci},{f1_upper_ci},{roc_auc},{roc_lower_ci},{roc_upper_ci},{pr_auc},{pr_lower_ci},{pr_upper_ci}\n")