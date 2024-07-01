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

    list_dir = os.listdir(args.inputDir)
    list_dir.remove("logs")

    for i,dir in enumerate(list_dir):

        folder_name_splitted = dir.split("/")[-1].split(".")
        path_to_csv = f'{args.inputDir}{dir}/features_importance.csv'
        name=f"{'_'.join(folder_name_splitted[2].split('_')[2:])}.{folder_name_splitted[3]}.{folder_name_splitted[4]}.{folder_name_splitted[5]}"

        if i==0:

            if os.path.isfile(path_to_csv):
                data = pd.read_csv(path_to_csv)
                data = data.drop(columns=['stdev','median','total_score','mean_rank'])
                data.rename(columns={'mean':name}, inplace=True)

        else:
            if os.path.isfile(path_to_csv):
                data2 = pd.read_csv(path_to_csv)
                data2 = data2.drop(columns=['stdev','median','total_score','mean_rank'])
                data2.rename(columns={'mean':name}, inplace=True)
                data = data.merge(data2, how='outer', on='feature')

    data.to_csv(args.output)
