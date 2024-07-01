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
        output_file.write("classes,features,primers,subsampling,sensitivity,sensitivity_95CI_lower,sensitivity_95CI_upper\n")
        for dir in os.listdir(args.inputDir):
            
            folder_name_splitted = dir.split("/")[-1].split(".")            
            path_to_csv = f'{args.inputDir}{dir}/sensitivity_specificity.csv'

            if os.path.isfile(path_to_csv):
                output_file.write(f"{'_'.join(folder_name_splitted[1].split('_')[2:])},{folder_name_splitted[2]},{folder_name_splitted[3]},{folder_name_splitted[4]},")
                data = pd.read_csv(path_to_csv)
                sensitivity = float(data.loc[data["specificity"]==0.99,"sensitivity"])
                sensitivity_95CI_lower = float(data.loc[data["specificity"]==0.99,"sensitivity_95CI_lower"])
                sensitivity_95CI_upper = float(data.loc[data["specificity"]==0.99,"sensitivity_95CI_upper"])
                output_file.write(f"{sensitivity},{sensitivity_95CI_lower},{sensitivity_95CI_upper}\n")

