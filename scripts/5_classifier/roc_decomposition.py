"""
This script calculates the Area Under the Receiver Operating Characteristic Curve (AUC-ROC) 
for different cancer types and stages. It also plots the ROC curves for each cancer type.

Usage:
    python script_name.py -d input_data.csv -s sample_list.csv -r number_of_runs -oC output_csv.csv -c color_files.csv -oP output_plot_prefix
"""

import argparse
import pandas as pd
import sys
from sklearn.metrics import (roc_curve, auc)
import numpy as np
import matplotlib.pyplot as plt
import os
import re

def get_99sensitivity(fpr):
    """
 Find the index of the false positive rate (FPR) that is closest to 0.99 specificity.

 Parameters:
 fpr (array): Array of false positive rates.

 Returns:
 int: Index of the FPR that is closest to 0.99 specificity.
 """
    target_specificity = 0.99
    closest_index = None
    closest_diff = float('inf')
    for i, spec in enumerate(1 - fpr):
        if abs(spec - target_specificity) <= closest_diff:
            closest_diff = abs(spec - target_specificity)
            closest_index = i
    return closest_index

def get_auc_roc_per_run(K_samples, H_samples, method, K_type, nruns, K_tag):
    """
    Calculate the AUC-ROC for a given cancer type and stage.

    Parameters:
    K_samples (DataFrame): DataFrame containing the predictions for the cancer samples.
    H_samples (DataFrame): DataFrame containing the predictions for the healthy samples.
    method (str): Method to use for subsampling (not used in this script).
    K_type (str): Type of cancer (e.g. 'BRC M0', 'OVC M+', etc.).
    nruns (int): Number of runs to perform.
    K_tag (str): the name of the cancer category to search in blind_test_file if blind_test_file exists

    Returns:
    list: List of dictionaries containing the AUC-ROC and sensitivity at 99% specificity for each run.
    dict: Dictionary containing the average ROC curves for each blind status.
    """
    res = list()
    col_list = [f'proba_breast_cancer_plasma_run_{n}' for n in range(0, nruns)]
    sel_col = [col_name for col_name in col_list if K_samples[col_name].count() >= 2]
    avg_roc_curves = {'seen': {'fprs': [], 'tprs': []}, 'blind': {'fprs': [], 'tprs': []}}
    for col_name in sel_col:
        k_pred = K_samples[col_name]
        k_pred = k_pred[k_pred.isna() == False]  # remove NaN values
        k_status = np.full(len(k_pred), 1)

        h_pred = H_samples[col_name]
        h_pred = h_pred[h_pred.isna() == False]

        if method == 'subsample':
            h_pred = np.random.choice(h_pred, size=len(k_pred))
            h_status = np.full(len(h_pred), 0)
        else:
            h_status = np.full(len(h_pred), 0)

        pred = np.concatenate([k_pred, h_pred])
        status = np.concatenate([k_status, h_status])

        fpr, tpr, _ = roc_curve(status, pred)

        ### blind status
        run_n=int(re.search(r'proba_breast_cancer_plasma_run_(\d+)', col_name).group(1))
        if "blind_test_file" in globals() and K_tag in list(blind_test_file.columns) :
            blind_status = blind_test_file.loc[blind_test_file['run'] == run_n, K_tag].values[0]
        else :
            blind_status = 'seen'

        res.append({'condition': K_type, 'auc': auc(fpr, tpr), 'sensitivity_at_99_spec': tpr[get_99sensitivity(fpr)], 'blind_status' : blind_status})
        avg_roc_curves[blind_status]['fprs'].append(fpr)
        avg_roc_curves[blind_status]['tprs'].append(tpr)

    # Calculate the average ROC curve
    for blind_status in avg_roc_curves:
        if avg_roc_curves[blind_status]['fprs']:  # Check if the list is not empty
            n_points = 100
            avg_fpr = np.linspace(0, 1, n_points)
            avg_tpr = np.zeros(n_points)
            for fpr, tpr in zip(avg_roc_curves[blind_status]['fprs'], avg_roc_curves[blind_status]['tprs']):
                interp_tpr = np.interp(avg_fpr, fpr, tpr)
                avg_tpr += interp_tpr
            avg_tpr /= len(avg_roc_curves[blind_status]['fprs'])
            avg_roc_curves[blind_status] = (avg_fpr, avg_tpr)
        else:
            avg_roc_curves[blind_status] = (None, None)  # or any other default value

    return res, avg_roc_curves



def plot_roc_curves(full_res, avg_roc_curves, color_files, output_plots):
    """
    Plot ROC curves for each k_type in a single plot with colors defined by a color file.
    Only k_types that are present in both the color file and the full_res condition column are plotted.

    Parameters:
    full_res (pd.DataFrame): DataFrame containing the results of the ROC curve analysis
    avg_roc_curves (dict): Dictionary containing the average ROC curves for each k_type and blind status
    color_files (list): List of paths to files containing the color definitions for each k_type
    output_plots (str): Prefix for the output plot file names

    Returns:
    None
    """
    # Iterate over each color file
    for color_file in color_files:
        # Load color definitions from file
        colors = pd.read_csv(color_file)

        # Create a dictionary mapping each k_type to its corresponding color
        colors_dict = dict(zip(colors['k_type'], colors['color']))

        if "blind_test_file" in globals():
            fig_seen, ax_seen = plt.subplots(figsize=(7, 7), constrained_layout=True)
            ax_seen.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

            for k_type in colors['k_type']:
                if k_type in full_res['condition'].values:
                    seen_res = full_res[(full_res['condition'] == k_type) & (full_res['blind_status'] == 'seen')]
                    if not seen_res.empty:
                        auc_score = seen_res['auc'].mean()
                        auc_score = round(auc_score, 2)
                        ax_seen.plot(avg_roc_curves[k_type]['seen'][0], avg_roc_curves[k_type]['seen'][1], label=f'{k_type} (AUC = {auc_score})', color=colors_dict[k_type])

            ax_seen.set(xlabel="1 - specificity", ylabel="Sensitivity", title=f"Mean ROC curve (Seen)")
            ax_seen.axis("square")
           # ax_seen.tick_params(axis='x', labelrotation=-90)

            ax_seen.legend(loc="lower right")
           # plt.xticks( list(map(lambda x: x/100, range(0,101,1))))
          #  plt.yticks( list(map(lambda x: x/100, range(0,101,1))))
            plot_name_seen = f"{output_plots}_seen_{os.path.splitext(os.path.basename(color_file))[0]}.png"
            plt.savefig(plot_name_seen)

            fig_blind, ax_blind = plt.subplots(figsize=(7, 7), constrained_layout=True)
            ax_blind.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)
            #ax_blind.tick_params(axis='x', labelrotation=-90)
            #plt.xticks( list(map(lambda x: x/100, range(0,101,1))))
            #plt.yticks( list(map(lambda x: x/100, range(0,101,1))))

            for k_type in colors['k_type']:
                if k_type in full_res['condition'].values:
                    blind_res = full_res[(full_res['condition'] == k_type) & (full_res['blind_status'] == 'blind')]
                    if not blind_res.empty:
                        auc_score = blind_res['auc'].mean()
                        auc_score = round(auc_score, 2)
                        ax_blind.plot(avg_roc_curves[k_type]['blind'][0], avg_roc_curves[k_type]['blind'][1], label=f'{k_type} (AUC = {auc_score})', color=colors_dict[k_type])

            ax_blind.set(xlabel="1 - specificity", ylabel="Sensitivity", title=f"Mean ROC curve (Blind)")
            ax_blind.axis("square")
            ax_blind.legend(loc="lower right")

            plot_name_blind = f"{output_plots}_blind_{os.path.splitext(os.path.basename(color_file))[0]}.png"
            plt.savefig(plot_name_blind)
        else:
            # Load color definitions from file

            # Create a figure and axis
            fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
            ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

            # Iterate over each k_type in the color file
            for k_type in colors['k_type']:
                if k_type in full_res['condition'].values:
                    auc_score = full_res[full_res['condition'] == k_type]['auc'].mean()
                    auc_score = round(auc_score, 2)
                    ax.plot(avg_roc_curves[k_type]['seen'][0], avg_roc_curves[k_type]['seen'][1], label=f'{k_type} (AUC = {auc_score})', color=colors_dict[k_type])

            # Set plot title and labels
            ax.set(
                xlabel="1 - specificity",
                ylabel="Sensitivity",
                title=f"Mean ROC curve",
            )
            ax.axis("square")
            ax.legend(loc="lower right")
            #ax.tick_params(axis='x', labelrotation=-90)

            #plt.xticks( list(map(lambda x: x/100, range(0,101,1))))
            #plt.yticks( list(map(lambda x: x/100, range(0,101,1))))
            # Save plot
            plot_name = f"{output_plots}_{os.path.splitext(os.path.basename(color_file))[0]}.png"
            plt.savefig(plot_name)




if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='ROC DECOMPO', usage="""This script calculates the Area Under the Receiver Operating Characteristic Curve (AUC-ROC)                                      for different cancer types and stages. It also plots the ROC curves for each cancer type.                                     """)
    parser.add_argument("-d","--data", type=str, default=None, help="input full_predictions.csv of all model that need to be decompose.")
    parser.add_argument("-s", "--sample_list", type=str, default="/home/klaus/Documents/24_01_16_heatmap_by_cancer_type/S3paper.csv", help=" a sample list with metastatic, cancer types and stages status.")
    parser.add_argument("-r", "--runs", type=int, default=1000, help="number of runs.")
    parser.add_argument("-oC", "--output_csv", type=str, default=None, help="output csv file.")
    parser.add_argument("-c", "--color_files", type=str, nargs='+', help="one or more color files for plotting.")
    parser.add_argument("-oP", "--output_plots", type=str, default=None, help="output plot file path and prefix, None if you dont want to plot ROC curves.")
    parser.add_argument("-b","--blind_test_file",type=str,default=None,help="a file indicating for each run the status (blind/seen) of each cancer type in training set.")
    parser.add_argument("-D","--sep",type=str,default=';',help="file delimiteur.")
    args = parser.parse_args()
    
    
    # Load data
    test_data=pd.read_csv(args.data)
    sample_list = pd.read_csv(args.sample_list,sep=args.sep)
    full=pd.merge(sample_list,test_data,on='Sample_ID')

    if args.blind_test_file :
        blind_test_file=pd.read_csv(args.blind_test_file)
        print('blind test file loaded')
    
    # Initialize
    all_res=list()
    avg_roc_curves = {}
    healthy=full[full['Disease_status']=='healthy']
    ## BRC M0
    sub=full[((full['Disease_status']=='breast_cancer') | (full['Disease_status']=='early_breast_cancer')) & (full['Metastasis_status']!='M+')]
    result,  avg_roc_curves['BRC M0']  = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='BRC M0',nruns=args.runs,K_tag='early_breast_cancer')
    all_res.append(pd.DataFrame(result))
    
    
    ## breast_cancer M+
    sub=full[((full['Disease_status']=='breast_cancer') | (full['Disease_status']=='early_breast_cancer')) & (full['Metastasis_status']!='M0')]
    result,avg_roc_curves['BRC M+'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='BRC M+',nruns=args.runs,K_tag='breast_cancer')
    all_res.append(pd.DataFrame(result))

    
    ## OVC M0
    sub=full[(full['Disease_status']=='ovarian_cancer')  & (full['Metastasis_status']!='M+')]
    result,avg_roc_curves['OVC M0'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='OVC M0',nruns=args.runs,K_tag='early_ovarian_cancer')
    all_res.append(pd.DataFrame(result))
    
    ## OVC M+
    
    sub=full[((full['Disease_status']=='ovarian_cancer')  & (full['Metastasis_status']!='M0'))]
    result,  avg_roc_curves['OVC M+'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='OVC M+',nruns=args.runs,K_tag='ovarian_cancer')
    all_res.append(pd.DataFrame(result))
    
    ## CRC M + 
    
    sub=full[(full['Disease_status']=='colorectal_cancer')]
    result, avg_roc_curves['CRC M+']  = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='CRC M+',nruns=args.runs,K_tag='colorectal_cancer')
    all_res.append(pd.DataFrame(result))
    
    
    ## 
    ## GAC M+
    sub=full[((full['Disease_status']=='gastric_cancer')  & (full['Metastasis_status']!='M0'))]
    result,  avg_roc_curves['GAC M+'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='GAC M+',nruns=args.runs,K_tag='gastric_cancer')
    all_res.append(pd.DataFrame(result))
    
    sub=full[((full['Disease_status']=='gastric_cancer')  & (full['Metastasis_status']!='M+'))]
    result,avg_roc_curves['GAC M0'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='GAC M0',nruns=args.runs,K_tag='early_gastric_cancer')
    all_res.append(pd.DataFrame(result))

    ## UVM
    sub=full[(full['Disease_status']=='uveal_melanoma_cancer')]
    result, avg_roc_curves['UVM M+'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='UVM M+',nruns=args.runs,K_tag='uveal_melanoma_cancer')
    all_res.append(pd.DataFrame(result))

    ## LC M+
    sub=full[(full['Disease_status']=='lung_cancer')]
    result,  avg_roc_curves['LC M+']= get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='LC M+',nruns=args.runs,K_tag='lung_cancer')
    all_res.append(pd.DataFrame(result))
    
    ## Early stage
    sub=full[((full['Stage']==1) | (full['Stage']==2))]
    result, avg_roc_curves['Early']= get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='Early',nruns=args.runs,K_tag='early')
    all_res.append(pd.DataFrame(result))
    
    
    ## Adv stage
    sub=full[(full['Stage']==3)]
    result,  avg_roc_curves['Adv.'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='Adv.',nruns=args.runs,K_tag='advanced')
    all_res.append(pd.DataFrame(result))
    
    ## Metastasic stage
    sub=full[(full['Stage']==4)]
    result,avg_roc_curves['Meta'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='Meta',nruns=args.runs,K_tag='meta')
    all_res.append(pd.DataFrame(result))
    
    ## all
    sub=full[(full['Disease_status']!='healthy')]
    result, avg_roc_curves['all cancer plasma'] = get_auc_roc_per_run(sub,healthy,'no_subsample',K_type='all cancer plasma',nruns=args.runs,K_tag='all')
    all_res.append(pd.DataFrame(result))

    
    
    full_res=pd.concat(all_res)
    full_res.to_csv(args.output_csv,index=False)
    data = []

    # Iterate over the dictionary
    for biological_class, categories in avg_roc_curves.items():
        for category, values in categories.items():
            if values is not None and len(values) == 2 and all(x is not None for x in values):
                # Create a separate row for each pair of fpr and tpr values
                for fpr, tpr in zip(values[0], values[1]):
                    # Append the data to the list
                    data.append({
                        "biological_class": biological_class,
                        "Blind_or_not": category,
                        "fpr": fpr,
                        "tpr": tpr
                    })
    df = pd.DataFrame(data)

    df.to_csv(args.output_csv.replace("_AUC.csv", "_tpr_fpr.csv"),index=False)
   
    # Plot ROC curves
    if args.output_plots:
        plot_roc_curves(full_res, avg_roc_curves, args.color_files, args.output_plots)

