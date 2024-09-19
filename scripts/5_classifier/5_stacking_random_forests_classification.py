#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse # arguments management
import pandas as pd # handle dataframes
import os # creating output directories
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import StackingClassifier
from sklearn.model_selection import train_test_split
from tqdm import tqdm # progress bar
import numpy as np

from sklearn.metrics import (accuracy_score,
                             f1_score,
                             auc,
                             roc_curve,
                             precision_recall_curve)
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from itertools import cycle # for colors
import math # for CI
from scipy import stats
# GLOBAL VARIABLES

FIG_DPI = 300

# FUNCTIONS

def compute_confidence_intervals(arr, confidence_level=0.95):
    """
    Compute the confidence intervals for a given array of values.

    Parameters:
    arr (list or numpy array): The input array of values.
    confidence_level (float, optional): The desired confidence level. Defaults to 0.95.

    Returns:
    tuple: A tuple containing the lower and upper bounds of the confidence interval.
    """
    # get the z score value to compute the confidence interval
    alpha=1-confidence_level 
    z_score=stats.norm.ppf(1-alpha/2)

    if len(arr) > 0 and type(arr[0]) is np.float64:
        mean=np.mean(arr)
        std_err=np.std(arr) / np.sqrt(len(arr))
        ci_lower_value=mean - z_score * std_err
        ci_upper_value=mean + z_score * std_err
        return ci_lower_value, ci_upper_value
    elif len(arr) > 0 and type(arr[0]) is np.ndarray:
        # if arr is a list of array then, the confidence interval is compute on column as
        # the data structure is supposed to be the following (nrun,specificity threshold).
        # Here, the goal is to compute at a defined threshold the confidence interval across runs
        ci_lowers = list()
        ci_uppers = list()
        arr=np.array(arr)
        for i in range(arr.shape[1]):
            sub_arr=arr[:,i]
            mean=np.mean(sub_arr)
            std_err=np.std(sub_arr) / np.sqrt(len(sub_arr))
            ci_lower_value=mean - z_score * std_err
            ci_upper_value=mean + z_score * std_err
            ci_lowers.append(ci_lower_value)
            ci_uppers.append(ci_upper_value)
        return (ci_lowers, ci_uppers)

def print_sensitivity_specificity_table(tprs, labels, output_dir):
    """
     Print a table of sensitivity and specificity values to a CSV file.
    
     Parameters:
     tprs (list): A list of true positive rates.
     labels (list): A list of labels corresponding to the true positive rates.
     output_dir (str): The directory where the output CSV file will be written.
     """
    with open(f'{output_dir}/sensitivity_specificity.csv', 'w') as f:
        f.write(f"biological_class,specificity,sensitivity,sensitivity_95CI_lower,sensitivity_95CI_upper\n")

        if len(labels) == 2:
            iter = [(0,labels[1])]
        else:
            iter = enumerate(labels)

        for nc, label in iter:
            mean_tpr = np.mean(tprs[nc], axis=0)
            ci_tpr = compute_confidence_intervals(tprs[nc])
            tprs_upper = np.minimum(ci_tpr[0], 1)
            tprs_lower = np.maximum(ci_tpr[1], 0)
            for percent in range(100):
                f.write(f"{label},{(100 - percent) / 100},{mean_tpr[percent]},{tprs_lower[percent]},{tprs_upper[percent]}\n")

def print_scores_to_csv(accuracies, f1s, tprs, recalls, mean_fpr, roc_aucs, pr_aucs, labels, output_dir):
    """
     Print prediction results to a CSV file.
    
     Parameters:
     d_predictions (dict): A dictionary of prediction results.
     labels (list): A list of labels corresponding to the prediction results.
     output_dir (str): The directory where the output CSV file will be written.
     """
    with open(f'{output_dir}/scores.csv', 'w') as f:
        header = "metric,score,lower_ci,upper_ci\n"
        f.write(header)
        ## accuracy
        accuracy = np.mean(accuracies)
        ci = compute_confidence_intervals(accuracies)
        print(accuracies)
        print(ci)
        f.write(f"accuracy,{accuracy},{ci[0]},{ci[1]}\n")
        ## f1 score
        f1 = np.mean(f1s)
        ci = compute_confidence_intervals(f1s)
        f.write(f"F1_score,{f1},{ci[0]},{ci[1]}\n")
        ## aucs
        if len(labels) == 2:
            iter = [(0,labels[1])]
        else:
            iter = enumerate(labels)
        for nc, label in iter:
            mean_tpr = np.mean(tprs[nc], axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = np.mean(roc_aucs[nc])
            ci = compute_confidence_intervals(roc_aucs[nc])
            f.write(f"ROC AUC {label},{mean_auc},{ci[0]},{ci[1]}\n")
            mean_recall = np.mean(recalls[nc], axis=0)
            mean_recall[-1] = 0.0
            mean_auc = np.mean(pr_aucs[nc])
            ci = compute_confidence_intervals(pr_aucs[nc])
            f.write(f"PR AUC {label},{mean_auc},{ci[0]},{ci[1]}\n")

def print_predictions_to_csv(d_predictions, labels, output_dir):
    """
    Print mean prediction results to a CSV file.
    
    Parameters:
    d_predictions (dict): A dictionary of prediction results.
    labels (list): A list of labels corresponding to the prediction results.
    output_dir (str): The directory where the output CSV file will be written.
    """
    # compute mean
    for sample in d_predictions:
        for j in range(len(labels)):
            if d_predictions[sample]['n'] != 0: # prevent division by zero
                d_predictions[sample][f'proba_{labels[j]}'] /= d_predictions[sample]['n']
                d_predictions[sample][f'pred_{labels[j]}'] /= d_predictions[sample]['n']
    # print
    with open(f'{output_dir}/predictions.csv', 'w') as f:
        ## header
        header="sample,biological_class,n,"
        for i in range(len(labels)):
            header += f"proba_{labels[i]},pred_{labels[i]},"
        header = header[:-1]+"\n"
        f.write(header)
        ## body
        for sample in d_predictions:
            line = f"{sample},{d_predictions[sample]['biological_class']},{d_predictions[sample]['n']},"
            for i in range(len(labels)):
                line += f"{d_predictions[sample][f'proba_{labels[i]}']},{d_predictions[sample][f'pred_{labels[i]}']},"
            line = line[:-1]+"\n"
            f.write(line)

def print_full_predictions_to_csv(f_predictions, output_dir):
    """
    Print the prediction results of all of the runs to a CSV file.

    Parameters:
    f_predictions (dict): A dictionary of full prediction results.
    output_dir (str): The directory where the output CSV file will be written.
    """  
    full_data=pd.DataFrame(f_predictions).T
    full_data.to_csv(f'{output_dir}/full_predictions.csv',index=False)

def print_blind_test_to_csv(blind_test, output_dir):
    """
    Print the blind test of all runs to a csv file

    Parameters:
    blind_test (dict) : A dictionary of cancer list seen and non seen in training set
    output_dir (str): The directory where the output CSV file will be written.
    """  
    full_data=pd.DataFrame(blind_test).T
    full_data.to_csv(f'{output_dir}/blind_test.csv',index=False)


def plot_roc_curve(tprs, mean_fpr, aucs, labels, output_dir, color_table_path):
    """
    Plot the ROC curve.
    
    Parameters:
    tprs (list): A list of true positive rates.
    mean_fpr (numpy array): The mean false positive rate.
    aucs (list): A list of AUC values.
    labels (list): A list of labels corresponding to the ROC curve.
    output_dir (str): The directory where the output plot will be written.
    color_table_path (str): The path to the color table CSV file.
    """
    fig, ax = plt.subplots(figsize=(7, 6.5), constrained_layout=True)
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    # check colors
    if color_table_path:
        color_table = pd.read_csv(color_table_path)
        colors = []
        for l in labels:
            if l in np.array(color_table.iloc[:,0]):
                col = color_table.loc[color_table.iloc[:,0] == l].iloc[:,1].item()
                colors.append(col)
            else:
                colors.append('NaN')
        if 'NaN' in colors:
            colors = cycle(mcolors.TABLEAU_COLORS)
    else:
        colors = cycle(mcolors.TABLEAU_COLORS)

    # iteration for n_class=2 or n_class>2
    zip_iter = list(zip(range(len(labels)), labels, colors))
    if len(labels) == 2:
        zip_iter = [(0,zip_iter[1][1],zip_iter[1][2])]

    for nc, label, color in zip_iter:

        mean_tpr = np.mean(tprs[nc], axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = np.mean(aucs[nc])
        ci_auc = compute_confidence_intervals(aucs[nc])
        auc_label = r"%s (AUC = %0.2f; 95%%CI = %0.2f–%0.2f)" % (label.replace('_', ' '), mean_auc, ci_auc[0], ci_auc[1])

        ax.plot(
            mean_fpr,
            mean_tpr,
            color=color,
            label=auc_label,
            lw=2,
            alpha=0.8,
        )

        ax.set(
            xlim=[-0.05, 1.05],
            ylim=[-0.05, 1.05],
            xlabel="1 - specificity",
            ylabel="Sensitivity",
            title=f"Mean ROC curve",
        )
        ax.axis("square")
        ax.legend(loc="lower right")

    plt.savefig(os.path.join(output_dir, "roc_curve.png"), dpi=FIG_DPI)

    for nc, label, color in zip_iter:
        ci_tpr = compute_confidence_intervals(tprs[nc])
        tprs_upper = np.minimum(ci_tpr[0], 1)
        tprs_lower = np.maximum(ci_tpr[1], 0)
        ax.fill_between(
            mean_fpr,
            tprs_lower,
            tprs_upper,
            color=color,
            alpha=0.2,
            label=None,
        )
    plt.savefig(os.path.join(output_dir, "roc_curve_ci.png"), dpi=FIG_DPI)

def plot_pr_curve(recalls, mean_fpr, aucs, labels, output_dir, color_table_path):
    """
    Plot the Presicion-Recall curve.

    Parameters:
    recalls (list): A list of recall values.
    mean_fpr (numpy array): The mean false positive rate.
    aucs (list): A list of AUC values.
    labels (list): A list of labels corresponding to the PR curve.
    output_dir (str): The directory where the output plot will be written.
    color_table_path (str): The path to the color table CSV file.
    """
    fig, ax = plt.subplots(figsize=(7, 6.5), constrained_layout=True)
    ax.plot([0, 1], [0,0], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    # check colors
    if color_table_path:
        color_table = pd.read_csv(color_table_path)
        colors = []
        for l in labels:
            if l in np.array(color_table.iloc[:,0]):
                col = color_table.loc[color_table.iloc[:,0] == l].iloc[:,1].item()
                colors.append(col)
            else:
                colors.append('NaN')
        if 'NaN' in colors:
            colors = cycle(mcolors.TABLEAU_COLORS)
    else:
        colors = cycle(mcolors.TABLEAU_COLORS)

    # iteration for n_class=2 or n_class>2
    zip_iter = list(zip(range(len(labels)), labels, colors))
    if len(labels) == 2:
        zip_iter = [(0,zip_iter[1][1],zip_iter[1][2])]

    for nc, label, color in zip_iter:

        mean_recall = np.mean(recalls[nc], axis=0)
        mean_recall[-1] = 0.0
        mean_auc = np.mean(aucs[nc])
        ci_auc = compute_confidence_intervals(aucs[nc])
        auc_label = r"%s (AUC = %0.2f; 95%%CI = %0.2f–%0.2f)" % (label.replace('_', ' '), mean_auc, ci_auc[0], ci_auc[1])

        ax.plot(
            mean_fpr,
            mean_recall,
            color=color,
            label=auc_label,
            lw=2,
            alpha=0.8,
        )

        ax.set(
            xlim=[-0.05, 1.05],
            ylim=[-0.05, 1.05],
            xlabel="Recall",
            ylabel="Precision",
            title=f"Mean Precision-Recall curve",
        )
        ax.axis("square")
        ax.legend(loc="lower right")

    plt.savefig(os.path.join(output_dir, "pr_curve.png"), dpi=FIG_DPI)

    for nc, label, color in zip_iter:
        ci_tpr = compute_confidence_intervals(recalls[nc])
        recalls_upper = np.minimum(ci_tpr[0], 1)
        recalls_lower = np.maximum(ci_tpr[1], 0)
        ax.fill_between(
            mean_fpr,
            recalls_lower,
            recalls_upper,
            color=color,
            alpha=0.2,
            label=None,
        )
    plt.savefig(os.path.join(output_dir, "pr_curve_ci.png"), dpi=FIG_DPI)

def run_classification(discovery_set_file, output_dir, runs, color_table_path,blind_size,minimumKsample2beInTrain):
    """
    Run the classification pipeline.

    Parameters:
    discovery_set_file (str): The path to the discovery set CSV file.
    output_dir (str): The directory where the output files will be written.
    runs (int): The number of runs to perform.
    color_table_path (str): The path to the color table CSV file.

    WARNING: Both data files must be CSVs with the control class (e.g., healthy plasma) before the case class (e.g., ovarian cancer plasma) to have the correct TN, FP, FN, TP order (see sklearn confusion matrix documentation).
    """

    # Load discovery set and extract information
    discovery_data = pd.read_csv(discovery_set_file)
    features = discovery_data.columns[2:]  # Extract feature columns
    y = discovery_data["biological_class"]  # Extract target variable
    _, labels = pd.factorize(y)  # Factorize target variable
    labelsHD_Cancer = ["healthy_plasma", "cancer_plasma"]  # Define labels for healthy and cancer plasma

    # Create output directory if it doesn't exist
    if output_dir:
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError:
                print(f"Creation of directory '{output_dir}' failed.")

    # Initialize output dictionaries
    data = discovery_data
    d_predictions = {}  # Dictionary to store prediction results
    f_predictions = {}  # Dictionary to store full prediction results
    for i, sample in enumerate(data["sample"]):
        d_predictions[sample] = {}
        d_predictions[sample]['biological_class'] = data.iloc[i].loc['biological_class']
        d_predictions[sample]['n'] = 0

        f_predictions[sample] = {}
        f_predictions[sample]['Sample_ID'] = sample

        for j in range(2):
            d_predictions[sample][f'proba_{labelsHD_Cancer[j]}'] = 0
            d_predictions[sample][f'pred_{labelsHD_Cancer[j]}'] = 0

    # Initialize variables to store performance metrics
    mean_fpr = np.linspace(0, 1, 100)
    tprs = {}
    recalls = {}
    roc_aucs = {}
    pr_aucs = {}
    accuracies = []
    f1s = []
    if blind_size > 0 :
        blind_test = {}
    
    
    tprs[0] = []
    recalls[0] = []
    roc_aucs[0] = []
    pr_aucs[0] = []
    
    cw =  None





    # Perform multiple runs
    for run_i in tqdm(range(0, runs)):
        # Split data into training and testing sets
      
       
        new_y = y
        _,labels=pd.factorize(new_y)
        X_train = discovery_data
        
        # Split healthy plasma data into training and testing sets
        y_train_hd, y_test_hd, X_train_hd, X_test_hd = train_test_split(new_y[new_y == "healthy_plasma"], X_train[new_y == "healthy_plasma"], test_size=0.33, random_state=None)
        
        # Split healthy plasma data into training sets for experts models and the stack
        y_train_expert_hd, y_train_stack_hd, X_train_expert_hd, X_train_stack_hd = train_test_split(y_train_hd, X_train_hd, test_size=0.5, random_state=None)
        
        # Initialize list to store classifiers
        classifiers = []
        
        # Initialize variables to store training stack cancer data
        X_stack_cancer_train = pd.DataFrame()
        y_stack_cancer_train = pd.Series()
        X_stack_cancer_test = pd.DataFrame()
        y_stack_cancer_test = pd.Series()


        # If minimumKsample2beInTrain > 0 : remove from training set every categories 
        # cancer that have less than minimumKsample2beInTrain value

        if minimumKsample2beInTrain > 0 :
            rejected_label=[]
            for label in labels :
                count=len(new_y[new_y==label])
                if count < minimumKsample2beInTrain :
                    print(f'{label} have been rejected from training set and moved to test set due to insufisent amount of sample. \n Nb sample={count} ; Minimum threshold {minimumKsample2beInTrain}')
                    labels=np.delete(labels,np.isin(labels,label))
                    rejected_label.append(label)


        # If blind size > 0 : remove X categories of cancers from training set

        if blind_size > 0 :
            drop_index=np.random.choice(labels[1:len(labels)+1],blind_size,replace=False)
            labels=np.delete(labels,np.isin(labels,drop_index))

            blind_test[run_i]={}
            blind_test[run_i]['run']=run_i

            ## because there two source of cancer cat removing (blind test or minimumKsample2beInTrain), here labs need to be recomputed, otherwise some data is lost
            _,all_lab=pd.factorize(y)
            seen_labs=np.delete(all_lab,np.isin(all_lab,labels,invert=True))
            blind_labs=np.delete(all_lab,np.isin(all_lab,labels))
            for label in seen_labs :
                blind_test[run_i][label]='seen'
            for drop_i in blind_labs :
                blind_test[run_i][drop_i]='blind'

            

        # Train expert classifiers for each cancer type
        for i in range(len(labels) - 1):
            expert = RandomForestClassifier(n_estimators=args.number_trees, criterion='gini', class_weight=cw)
            print(i)
        
            # Split cancer data into training and testing sets
            y_train_cancer, y_test_cancer, X_train_cancer, X_test_cancer = train_test_split(new_y[new_y == labels[i + 1]], X_train[new_y == labels[i + 1]], test_size=0.33, random_state=None)
        
            # Split cancer data into expert and stack sets
            y_train_expert_cancer, y_train_rest_cancer, X_train_expert_cancer, X_train_rest_cancer = train_test_split(y_train_cancer, X_train_cancer, test_size=0.5, random_state=None)
        
            # Train expert classifier
            X_train_i = pd.concat([X_train_expert_hd, X_train_expert_cancer])
            y_train_i = pd.concat([y_train_expert_hd, y_train_expert_cancer])
            y_train_i = (y_train_i == "healthy_plasma").astype(int)
            expert.fit(X_train_i[features], y_train_i)
        
            # Store expert classifier
            classifiers.append((labels[i + 1], expert))
        
            # Store cancer data for stack training
            X_stack_cancer_train = pd.concat([X_stack_cancer_train, X_train_rest_cancer])
            y_stack_cancer_train = pd.concat([y_stack_cancer_train, y_train_rest_cancer])
            X_stack_cancer_test = pd.concat([X_stack_cancer_test, X_test_cancer])
            y_stack_cancer_test = pd.concat([y_stack_cancer_test, y_test_cancer])
       
       # Consequently, every non seen in training K type is added to test set (by blind test or due to minimumKsample2beInTrain)
        if minimumKsample2beInTrain > 0 and len(rejected_label) > 0 :
            for label in rejected_label :
                X_stack_cancer_test=pd.concat([X_stack_cancer_test,X_train[new_y==label]])
                y_stack_cancer_test=pd.concat([y_stack_cancer_test,new_y[new_y==label]])
        if blind_size > 0 :
            for drop_i in drop_index :
                X_stack_cancer_test=pd.concat([X_stack_cancer_test,X_train[new_y==drop_i]])
                y_stack_cancer_test=pd.concat([y_stack_cancer_test,new_y[new_y==drop_i]])

        
        # Train stacked classifier
        stack = StackingClassifier(classifiers, RandomForestClassifier(n_estimators=args.number_trees, criterion='gini'), passthrough=True, stack_method="predict_proba", cv="prefit")
        X_train_stack = pd.concat([X_train_expert_hd, X_stack_cancer_train])
        y_train_stack = pd.concat([y_train_expert_hd, y_stack_cancer_train])
        y_train_stack = (y_train_stack == "healthy_plasma").astype(int)
        stack.fit(X_train_stack[features], y_train_stack)
        
        # Get prediction probabilities and predicted classes
        X_test = pd.concat([X_test_hd, X_stack_cancer_test])
        y_test = pd.concat([y_test_hd, y_stack_cancer_test])
        y_test = (y_test == "healthy_plasma").astype(int)
        prediction_proba = stack.predict_proba(X_test[features])
        y_pred = stack.predict(X_test[features])  # as integers
        y_pred_labeled = np.where(y_pred, 'healthy_plasma', 'cancer_plasma')  # as labels
                
        for i, pred in enumerate(prediction_proba):
            sample = X_test.iloc[i].loc['sample']
            
            d_predictions[sample][f'proba_{labelsHD_Cancer[1]}_run_{run_i}'] = pred[0]
            d_predictions[sample][f'proba_{labelsHD_Cancer[0]}_run_{run_i}'] = pred[1]
    
            
            d_predictions[sample][f'pred_{y_pred_labeled[i]}'] += 1
            d_predictions[sample]['n'] += 1
  
    
            f_predictions[sample][f'proba_{labelsHD_Cancer[1]}_run_{run_i}'] = pred[0]
            f_predictions[sample][f'proba_{labelsHD_Cancer[0]}_run_{run_i}'] = pred[1]
    

        accuracies.append(accuracy_score(y_test, y_pred))
        f1s.append(f1_score(y_test, y_pred, average = 'weighted'))
       
    
        fpr, tpr, _ = roc_curve(y_test, prediction_proba[:, 1], pos_label=1) # tpr = sensitivity = recall / fpr = specificity
        interp_tpr = np.interp(mean_fpr, fpr, tpr)
        interp_tpr[0] = 0.0
        tprs[0].append(interp_tpr)
        roc_aucs[0].append(auc(fpr, tpr))
    
        precision, recall, _ = precision_recall_curve(y_test, prediction_proba[:, 1], pos_label=1) # recall = tpr = sensitivity / precision = tp/(fp+tp)
        interp_recall = np.interp(mean_fpr, precision, recall)
        interp_recall[0] = 1.0
        recalls[0].append(interp_recall)
        pr_aucs[0].append(auc(recall,precision))
    
       # end of multiple runs

    # print classification scores
    #print_scores_to_csv(accuracies, f1s, tprs, recalls, mean_fpr, roc_aucs, pr_aucs, labelsHD_Cancer, output_dir)

    # print prediction results
    print_predictions_to_csv(d_predictions, labelsHD_Cancer, output_dir)

    # print full prediction results
    print_full_predictions_to_csv(f_predictions, output_dir)
    # print blind test
    if blind_size > 0 :
        print_blind_test_to_csv(blind_test,output_dir)
    # print sensitivity/specificity
    print_sensitivity_specificity_table(tprs,labelsHD_Cancer, output_dir)

    # plot "mean" ROC curve
    plot_roc_curve(tprs, mean_fpr, roc_aucs, labelsHD_Cancer, output_dir, color_table_path)

    # plot "mean" PR curve
    plot_pr_curve(recalls, mean_fpr, pr_aucs, labelsHD_Cancer, output_dir, color_table_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--discovery_set_file", type=str, default=None, help="input cg_methyl/haplotypes data file for discovery training")
    parser.add_argument("-r", "--runs", type=int, default=5000, help="number of runs")
    parser.add_argument("-t", "--number_trees", type=int, default=300, help="number of trees in the RF")
    parser.add_argument("-m", "--model_meta", type=str, default=None, help="Additionaly to cancer types expert add one or more model training on metastatic status. Binary: add a m0 vs HD and a m+ vs HD models. multiclass: add a multiclas HD vs M0 vs M+ model. None for no additional model ")
    parser.add_argument("-s", "--true_split_expert", type=bool, default=False, help="Do expert cancer types models learn on half or all of cancer. ")
    parser.add_argument("-o", "--output_dir", type=str, default=None, help="output directory")
    parser.add_argument("-c","--color_table_path", type=str, default=None, help="color table")
    parser.add_argument("-z","--blind_size",type=int,default=0,help="number of cancer types non seen in training set")
    parser.add_argument("-N","--minimumKsample2beInTrain",type=int,default=0,help="minimum amount of sample for a cancer type/stage to be accepted in train set.")


    args = parser.parse_args()

    if not args.discovery_set_file:
        print("a training set is required!")
        exit(1)




    if args.color_table_path == "None":
        args.color_table_path = None


    print(f"Classification has been run with the following arguments:\n- Training set: {args.discovery_set_file}\n- Number of runs: {args.runs}\n- Number of trees: {args.number_trees}\n- Ouput directory: {args.output_dir}\n- Color table: {args.color_table_path}\n")

    run_classification(discovery_set_file=args.discovery_set_file, output_dir=args.output_dir, runs=args.runs, color_table_path=args.color_table_path, blind_size=args.blind_size,minimumKsample2beInTrain=args.minimumKsample2beInTrain)

