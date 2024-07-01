import argparse # arguments management
import pandas as pd # handle dataframes
import os # creating output directories
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from tqdm import tqdm # progress bar
import numpy as np
from imblearn.over_sampling import SMOTE # oversampling
from imblearn.under_sampling import RandomUnderSampler # undersampling
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

TEST_DATASET_SIZE = .4
FIG_DPI = 300

# FUNCTIONS

def compute_confidence_intervals(arr, confidence_level=0.95):
    """Expects a list of either numpy arrays or numpy float64."""
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
    with open(f'{output_dir}/scores.csv', 'w') as f:
        header = "metric,score,lower_ci,upper_ci\n"
        f.write(header)
        ## accuracy
        accuracy = np.mean(accuracies)
        ci = compute_confidence_intervals(accuracies)
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
    # compute mean
    for sample in d_predictions:
        for j in range(len(labels.unique())):
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
    # print
    full_data=pd.DataFrame(f_predictions).T
    full_data.to_csv(f'{output_dir}/full_predictions.csv',index=False)

def print_importances_to_csv(d_features_imp, output_dir):
    # compute mean
    for feature in d_features_imp.keys():
        d_features_imp[feature]["mean"] = np.mean(d_features_imp[feature]["scores"])
        d_features_imp[feature]["median"] = np.median(d_features_imp[feature]["scores"])
        d_features_imp[feature]["stdev"] = np.std(d_features_imp[feature]["scores"])
        d_features_imp[feature]["total_score"] = sum(d_features_imp[feature]["scores"])
        d_features_imp[feature]["mean_rank"] = np.mean(d_features_imp[feature]["ranks"])
        del d_features_imp[feature]["scores"]
        del d_features_imp[feature]["ranks"]
    # print
    with open(f'{output_dir}/features_importance.csv', 'w') as f:
        ## header
        header="feature,mean,median,stdev,total_score,mean_rank\n"
        f.write(header)
        ## body
        for feature in d_features_imp:
            line = f"{feature},{d_features_imp[feature]['mean']},{d_features_imp[feature]['median']},{d_features_imp[feature]['stdev']},{d_features_imp[feature]['total_score']},{d_features_imp[feature]['mean_rank']}\n"
            f.write(line)

def plot_roc_curve(tprs, mean_fpr, aucs, labels, output_dir, color_table_path):

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

def run_classification(discovery_set_file, validation_set_file, output_dir, runs, balance_method, color_table_path):
    """
    data_file must be a csv with control class (e.g. healthy plasma) before
    case class (e.g. ovarian cancer plasma) in order to have correct tn, fp,
    fn, tp order (see sklearn confusion matrix documentation).
    """

    # Load discovery set and extract infos
    discovery_data = pd.read_csv(discovery_set_file)
    features = discovery_data.columns[2:]
    y, labels = pd.factorize(discovery_data['biological_class'])
    n_class = len(discovery_data["biological_class"].unique())

    # Load validation set (if any)
    if validation_set_file:
        validation_data = pd.read_csv(validation_set_file)
        ## check consistency of the features
        if not np.array_equal(validation_data.columns[2:], features):
            print("Features aren't identical between train and test data.")

    # create output directory
    if output_dir:
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError:
                print(f"Creation of directory '{output_dir}' failed.")

    # Initialize output dictionaries
    data = validation_data if validation_set_file else discovery_data
    ## prediction probabilities and classes
    d_predictions = {}
    f_predictions = {} # full prediction
    for i,sample in enumerate(data["sample"]):
        d_predictions[sample] = {}
        d_predictions[sample]['biological_class'] = data.iloc[i].loc['biological_class']
        d_predictions[sample]['n'] = 0

        f_predictions[sample] = {}
        f_predictions[sample]['Sample_ID']=sample

        for j in range(n_class):
            d_predictions[sample][f'proba_{labels[j]}'] = 0
            d_predictions[sample][f'pred_{labels[j]}'] = 0
    ## features importance
    d_features_imp = {}
    for feature in features:
        d_features_imp[feature] = {}
        d_features_imp[feature]["scores"] = []
        d_features_imp[feature]["ranks"] = []
    ## roc and pr curves info
    mean_fpr = np.linspace(0, 1, 100)
    tprs = {}
    recalls = {}
    roc_aucs = {}
    pr_aucs = {}
    accuracies = []
    f1s = []
    if n_class == 2:
        tprs[0] = []
        recalls[0] = []
        roc_aucs[0] = []
        pr_aucs[0] = []
    else:
        for i in range(n_class):
            tprs[i] = []
            recalls[i] = []
            roc_aucs[i] = []
            pr_aucs[i] = []

    # Create a random forest classifier
    cw = "balanced_subsample" if balance_method=="class_weight" else None
    clf = RandomForestClassifier(n_estimators=300, criterion='gini',class_weight=cw)

    # Multiple runs

    for run_i in tqdm(range(0, runs)):

        if not validation_set_file:
            if balance_method=="undersampling":
                new_data, new_y = RandomUnderSampler(random_state=run_i).fit_resample(data, y)
            else:
                new_data = data
                new_y = y
            # split train/test sets
            X_train, X_test, y_train, y_test = train_test_split(new_data, new_y, test_size=TEST_DATASET_SIZE, stratify=new_y)
        else:
            if balance_method=="undersampling":
                new_data, new_y = RandomUnderSampler(random_state=run_i).fit_resample(discovery_data, y)
            else:
                new_data = discovery_data
                new_y = y
            X_train = new_data
            X_test = validation_data
            y_train = new_y
            y_test, _ = pd.factorize(validation_data['biological_class'])

        # training
        clf.fit(X_train[features], y_train)

        # get prediction probabilities and predicted classes
        prediction_proba = clf.predict_proba(X_test[features])
        y_pred = clf.predict(X_test[features]) # as integers
        y_pred_labeled = labels[y_pred] # as labels
        for i, pred in enumerate(prediction_proba):
            sample = X_test.iloc[i].loc['sample']
            for j in range(n_class):
                d_predictions[sample][f'proba_{labels[j]}'] += pred[j]
            d_predictions[sample][f'pred_{y_pred_labeled[i]}'] += 1
            d_predictions[sample]['n'] += 1

            f_predictions[sample][f'proba_{labels[0]}_run_{run_i}'] = pred[0]
            f_predictions[sample][f'proba_{labels[1]}_run_{run_i}'] = pred[1]

        # get features importance
        feature_importances = clf.feature_importances_
        ranks = pd.factorize(-feature_importances, sort=True)[0] + 1
        for i in range(len(features)):
            d_features_imp[features[i]]["scores"].append(feature_importances[i])
            d_features_imp[features[i]]["ranks"].append(ranks[i])

        # get "single" roc curve infos and precision recall
        accuracies.append(accuracy_score(y_test, y_pred))
        f1s.append(f1_score(y_test, y_pred, average = 'weighted'))
        if n_class == 2:

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

        else:

            for i in range(n_class):

                fpr, tpr, _ = roc_curve(y_test, prediction_proba[:, i], pos_label=i)
                interp_tpr = np.interp(mean_fpr, fpr, tpr)
                interp_tpr[0] = 0.0
                tprs[i].append(interp_tpr)
                roc_aucs[i].append(auc(fpr, tpr))

                precision, recall, _ = precision_recall_curve(y_test, prediction_proba[:, i], pos_label=i)
                interp_recall = np.interp(mean_fpr, precision, recall)
                interp_recall[0] = 1.0
                recalls[i].append(interp_recall)
                pr_aucs[i].append(auc(recall,precision))

    # end of multiple runs

    # print classification scores
    print_scores_to_csv(accuracies, f1s, tprs, recalls, mean_fpr, roc_aucs, pr_aucs, labels, output_dir)

    # print prediction results
    print_predictions_to_csv(d_predictions, labels, output_dir)

    # print full prediction results
    print_full_predictions_to_csv(f_predictions, output_dir)

    # print features importance results
    print_importances_to_csv(d_features_imp, output_dir)

    # print sensitivity/specificity
    print_sensitivity_specificity_table(tprs,labels, output_dir)

    # plot "mean" ROC curve
    plot_roc_curve(tprs, mean_fpr, roc_aucs, labels, output_dir, color_table_path)

    # plot "mean" PR curve
    plot_pr_curve(recalls, mean_fpr, pr_aucs, labels, output_dir, color_table_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d","--discovery_set_file", type=str, default=None, help="input cg_methyl/haplotypes data file for discovery training")
    parser.add_argument("-v","--validation_set_file", type=str, default=None, help="input cg_methyl/haplotypes data file for validation")
    parser.add_argument("-r", "--runs", type=int, default=5000, help="number of runs")
    parser.add_argument("-b", "--balance_method", type=str, default=None, help="balancing method for imbalanced classes [undersampling, oversampling, class_weight, None]")
    parser.add_argument("-o", "--output_dir", type=str, default=None, help="output directory")
    parser.add_argument("-c","--color_table_path", type=str, default=None, help="color table")
    args = parser.parse_args()

    if not args.discovery_set_file:
        print("Testing set is optional but a training set is required!")
        exit(1)

    if args.validation_set_file == "None":
        args.validation_set_file = None

    if args.color_table_path == "None":
        args.color_table_path = None

    if args.balance_method != "undersampling" and args.balance_method != "oversampling" and args.balance_method != "class_weight" and args.balance_method != "None":
        print("Invalid balance method. Choose between undersampling, oversampling, class_weight, or None.")
        exit(1)

    print(f"Classification has been run with the following arguments:\n- Training set: {args.discovery_set_file}\n- Testing set: {args.validation_set_file}\n- Number of runs: {args.runs}\n- Balance method: {args.balance_method}\n- Ouput directory: {args.output_dir}\n- Color table: {args.color_table_path}\n")

    run_classification(discovery_set_file=args.discovery_set_file, validation_set_file=args.validation_set_file, output_dir=args.output_dir, runs=args.runs, balance_method=args.balance_method, color_table_path=args.color_table_path)
