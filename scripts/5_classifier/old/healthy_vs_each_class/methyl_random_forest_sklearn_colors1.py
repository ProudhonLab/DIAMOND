#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Methylation analysis using random forest.

@author: Marc Michel
@email: marc.michel@curie.fr
@project: https://github.com/michel-m
"""
import argparse
import itertools
import logging
import math
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys

from itertools import cycle
from math import sqrt
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (auc,
                             confusion_matrix,
                             roc_curve,
                             RocCurveDisplay)
from sklearn.model_selection import train_test_split
from statistics import mean, median, stdev
from tabulate import tabulate
from tqdm import tqdm


# Switch to 'Agg' for headless graphics generation (such as on a cluster)
mpl.use('Agg')
mpl.rcParams['figure.dpi'] = 80
if plt.get_backend() == 'Qt5Agg':
    from matplotlib.backends.qt_compat import QtWidgets
    qApp = QtWidgets.QApplication(sys.argv)
    plt.matplotlib.rcParams['figure.dpi'] = qApp.desktop().physicalDpiX()

TEST_DATASET_SIZE = .4
TOP_FEATURES_PERC = 10
FIG_DPI = 300
SAMPLES_INFO_FILE = '/home/genouest/cnrs_umr6074/kdasilva/proudhon_lab/psl1_meth/scripts/5_classifier/old/healthy_vs_each_class/sample_biological-class_annotation.csv'

def load_samples_info(samples_info_file_path):
    return pd.read_csv(samples_info_file_path, index_col='sample')


def compute_confidence_intervals(arr, alpha=0.95):
    """Expects a list of either numpy arrays or numpy float64."""
    if len(arr) > 0 and type(arr[0]) is np.float64:
        arr.sort()
        ci_lower_index = math.ceil(len(arr) * (1 - alpha)) - 1
        ci_upper_index = math.floor(len(arr) * alpha) - 1
        return arr[ci_lower_index], arr[ci_upper_index]
    elif len(arr) > 0 and type(arr[0]) is np.ndarray:
        ci_lowers = list()
        ci_uppers = list()
        for i in range(len(arr[0])):
            tprs_at_each_tick = list()
            ci_lower_index = math.ceil(len(arr) * (1 - alpha)) - 1
            ci_upper_index = math.floor(len(arr) * alpha) - 1
            for sub_arr in arr:
                tprs_at_each_tick.append(sub_arr[i])
            tprs_at_each_tick.sort()
            ci_lowers.append(tprs_at_each_tick[ci_lower_index])
            ci_uppers.append(tprs_at_each_tick[ci_upper_index])
        return (ci_lowers, ci_uppers)


def plot_roc_curve_multiclass_crossval(ax, tprs, mean_fpr, aucs, labels, ci_type="ci"):
    plot_ci = True if len(labels) < 3 else False

    ax.plot([0, 1], [0, 1], linestyle="--", lw=2,
            color="r", label="Chance", alpha=0.8)

    colors = cycle(mcolors.TABLEAU_COLORS)

    for nc, label, color in zip(range(len(labels)), labels, colors):
        for i in range(len(tprs[nc])):
            tprs[nc][i][0] = 0.0
        mean_tpr = np.mean(tprs[nc], axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)

        ## temporary block for colors
        if label=="healthy_plasma":
            color="grey"
        if label=="gastric_cancer_plasma":
            color="#09BE68"
        if label=="ovarian_cancer_plasma":
            color="#A52A2B"
        if label=="colorectal_cancer_plasma":
            color="#CD9401"
        if label=="breast_cancer_plasma":
            color="#004096"
        if label=="uveal_melanoma_cancer_plasma":
            color="#FF61CC"
        if label=="lung_cancer_plasma":
            color="#C67CFE"
        if label=="early_breast_cancer_plasma":
            color="#00A9FF"

        if ci_type == "stdev":
            ci_auc = np.std(aucs[nc], ddof=1)
            auc_label = r"%s (AUC = %0.2f $\pm$ %0.2f)" % (label.replace('_', ' '), mean_auc, ci_auc)
        elif ci_type == "ci":
            ci_auc = compute_confidence_intervals(aucs[nc])
            auc_label = r"%s (AUC = %0.2f; 95%%CI = %0.2f–%0.2f)" % (
                label.replace('_', ' '), mean_auc, ci_auc[0], ci_auc[1])
        ax.plot(
            mean_fpr,
            mean_tpr,
            color=color,
            label=auc_label,
            lw=2,
            alpha=0.8,
        )

        if plot_ci is True:
            if ci_type == "stdev":
                ci_tpr = np.std(tprs[nc], axis=0, ddof=1)
                tprs_upper = np.minimum(mean_tpr + ci_tpr, 1)
                tprs_lower = np.maximum(mean_tpr - ci_tpr, 0)
                tprs_label = r"$\pm$ std. dev."
            elif ci_type == "ci":
                ci_tpr = compute_confidence_intervals(tprs[nc])
                tprs_upper = np.minimum(ci_tpr[0], 1)
                tprs_lower = np.maximum(ci_tpr[1], 0)
                tprs_label = r"$\pm$ 95%CI"
            ax.fill_between(
                mean_fpr,
                tprs_lower,
                tprs_upper,
                # color="grey",
                color=color,
                alpha=0.2,
                # label=tprs_label,
                label=None,
            )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title="Receiver operating characteristic",
    )
    ax.set_xlabel('1 - specificity')
    ax.set_ylabel('Sensitivity')

    ax.legend(loc="lower right")


def plot_roc_curve_multiclass(all_class_fprs, all_class_tprs, labels, ci_type='ci'):
    plot_specificity_dots = False
    plot_stdev = False

    f, ax = plt.subplots(figsize=(7, 6.5), constrained_layout=True)

    # Plot the chance line.
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
            alpha=.8)  # , label='Chance')

    colors = cycle(mcolors.TABLEAU_COLORS)

    for nc, label, color in zip(range(len(labels)), labels, colors):
        fprs = all_class_fprs[nc]
        tprs = all_class_tprs[nc]

        ## temporary block for colors
        if label=="healthy_plasma":
            color="grey"
        if label=="gastric_cancer_plasma":
            color="#09BE68"
        if label=="ovarian_cancer_plasma":
            color="#A52A2B"
        if label=="colorectal_cancer_plasma":
            color="#CD9401"
        if label=="breast_cancer_plasma":
            color="#004096"
        if label=="uveal_melanoma_cancer_plasma":
            color="#FF61CC"
        if label=="lung_cancer_plasma":
            color="#C67CFE"
        if label=="early_breast_cancer_plasma":
            color="#00A9FF"

        # Initialize useful lists + the plot axes.
        tprs_interp = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        # Plot ROC for each K-Fold + compute AUC scores.
        for i, (fpr, tpr) in enumerate(zip(fprs, tprs)):
            tprs_interp.append(np.interp(mean_fpr, fpr, tpr))
            tprs_interp[-1][0] = 0.0
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
            # ax.plot(fpr, tpr, lw=1, alpha=0.3,
            #         label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

        # Plot the mean ROC.
        mean_tpr = np.mean(tprs_interp, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        if ci_type == 'stdev':
            ci_auc = np.std(aucs, ddof=1)
            auc_label = r"%s (AUC = %0.2f $\pm$ %0.2f)" % (label.replace('_', ' '), mean_auc, ci_auc)
        elif ci_type == 'ci':
            ci_auc = compute_confidence_intervals(aucs)
            auc_label = r"%s (AUC = %0.2f; 95%%CI = %0.2f–%0.2f)" % (
                label.replace('_', ' '), mean_auc, ci_auc[0], ci_auc[1])

        ax.plot(mean_fpr, mean_tpr, color=color,
                label=auc_label,
                lw=2, alpha=.8)

        if plot_stdev is True:
            # Plot the standard deviation around the mean ROC.
            std_tpr = np.std(tprs_interp, axis=0, ddof=1)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey',
                            alpha=.2, label=r'$\pm$ std. dev.')

        ax.set_xlim([-0.05, 1.05])
        ax.set_ylim([-0.05, 1.05])
        ax.set_xlabel('1 - specificity')
        ax.set_ylabel('Sensitivity')
        ax.set_title('Receiver operating characteristic')
        ax.legend(loc="lower right")

        if plot_specificity_dots is True:
            specificity_texts = list()
            for specificity_perc in [94, 98, 99]:
                if plot_stdev is True:
                    std_fmt = r"($\pm$ " + f"{std_tpr[100 - specificity_perc] * 100:.1f}" + ")"
                    specificity_texts.append(f"{mean_tpr[100 - specificity_perc] * 100:.1f}% {std_fmt} sensitivity at {specificity_perc}% specificity")
                else:
                    specificity_texts.append(f"{mean_tpr[100 - specificity_perc] * 100:.1f}% sensitivity at {specificity_perc}% specificity")
                ax.plot(1 - (specificity_perc / 100), mean_tpr[100 - specificity_perc], 'bo')
            # plt.gcf().text(0.05, 0.01, "\n".join(specificity_texts), fontsize=11)
            ax.annotate("\n".join(specificity_texts), xy=(100, 20),
                        xycoords='figure pixels', fontsize=9)
            # , horizontalalignment='left', verticalalignment='top')
            # plt.subplots_adjust(bottom=0.05 * len(specificity_perc))
            # ax.text(0, 0, "\n".join(specificity_texts))

        f.set_constrained_layout_pads(w_pad=.15, h_pad=.22)


def sensitivity_specificity_table(tprs, labels):
    sensitivity_specificity_output = list()
    sensitivity_specificity_output.append(
        "biological_class,specificity,sensitivity,sensitivity_95CI_lower,sensitivity_95CI_upper")

    for nc, label in enumerate(labels):
        mean_tpr = np.mean(tprs[nc], axis=0)
        ci_tpr = compute_confidence_intervals(tprs[nc])
        tprs_upper = np.minimum(ci_tpr[0], 1)
        tprs_lower = np.maximum(ci_tpr[1], 0)

        for percent in range(100):
            sensitivity_specificity_output.append(
                f"{label},{(100 - percent) / 100},{mean_tpr[percent]},{tprs_lower[percent]},{tprs_upper[percent]}")

    return sensitivity_specificity_output


def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          ci=None,
                          ci_type="ci",
                          n_runs=0,
                          title='Confusion matrix',
                          cmap=plt.cm.Blues):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        if ci is not None:
            cm = cm.astype('float') / n_runs
        else:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    classes = [c.replace('_', ' ') for c in classes]

    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, ha="right",
               rotation_mode="anchor")
    plt.yticks(tick_marks, classes)

    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    if ci is None:
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            plt.text(j, i, format(cm[i, j], fmt),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")
    else:
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            if ci_type == "ci":
                label = f"{cm[i, j]:.2f}\n(95%CI\n{ci[i, j][0]:.2f}–{ci[i, j][1]:.2f})"
            elif ci_type == "stdev":
                label = f"{cm[i, j]:.2f} ($\\pm$ {ci[i, j]:.2f})"
            plt.text(j, i, label,
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")

    plt.ylabel('True label')
    plt.xlabel('Predicted label')


def output_prediction_scores(predictions_dict, labels, output_dir):
    """
    TODO: add 95% CI computation to output file
    """
    info = load_samples_info(SAMPLES_INFO_FILE)
    out = list()
    for k in predictions_dict.keys():
        for i in range(len(labels)):
            predictions_dict[k][f'pred_{labels[i]}'] /= predictions_dict[k]['n']
    out.append("sample,biological_class,")
    for i in range(len(labels)):
        out.append(f"pred_{labels[i]},")
    out.append("n,maf\n")
    for k in sorted(predictions_dict.keys()):
        bio_class = predictions_dict[k]['biological_class']
        n = predictions_dict[k]['n']
        try:
            maf = "NA" if np.isnan(info.loc[k]['maf']) else info.loc[k]['maf']
        except KeyError:
            maf = "NA"
        out.append(f"{k},{bio_class},")
        for i in range(len(labels)):
            pred_label = predictions_dict[k][f'pred_{labels[i]}']
            out.append(f"{pred_label:.3f},")
        out.append(f"{n},{maf}\n")

    if output_dir is None:
        print("".join(out))
    else:
        with open(os.path.join(output_dir, 'prediction_scores.csv'),
                  mode='w') as prediction_scores_output_file:
            prediction_scores_output_file.write("".join(out))


def run_classification(data_file, output_dir, runs, subsampling, plot_suffix, seed=None):
    """
    data_file must be a csv with control class (e.g. healthy plasma) before
    case class (e.g. ovarian cancer plasma) in order to have correct tn, fp,
    fn, tp order (see sklearn confusion matrix documentation).
    """
    # np.random.seed(seed)

    data = pd.read_csv(data_file)

    # Create a list of the feature column's names
    features = data.columns[2:]

    n_top_features = round(len(features) * TOP_FEATURES_PERC / 100)

    if output_dir is not None:
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError:
                print(f"Creation of directory '{output_dir}' failed.")

    # Create a random forest classifier
    clf = RandomForestClassifier(n_estimators=300, criterion='gini')

    misclassified_samples = dict()
    features_importance = dict()
    top_features = dict()
    feature_scores = dict()
    cm = None
    cm_list = list()
    predictions_dict = dict()

    unique_biological_classes = data.biological_class.unique()

    n_class = len(unique_biological_classes)
    fprs = dict()
    tprs = dict()
    for i in range(n_class):
        fprs[i] = list()
        tprs[i] = list()

    cross_val_tprs = dict()
    cross_val_aucs = dict()
    for i in range(n_class):
        cross_val_tprs[i] = list()
        cross_val_aucs[i] = list()
    cross_val_mean_fpr = np.linspace(0, 1, 100)
    cross_val_fig, cross_val_ax = plt.subplots(
        figsize=(7, 6.5), constrained_layout=True)

    samples_info = load_samples_info(SAMPLES_INFO_FILE)
    samples_info = samples_info[samples_info.index.isin(data['sample'])]
    samples_info_default_index = samples_info.copy().reset_index()
    data_true_bc = pd.merge(data['sample'], samples_info_default_index[['sample', 'biological_class']], how='left', on='sample')
    real_labels = pd.DataFrame(data={'label': data_true_bc['biological_class'].unique()})
    fprs = dict()
    tprs = dict()
    for i in range(len(real_labels['label'])):
        fprs[i] = list()
        tprs[i] = list()
        cross_val_tprs[i] = list()
        cross_val_aucs[i] = list()

    for run_n in tqdm(range(0, runs)):
        if subsampling is True:
            
            biological_classes_count = list()
            for biological_class in unique_biological_classes:
                biological_classes_count.append(len(data[data['biological_class'] == biological_class]))
                
            n_samples_min = min(biological_classes_count)
            new_data = data[data['biological_class'] == unique_biological_classes[0]].sample(n_samples_min)
            for i in range(1, len(unique_biological_classes)):
                new_data = pd.concat([new_data,data[data['biological_class'] == unique_biological_classes[i]].sample(n_samples_min)])
        else:
            new_data = data
        y, labels = pd.factorize(new_data['biological_class'])

        X_train, X_test, y_train, y_test = train_test_split(new_data, y, test_size=TEST_DATASET_SIZE, stratify=y)

        clf.fit(X_train[features], y_train)

        predictions = clf.predict_proba(X_test[features])
        for i, pred in enumerate(predictions):
            pred_sample = X_test.iloc[i].loc['sample']
            if pred_sample not in predictions_dict.keys():
                predictions_dict[pred_sample] = dict()
                predictions_dict[pred_sample]['biological_class'] = X_test.iloc[i].loc['biological_class']
                for j in range(n_class):
                    predictions_dict[pred_sample][f'pred_{labels[j]}'] = 0
                predictions_dict[pred_sample]['n'] = 0
            for j in range(n_class):
                predictions_dict[pred_sample][f'pred_{labels[j]}'] += pred[j]
            predictions_dict[pred_sample]['n'] += 1

        y_pred = clf.predict(X_test[features])

        y_pred_labeled = labels[y_pred]

        feature_importances_df = pd.DataFrame({'features': features, 'importance': clf.feature_importances_})

        for feature, importance in feature_importances_df.values:
            if feature not in features_importance.keys():
                features_importance[feature] = list()
            features_importance[feature].append(importance)

        # Aggregate top features
        feature_importances_df_sorted = feature_importances_df.sort_values('importance', ascending=False)
        for feature in feature_importances_df_sorted.iloc[0:n_top_features]['features']:
            if feature not in top_features.keys():
                top_features[feature] = 0
            top_features[feature] += 1

        feature_importances_df_sorted = feature_importances_df.sort_values('importance', ascending=True)
        for i in range(0, len(feature_importances_df_sorted) - 1):
            feature = feature_importances_df_sorted.iloc[i]['features']
            if feature not in feature_scores.keys():
                feature_scores[feature] = 0
            feature_scores[feature] += i + 1

        # Build misclassified items dictionary
        misclassified = X_test['biological_class'] == y_pred_labeled
        for i in misclassified.index:
            if X_test['sample'][i] not in misclassified_samples.keys():
                misclassified_samples[X_test['sample'][i]] = 0
            if misclassified[i] == False:  # noqa
                misclassified_samples[X_test['sample'][i]] += 1

        if n_class > 2:
            for i in range(n_class):
                fpr_RF, tpr_RF, _ = roc_curve(
                    y_test, predictions[:, i], pos_label=i)
                fprs[i].append(fpr_RF)
                tprs[i].append(tpr_RF)

                logging.disable(logging.WARNING)
                cross_val_viz = RocCurveDisplay.from_predictions(y_test, predictions[:, i], pos_label=i, label=None, alpha=0.0, lw=1, ax=cross_val_ax)
                logging.disable(logging.NOTSET)
                cross_val_interp_tpr = np.interp(cross_val_mean_fpr, cross_val_viz.fpr, cross_val_viz.tpr)
                # cross_val_interp_tpr[0] = 0.0
                cross_val_tprs[i].append(cross_val_interp_tpr)
                cross_val_aucs[i].append(cross_val_viz.roc_auc)
        elif samples_info['biological_class'].unique().size > 2:
            # associate true biological classes to X_test
            X_test_true_bc = pd.merge(X_test['sample'], samples_info_default_index[['sample', 'biological_class']], how='left', on='sample')
            y_test_real_bc = y_test.copy()
            for i in range(len(y_test)):
                y_test_real_bc[i] = real_labels[real_labels['label'] == X_test_true_bc.iloc[i]['biological_class']].index.item()
            y_pred_real_bc = y_pred.copy()
            for i in range(len(y_pred)):
                if y_pred[i] == 1 and y_test_real_bc[i] != 0:
                    y_pred_real_bc[i] = y_test_real_bc[i]
            for i in range(len(real_labels['label'])):
                fpr_RF, tpr_RF, _ = roc_curve(
                    y_test_real_bc, y_pred_real_bc, pos_label=i)
                if np.isnan(fpr_RF).any() or np.isnan(tpr_RF).any():
                    continue
                fprs[i].append(fpr_RF)
                tprs[i].append(tpr_RF)

                logging.disable(logging.WARNING)
                cross_val_viz = RocCurveDisplay.from_predictions(
                    y_test_real_bc, y_pred_real_bc, pos_label=i,
                    label=None,
                    alpha=0.0,
                    lw=1,
                    ax=cross_val_ax)
                logging.disable(logging.NOTSET)
                cross_val_interp_tpr = np.interp(
                    cross_val_mean_fpr, cross_val_viz.fpr, cross_val_viz.tpr)
                # cross_val_interp_tpr[0] = 0.0
                cross_val_tprs[i].append(cross_val_interp_tpr)
                cross_val_aucs[i].append(cross_val_viz.roc_auc)
        else:
            fpr_RF, tpr_RF, _ = roc_curve(
                y_test, predictions[:, 1], pos_label=1)  # /!\ control class MUST come first in the data file
            fprs[0].append(fpr_RF)
            tprs[0].append(tpr_RF)

            logging.disable(logging.WARNING)
            cross_val_viz = RocCurveDisplay.from_predictions(
                y_test, predictions[:, 1], pos_label=1,
                label=None,
                alpha=0.0,
                lw=1,
                ax=cross_val_ax)
            logging.disable(logging.NOTSET)
            cross_val_interp_tpr = np.interp(
                cross_val_mean_fpr, cross_val_viz.fpr, cross_val_viz.tpr)
            # cross_val_interp_tpr[0] = 0.0
            cross_val_tprs[0].append(cross_val_interp_tpr)
            cross_val_aucs[0].append(cross_val_viz.roc_auc)

        if cm is None:
            cm = confusion_matrix(X_test['biological_class'],
                                  y_pred_labeled, labels=labels)
        else:
            cm += confusion_matrix(X_test['biological_class'],
                                   y_pred_labeled, labels=labels)
        cm_list.append(confusion_matrix(X_test['biological_class'],
                                        y_pred_labeled, labels=labels))

    print("Generating outputs...")

    roc_labels = labels if n_class > 2 else [labels[1]]

    # Plot ROC curve crossval
    cross_val_fig.set_constrained_layout_pads(w_pad=.15, h_pad=.26)
    plot_roc_curve_multiclass_crossval(cross_val_ax, cross_val_tprs,
                                       cross_val_mean_fpr, cross_val_aucs,
                                       roc_labels, ci_type="ci")
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, f"roc_multiclass_crossval_ci{plot_suffix}.png"),
                    dpi=FIG_DPI)

    sensitivity_specificity_output = sensitivity_specificity_table(
        cross_val_tprs, roc_labels)
    if output_dir is None:
        print('\n'.join(sensitivity_specificity_output))
    else:
        with open(os.path.join(output_dir, 'sensitivity_specificity.csv'),
                  mode='w') as sensitivity_specificity_output_file:
            sensitivity_specificity_output_file.write(
                '\n'.join(sensitivity_specificity_output))

    # Output ROC multiclass
    plot_roc_curve_multiclass(fprs, tprs, roc_labels)
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir,
                                 f"roc_multiclass{plot_suffix}.png"),
                    dpi=FIG_DPI)

    # Output prediction scores
    output_prediction_scores(predictions_dict, labels, output_dir)

    # Output misclassified samples
    s = [(k, misclassified_samples[k]) for k in
         sorted(misclassified_samples, key=misclassified_samples.get,
                reverse=True)]
    misclassified_samples_output = list()
    misclassified_samples_output.append(
        "sample,biological_class,misclassified_n")
    for k, v in s:
        try:
            biological_class = samples_info.loc[k]['biological_class']
        except KeyError:
            biological_class = "unknown_biological_class"
        misclassified_samples_output.append(f"{k},{biological_class},{v}")
    if output_dir is None:
        print('\n'.join(misclassified_samples_output))
    else:
        with open(os.path.join(output_dir, 'misclassified_samples.csv'),
                  mode='w') as misclassified_samples_output_file:
            misclassified_samples_output_file.write(
                '\n'.join(misclassified_samples_output))

    np.set_printoptions(precision=2)

    # Plot non-normalized confusion matrix
    plt.figure(constrained_layout=True)
    plot_confusion_matrix(cm, classes=labels,
                          title='Confusion matrix, without normalization')
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, f"cm_raw{plot_suffix}.png"),
                    dpi=FIG_DPI)

    # Plot normalized confusion matrix
    plt.figure(constrained_layout=True)
    plot_confusion_matrix(cm, classes=labels, normalize=True,
                          title='Normalized confusion matrix')
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, f"cm_norm{plot_suffix}.png"),
                    dpi=FIG_DPI)

    if n_class == 2:
        for i in range(0, len(cm_list)):
            div_factor = cm_list[i].sum(axis=1)[:, np.newaxis]
            cm_list[i] = cm_list[i].astype('float') / div_factor
        cm_list = np.array(cm_list)
        cm_list_by_axis = [[], [], [], []]
        for i in range(0, len(cm_list)):
            cm_list_by_axis[0].append(cm_list[i][0][0])
            cm_list_by_axis[1].append(cm_list[i][0][1])
            cm_list_by_axis[2].append(cm_list[i][1][0])
            cm_list_by_axis[3].append(cm_list[i][1][1])
        cm_ci = np.ndarray([2, 2, 2])
        cm_ci[0, 0] = compute_confidence_intervals(cm_list_by_axis[0])
        cm_ci[0, 1] = compute_confidence_intervals(cm_list_by_axis[1])
        cm_ci[1, 0] = compute_confidence_intervals(cm_list_by_axis[2])
        cm_ci[1, 1] = compute_confidence_intervals(cm_list_by_axis[3])
        cm = sum(cm_list)
        plt.figure(constrained_layout=True)
        plot_confusion_matrix(cm, classes=labels, normalize=True, n_runs=runs,
                              ci_type="ci", ci=cm_ci,
                              title='Normalized confusion matrix')
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir,
                                     f"cm_norm_95CI{plot_suffix}.png"),
                        dpi=FIG_DPI)

    if output_dir is None:
        plt.show()

    # Output features by importance
    features_importance = pd.DataFrame(features_importance.items(),
                                       columns=['feature', 'importance'])
    features_importance['mean'] = features_importance.importance.apply(mean)
    features_importance['median'] = features_importance.importance.apply(
        median)
    features_importance['stdev'] = features_importance.importance.apply(stdev)
    if output_dir is None:
        print(tabulate(
            features_importance[['feature', 'mean', 'stdev']].sort_values(
                'mean', ascending=False), headers='keys', tablefmt='psql'))
    else:
        with open(os.path.join(output_dir, 'features_by_importance.csv'),
                  mode='w') as features_importance_output_file:
            features_importance_output_file.write(
                features_importance[['feature', 'mean', 'stdev']].sort_values(
                    'mean', ascending=False).to_csv(index=False))

    # Output top features
    s = [(k, top_features[k]) for k in
         sorted(top_features, key=top_features.get,
                reverse=True)]
    if output_dir is None:
        print(f"Most important {TOP_FEATURES_PERC}% of {len(features)} "
              "features (feature: count +1 per run when in the top "
              f"{TOP_FEATURES_PERC}%):")
        for k, v in s:
            print(f"{k}: {v / runs}")
    else:
        with open(os.path.join(output_dir, 'top_features.csv'),
                  mode='w') as top_features_output_file:
            top_features_output_file.write(
                f"Most important {TOP_FEATURES_PERC}% of {len(features)} "
                "features,(feature: count +1 per run when in the top "
                f"{TOP_FEATURES_PERC}%):")
            for k, v in s:
                top_features_output_file.write(f"\n{k},{v / runs}")

    # Output top features total score
    s = [(k, feature_scores[k]) for k in
         sorted(feature_scores, key=feature_scores.get,
                reverse=True)]
    if output_dir is None:
        print("Features total score (feature: count):")
        for k, v in s:
            print(f"{k}: {v / runs / len(features)}")
    else:
        with open(os.path.join(output_dir, 'feature_scores.csv'),
                  mode='w') as feature_scores_output_file:
            feature_scores_output_file.write("Features,total score")
            for k, v in s:
                feature_scores_output_file.write(
                    f"\n{k},{v / runs / len(features)}")



def run_validation(data_file, output_dir, runs, subsampling,
                   plot_suffix, test_data_file):
    """
    data_file must be a csv with control class (e.g. healthy plasma) before
    case class (e.g. ovarian cancer plasma) in order to have correct tn, fp,
    fn, tp order (see sklearn confusion matrix documentation).
    """
    # np.random.seed(0)

    data = pd.read_csv(data_file)
    test_data = pd.read_csv(test_data_file)

    # Create a list of the feature column's names
    features = data.columns[2:]
    n_top_features = round(len(features) * TOP_FEATURES_PERC / 100)

    if test_data_file is not None:
        if not np.array_equal(test_data.columns[2:], features):
            print("Features aren't identical between train and test data.")
            return

    # Create a random forest classifier
    clf = RandomForestClassifier(n_estimators=300, criterion='gini')

    unique_biological_classes = data.biological_class.unique()
    n_class = len(unique_biological_classes)

    predictions_dict = dict()
    misclassified_samples = dict()
    features_importance = dict()
    top_features = dict()
    feature_scores = dict() 
    fprs = dict()
    tprs = dict()
    for i in range(n_class):
        fprs[i] = list()
        tprs[i] = list()

    cross_val_tprs = dict()
    cross_val_aucs = dict()
    for i in range(n_class):
        cross_val_tprs[i] = list()
        cross_val_aucs[i] = list()
    cross_val_mean_fpr = np.linspace(0, 1, 100)
    if len(test_data.biological_class.unique()) > 1:
        cross_val_fig, cross_val_ax = plt.subplots(
            figsize=(7, 6.5), constrained_layout=True)
    cm = None
    cm_list = list()

    for run_n in tqdm(range(0, runs)):
        if subsampling is True:
            biological_classes_count = list()
            for biological_class in unique_biological_classes:
                biological_classes_count.append(
                    len(data[data['biological_class'] == biological_class]))
            n_samples_min = min(biological_classes_count)
            train_data = data[data['biological_class']
                              == unique_biological_classes[0]].sample(
                                n_samples_min)
            for i in range(1, len(unique_biological_classes)):
                train_data = train_data.append(
                    data[data['biological_class']
                         == unique_biological_classes[i]].sample(
                             n_samples_min))
        else:
            train_data = data
        y_train, labels = pd.factorize(train_data['biological_class'])

        clf.fit(train_data[features], y_train)

        y_test_raw, labels_test = pd.factorize(test_data['biological_class'])
        y_test = y_test_raw
        for i in range(len(y_test_raw)):
            y_test[i] = labels.get_loc(labels_test[y_test_raw[i]])

        predictions = clf.predict_proba(test_data[features])
        for i, pred in enumerate(predictions):
            pred_sample = test_data.iloc[i].loc['sample']
            if pred_sample not in predictions_dict.keys():
                predictions_dict[pred_sample] = dict()
                predictions_dict[pred_sample]['biological_class'] = \
                    test_data.iloc[i].loc['biological_class']
                for j in range(n_class):
                    predictions_dict[pred_sample][f'pred_{labels[j]}'] = 0
                predictions_dict[pred_sample]['n'] = 0
            for j in range(n_class):
                predictions_dict[pred_sample][f'pred_{labels[j]}'] += pred[j]
            predictions_dict[pred_sample]['n'] += 1

        y_pred = clf.predict(test_data[features])

        y_pred_labeled = labels[y_pred]

        feature_importances_df = pd.DataFrame({'features': features, 'importance': clf.feature_importances_})

        for feature, importance in feature_importances_df.values:
            if feature not in features_importance.keys():
                features_importance[feature] = list()
            features_importance[feature].append(importance)

        # Aggregate top features
        feature_importances_df_sorted = feature_importances_df.sort_values('importance', ascending=False)
        for feature in feature_importances_df_sorted.iloc[0:n_top_features]['features']:
            if feature not in top_features.keys():
                top_features[feature] = 0
            top_features[feature] += 1

        feature_importances_df_sorted = feature_importances_df.sort_values('importance', ascending=True)
        for i in range(0, len(feature_importances_df_sorted) - 1):
            feature = feature_importances_df_sorted.iloc[i]['features']
            if feature not in feature_scores.keys():
                feature_scores[feature] = 0
            feature_scores[feature] += i + 1

        # Build misclassified items dictionary
        misclassified = test_data['biological_class'] == y_pred_labeled
        for i in misclassified.index:
            if test_data['sample'][i] not in misclassified_samples.keys():
                misclassified_samples[test_data['sample'][i]] = 0
            if misclassified[i] == False:  # noqa
                misclassified_samples[test_data['sample'][i]] += 1

        if len(test_data.biological_class.unique()) > 2:
            for i in range(n_class):
                fpr_RF, tpr_RF, _ = roc_curve(
                    y_test, predictions[:, i], pos_label=i)
                fprs[i].append(fpr_RF)
                tprs[i].append(tpr_RF)

                logging.disable(logging.WARNING)
                cross_val_viz = RocCurveDisplay.from_predictions(
                    y_test, predictions[:, i], pos_label=i,
                    label=None,
                    alpha=0.0,
                    lw=1,
                    ax=cross_val_ax)
                logging.disable(logging.NOTSET)
                cross_val_interp_tpr = np.interp(
                    cross_val_mean_fpr, cross_val_viz.fpr, cross_val_viz.tpr)
                # cross_val_interp_tpr[0] = 0.0
                cross_val_tprs[i].append(cross_val_interp_tpr)
                cross_val_aucs[i].append(cross_val_viz.roc_auc)
        else:
            fpr_RF, tpr_RF, _ = roc_curve(
                y_test, predictions[:, 1], pos_label=1)
            fprs[0].append(fpr_RF)
            tprs[0].append(tpr_RF)

            logging.disable(logging.WARNING)
            cross_val_viz = RocCurveDisplay.from_predictions(
                y_test, predictions[:, 1], pos_label=1,
                label=None,
                alpha=0.0,
                lw=1,
                ax=cross_val_ax)
            logging.disable(logging.NOTSET)
            cross_val_interp_tpr = np.interp(
                cross_val_mean_fpr, cross_val_viz.fpr, cross_val_viz.tpr)
            # cross_val_interp_tpr[0] = 0.0
            cross_val_tprs[0].append(cross_val_interp_tpr)
            cross_val_aucs[0].append(cross_val_viz.roc_auc)

        if cm is None:
            cm = confusion_matrix(test_data['biological_class'],
                                  y_pred_labeled, labels=labels)
        else:
            cm += confusion_matrix(test_data['biological_class'],
                                   y_pred_labeled, labels=labels)
        cm_list.append(confusion_matrix(test_data['biological_class'],
                                        y_pred_labeled, labels=labels))

    print("Generating outputs...")

    if output_dir is not None:
        if not os.path.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except OSError:
                print(f"Creation of directory '{output_dir}' failed.")

    # Output prediction scores
    output_prediction_scores(predictions_dict, labels, output_dir)

    # Output misclassified samples
    samples_info = load_samples_info(SAMPLES_INFO_FILE)
    s = [(k, misclassified_samples[k]) for k in
         sorted(misclassified_samples, key=misclassified_samples.get,
                reverse=True)]
    misclassified_samples_output = list()
    misclassified_samples_output.append(
        "sample,biological_class,misclassified_n")
    for k, v in s:
        try:
            biological_class = samples_info.loc[k]['biological_class']
        except KeyError:
            biological_class = "unknown_biological_class"
        misclassified_samples_output.append(f"{k},{biological_class},{v}")
    if output_dir is None:
        print('\n'.join(misclassified_samples_output))
    else:
        with open(os.path.join(output_dir, 'misclassified_samples.csv'),
                  mode='w') as misclassified_samples_output_file:
            misclassified_samples_output_file.write(
                '\n'.join(misclassified_samples_output))

    if len(test_data.biological_class.unique()) > 1:
        roc_labels = [labels[1]] if len(test_data.biological_class.unique()) == 2 else labels
        # Plot ROC curve crossval
        cross_val_fig.set_constrained_layout_pads(w_pad=.15, h_pad=.26)
        plot_roc_curve_multiclass_crossval(cross_val_ax, cross_val_tprs,
                                           cross_val_mean_fpr, cross_val_aucs,
                                           roc_labels, ci_type="ci")
        if output_dir is not None:
            plt.savefig(
                os.path.join(output_dir,
                             f"roc_multiclass_crossval_ci{plot_suffix}.png"),
                dpi=FIG_DPI)

        # Output ROC multiclass
        plot_roc_curve_multiclass(fprs, tprs, roc_labels)
        if output_dir is not None:
            plt.savefig(os.path.join(output_dir,
                                     f"roc_multiclass{plot_suffix}.png"),
                        dpi=FIG_DPI)

        sensitivity_specificity_output = sensitivity_specificity_table(
            cross_val_tprs, roc_labels)
        if output_dir is None:
            print('\n'.join(sensitivity_specificity_output))
        else:
            with open(os.path.join(output_dir, 'sensitivity_specificity.csv'),
                    mode='w') as sensitivity_specificity_output_file:
                sensitivity_specificity_output_file.write(
                    '\n'.join(sensitivity_specificity_output))

    np.set_printoptions(precision=2)

    # Plot non-normalized confusion matrix
    plt.figure(constrained_layout=True)
    plot_confusion_matrix(cm, classes=labels,
                          title='Confusion matrix, without normalization')
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, f"cm_raw{plot_suffix}.png"),
                    dpi=FIG_DPI)

    # Plot normalized confusion matrix
    plt.figure(constrained_layout=True)
    plot_confusion_matrix(cm, classes=labels, normalize=True,
                          title='Normalized confusion matrix')
    if output_dir is not None:
        plt.savefig(os.path.join(output_dir, f"cm_norm{plot_suffix}.png"),
                    dpi=FIG_DPI)

    if output_dir is None:
        plt.show()

    # Output features by importance
    features_importance = pd.DataFrame(features_importance.items(),
                                       columns=['feature', 'importance'])
    features_importance['mean'] = features_importance.importance.apply(mean)
    features_importance['median'] = features_importance.importance.apply(
        median)
    features_importance['stdev'] = features_importance.importance.apply(stdev)
    if output_dir is None:
        print(tabulate(
            features_importance[['feature', 'mean', 'stdev']].sort_values(
                'mean', ascending=False), headers='keys', tablefmt='psql'))
    else:
        with open(os.path.join(output_dir, 'features_by_importance.csv'),
                  mode='w') as features_importance_output_file:
            features_importance_output_file.write(
                features_importance[['feature', 'mean', 'stdev']].sort_values(
                    'mean', ascending=False).to_csv(index=False))

    # Output top features
    s = [(k, top_features[k]) for k in
         sorted(top_features, key=top_features.get,
                reverse=True)]
    if output_dir is None:
        print(f"Most important {TOP_FEATURES_PERC}% of {len(features)} "
              "features (feature: count +1 per run when in the top "
              f"{TOP_FEATURES_PERC}%):")
        for k, v in s:
            print(f"{k}: {v / runs}")
    else:
        with open(os.path.join(output_dir, 'top_features.csv'),
                  mode='w') as top_features_output_file:
            top_features_output_file.write(
                f"Most important {TOP_FEATURES_PERC}% of {len(features)} "
                "features,(feature: count +1 per run when in the top "
                f"{TOP_FEATURES_PERC}%):")
            for k, v in s:
                top_features_output_file.write(f"\n{k},{v / runs}")

    # Output top features total score
    s = [(k, feature_scores[k]) for k in
         sorted(feature_scores, key=feature_scores.get,
                reverse=True)]
    if output_dir is None:
        print("Features total score (feature: count):")
        for k, v in s:
            print(f"{k}: {v / runs / len(features)}")
    else:
        with open(os.path.join(output_dir, 'feature_scores.csv'),
                  mode='w') as feature_scores_output_file:
            feature_scores_output_file.write("Features,total score")
            for k, v in s:
                feature_scores_output_file.write(
                    f"\n{k},{v / runs / len(features)}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file", help="input cg_methyl/haplotypes data file")
    parser.add_argument("--test",
                        type=str,
                        default=None,
                        help="test cg_methyl/haplotypes data file")
    parser.add_argument("-r", "--runs",
                        type=int,
                        default=3000,
                        help="number of runs")
    parser.add_argument("-s", "--subsampling",
                        default=False,
                        action='store_true',
                        help="enable subsampling")
    parser.add_argument("-o", "--output_dir",
                        type=str,
                        default=None,
                        help="output directory")
    parser.add_argument("--plot_suffix",
                        type=str,
                        default="",
                        help="plot suffix")
    args = parser.parse_args()
    if args.test is not None:
        run_validation(args.data_file,
                       output_dir=args.output_dir,
                       runs=args.runs, subsampling=args.subsampling,
                       plot_suffix=args.plot_suffix,
                       test_data_file=args.test)
    else:
        run_classification(args.data_file,
                           output_dir=args.output_dir,
                           runs=args.runs, subsampling=args.subsampling,
                           plot_suffix=args.plot_suffix)


if __name__ == '__main__':
    main()
