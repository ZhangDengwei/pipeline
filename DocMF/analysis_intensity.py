#!/usr/bin/python
# -*- coding: ascii -*-

import re
import gzip
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def openfile(file):
    if re.search('.gz$', file):
        openfile = gzip.open(file, 'rb')
        with openfile as f:
            while True:
                line = f.readline().decode().rstrip('\n')
                if not line:
                    break
                yield line
    else:
        openfile = open(file, 'r')
        with openfile as f:
            while True:
                line = f.readline().rstrip('\n')
                if not line:
                    break
                yield line

def convertseq(target, sequence):
    regex = {'N': '[ATCG]', 'Y': ['TC'], 'R': ['AG'], 'W': ['AT'],
             'S': ['CG'], 'M': ['AC'], 'K': ['GT'],
             'H': ['ATC'], 'B': ['GTC'], 'V': ['GAC'], 'D': ['GAT'], 
             'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G'}
    pattern = ''
    for i in target:
        pattern += regex[i]
    if re.search(pattern, sequence):
        return True
    else:
        return False

def plot(target_sets, non_target_sets, title):
    fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,6))
    plt.title(str('Histogram of ' + title))
    sns_plot = sns.distplot(target_sets, kde=False)
    fig = sns_plot.get_figure()
    fig.savefig("target_hist.png")

    fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,6))
    plt.title(str('Histogram of non-' + title))
    sns_plot = sns.distplot(non_target_sets, kde=False)
    fig = sns_plot.get_figure()
    fig.savefig("non-target_hist.png")

    fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,6))
    plt.title(str('Kdeplot of ' + title))
    sns_plot = sns.kdeplot(target_sets, shade=True)
    fig = sns_plot.get_figure()
    fig.savefig("target_kdeplot.png")

    fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,6))
    plt.title(str('Kdeplot of non-' + title))
    sns_plot = sns.kdeplot(non_target_sets, shade=True)
    fig = sns_plot.get_figure()
    fig.savefig("non-target_kdeplot.png")

    fig,axes=plt.subplots(nrows=1,ncols=1,figsize=(8,6))
    sns_plot = sns.distplot(target_sets, hist=False, label=title, color='r')
    sns_plot = sns.distplot(non_target_sets, hist=False, label=str('non-'+title), color='b')
    fig = sns_plot.get_figure()
    fig.savefig("two_group_kdeplot.png")



def analysis_intensity(infile, splitseq, start, end):
    read_infile = openfile(infile) 
    target_seq_sets, non_target_seq_sets = [], []
    target_seq_counts, non_target_seq_sets_counts = 0, 0
    target_seq_total_intensity, non_target_seq_total_intensity = 0, 0

    with open("Different_type_intensity.log", "w") as fout:
        for r in iter(read_infile):
            sequence = r.split('\t')[1]
            intensity = r.split('\t')[2]
            if 'N' not in sequence:
                if convertseq(splitseq, sequence[start:end]):
                    target_seq_sets.append(float(intensity))
                    target_seq_counts += 1
                    target_seq_total_intensity += float(intensity)
                else:
                    non_target_seq_sets.append(float(intensity))
                    non_target_seq_sets_counts += 1
                    non_target_seq_total_intensity += float(intensity)
            else:
                continue

        print("target sequences counts: ", target_seq_counts, file=fout)
        print("total target sequences intensity: ", target_seq_total_intensity, file=fout)
        print("non_target sequences counts: ", non_target_seq_sets_counts, file=fout)
        print("total non_target sequences intensity: ", non_target_seq_total_intensity, file=fout)        

    plot(target_seq_sets, non_target_seq_sets, splitseq)

def main():
    parse = argparse.ArgumentParser(description="split intensity for partcular sequence")
    parse.add_argument("--infile", "-i", help="input file", required=True)
    parse.add_argument("--target_seq", "-t", help="target sequence, eg. NGG", required=True)
    parse.add_argument("--start", "-s", help="the beginning position, the 1st coordinate is 0", required=True, type=int)
    parse.add_argument("--end", "-e", help="the ending position", required=True, type=int)

    args = parse.parse_args()

    analysis_intensity(args.infile, args.target_seq, args.start, args.end)


if __name__ == "__main__":
    main()

