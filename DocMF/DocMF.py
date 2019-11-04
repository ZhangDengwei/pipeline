import re
import gzip
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
        with fastq as f:
            while True:
                l1 = f.readline().decode().rstrip('\n')
                if not l1:
                    break
                l2 = f.readline().decode().rstrip('\n')
                l3 = f.readline().decode().rstrip('\n')
                l4 = f.readline().decode().rstrip('\n')
                yield [l1, l2, l3, l4]
    else:
        fastq = open(file, 'r')
        with fastq as f:
            while True:
                l1 = f.readline().rstrip('\n')
                if not l1:
                    break
                l2 = f.readline().rstrip('\n')
                l3 = f.readline().rstrip('\n')
                l4 = f.readline().rstrip('\n')
                yield [l1, l2, l3, l4]

def plot(x):
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    sns_plot_1 = sns.distplot(x, kde=False)
    fig_1 = sns_plot_1.get_figure()
    fig_1.savefig("Histogram_of_Fluorescence_intensity_for_positive_data.png")

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    sns_plot_2 = sns.boxplot(x)
    fig_2 = sns_plot_2.get_figure()
    fig_2.savefig("Box_plot_of_Fluorescence_intensity_for_positive_data.png")

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    sns_plot_3 = sns.distplot(x, kde=True)
    fig_3 = sns_plot_3.get_figure()
    fig_3.savefig("Intensity_plot_of_Fluorescence_intensity_for_positive_data.png")

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    sns_plot_4 = sns.violinplot(x)
    fig_4 = sns_plot_4.get_figure()
    fig_4.savefig("Violin_plot_of_Fluorescence_intensity_for_positive_data.png")
 

def split(fastq_filename, target_sequence, start, end, threshold_foldchange, length, DocMF_type):
    fastq_file = fq(fastq_filename)
    target_complementary_reverse = target_sequence.translate(str.maketrans("ATCG","TAGC"))[::-1]
    read_number = 0
    pass_filter_reads_number = 0
    left_seq_set, right_seq_set = {}, {}
    data = []

    with gzip.open("left_intermediate.gz", 'wb') as out_f1, gzip.open("right_intermediate.gz", 'wb') as out_f2, open("log.txt", "w") as log:
        for r1 in iter(fastq_file):
            read_number += 1
            read_id = r1[0]
            read_sequence = r1[1]
            intensity_Before = float(r1[2].split()[1]) 
            intensity_After = float(r1[2].split()[2]) if float(r1[2].split()[2]) > 0 else 0
            if float(intensity_Before) > 0: # intensity_Before ought to be over 0, otherwise this read was FALSE.

                fold_change_of_intensity = intensity_After / intensity_Before # need to be verified
                if DocMF_type == "cutting":
                    # A smaller fold_change_of_intensity represents a higher cutting efficiency.

                    if fold_change_of_intensity <= float(threshold_foldchange):
                        data.append(fold_change_of_intensity)
                        if read_sequence[start:end] == target_sequence:
                            pass_filter_reads_number += 1
                            left_sequence = read_sequence[:start][-length:]
                            right_sequence = read_sequence[end:][:length]
                        elif read_sequence[start:end] == target_complementary_reverse:
                            pass_filter_reads_number += 1
                            new_sequence = read_sequence.translate(str.maketrans("ATCG","TAGC"))[::-1]
                            new_start = len(new_sequence) - end
                            new_end = len(target_sequence) + new_start
                            left_sequence = new_sequence[:new_start][-length:]
                            right_sequence = new_sequence[new_end:][:length]
                        else:
                            continue

                        # To counte the total number of each type of sequence
                        if "N" not in left_sequence:
                            if left_seq_set.get(left_sequence, 0) == 0:
                                left_seq_set[left_sequence] = [1, fold_change_of_intensity]
                            else:   
                                left_seq_set[left_sequence][0] += 1
                                left_seq_set[left_sequence][1] += fold_change_of_intensity
                        else:
                            pass
                        if "N" not in right_sequence:
                            if right_seq_set.get(right_sequence, 0) == 0:
                                right_seq_set[right_sequence] = [1, fold_change_of_intensity]
                            else:
                                right_seq_set[right_sequence][0] += 1
                                right_seq_set[right_sequence][1] += fold_change_of_intensity
                        else:
                            pass
                        write_left = read_id + "\t" + left_sequence + "\t" + str(fold_change_of_intensity) + "\n"
                        write_right = read_id + "\t" + right_sequence + "\t" + str(fold_change_of_intensity) + "\n"
                        out_f1.write(write_left.encode())
                        out_f2.write(write_right.encode())
                elif DocMF_type == "binding":
                    # One with bigger fold_change_of_intensity represents a higher binding efficiency.

                    if fold_change_of_intensity >= float(threshold_foldchange):
                        data.append(fold_change_of_intensity)
                        if read_sequence[start:end] == target_sequence:
                            pass_filter_reads_number += 1
                            left_sequence = read_sequence[:start][-length:]
                            right_sequence = read_sequence[end:][:length]
                        elif read_sequence[start:end] == target_complementary_reverse:
                            pass_filter_reads_number += 1
                            new_sequence = read_sequence.translate(str.maketrans("ATCG","TAGC"))[::-1]
                            new_start = len(new_sequence) - end
                            new_end = len(target_sequence) + new_start
                            left_sequence = new_sequence[:new_start][-length:]
                            right_sequence = new_sequence[new_end:][:length]
                        else:
                            continue

                        # To counte the total number of each type of sequence
                        if "N" not in left_sequence:
                            if left_seq_set.get(left_sequence, 0) == 0:
                                left_seq_set[left_sequence] = [1, fold_change_of_intensity]
                            else:   
                                left_seq_set[left_sequence][0] += 1
                                left_seq_set[left_sequence][1] += fold_change_of_intensity
                        else:
                            pass
                        if "N" not in right_sequence:
                            if right_seq_set.get(right_sequence, 0) == 0:
                                right_seq_set[right_sequence] = [1, fold_change_of_intensity]
                            else:
                                right_seq_set[right_sequence][0] += 1
                                right_seq_set[right_sequence][1] += fold_change_of_intensity
                        else:
                            pass
                        write_left = read_id + "\t" + left_sequence + "\t" + str(fold_change_of_intensity) + "\n"
                        write_right = read_id + "\t" + right_sequence + "\t" + str(fold_change_of_intensity) + "\n"
                        out_f1.write(write_left.encode())
                        out_f2.write(write_right.encode())
            else:
                continue
        print("Total reads before filter: ", str(read_number), file=log)
        print("total reads after filter: ", str(pass_filter_reads_number), file=log)
        
    plot(data)

    # output final report
    with open("left_report.txt", "w") as f1, open("right_report.txt", "w") as f2:
        print("Sequence\tCounts\ttotal_fluorescence_intensity", file=f1)
        print("Sequence\tCounts\ttotal_fluorescence_intensity", file=f2)
        for i,j in sorted(left_seq_set.items(), key=lambda item:item[1][1], reverse=True):
            print(i, "\t", str(j[0]), "\t", str(j[1]), file=f1)
        for i,j in sorted(right_seq_set.items(), key=lambda item:item[1][1], reverse=True):
            print(i, "\t", str(j[0]), "\t", str(j[1]), file=f2)


def main():
    parse = argparse.ArgumentParser(description="split sequencing reads on DocMF platform")
    parse.add_argument("--fastq", "-f", help="fastq filename with special format", required=True)
    parse.add_argument("--target_sequence", "-t", help="target sequence", required=True)
    parse.add_argument("--start", "-s", help="the beginning position of target sequence in entire sequence", required=True, type=int)
    parse.add_argument("--end", "-e", help="the ending position of target sequence in entire sequence", required=True, type=int)
    parse.add_argument("--threshold", "-c", help="the threshold of fluorescence intensity, default: 0.3", default=0.3, required=True)
    parse.add_argument("--length", "-L", help="the fetched length of randomized N sequence", required=True, type=int)
    parse.add_argument("--DocMF_type", "-D", help="the experimental type", choices = ['cutting', 'binding'], required=True)

    args = parse.parse_args()

    split(args.fastq, args.target_sequence, args.start, args.end, args.threshold, args.length, args.DocMF_type)


if __name__ == "__main__":
    main()

