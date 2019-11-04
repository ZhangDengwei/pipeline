import argparse
import os
import sys
import subprocess
import logging
import re
import shlex

'''
dsCir_seq.py pipeline is stemmed from circleseq.py based on a published paper - 
CIRCLE-seq: A highly sensitive in vitro screen for genome-wide CRISPR-Cas9 
nuclease off-targets, whose github address is "https://github.com/tsailabSJ/circleseq".

However, original circleseq pipeline would merge pair-end reads into one, then 
being aligned to reference. This dsCir-Seq extracts aligned result from the 
mapping result of pair-end reads, and judges the direction of library DNA whether
 it is circular or linear according to its mapping location.

Some steps in circleseq pipeline have been simplified and modified.
'''

'''
Procedure of dsCir_seq:
01. reads QC using fastp
02. align pair-end reads to reference using BWA
03. extract mapping information and generate report
04. summarize sequencing depth and uniformity
'''

'''
User must offer the two following file with requested format:
#1: sample.txt (seven columns separated by tab-delimited)
column 1: sample name
column 2: absolute path of read 1
column 3: absolute path of read 2
column 4: target site
column 5: chromosome
column 6: start
column 7: end

#2: adapter.fa (contains all adapter sequences with Fasta format)
>adapter_1
ATGCG...
>adapter_2
ACGTT...
>adapter_3
AGGTT...
'''


def summarize_offtarget(file):
    # extract the information of off-target in log file generated in step3
    with open(file, 'r') as f:
        all_line = str(f.readlines())
        contaminated_reads = re.search("(The counts of liner contaminate reads).+?(\d+)", all_line).group(2)
        dsCir_reads = re.search("(The counts of ds-circle reads).+?(\d+)", all_line).group(2)
        overlap_reads = re.search("(The counts of mapped reads with right).+?:.+?(\d+)", all_line).group(2)
        broken_points = re.search("(Total broken points).+?(\d+)", all_line).group(2)
        offtarget_sites = re.search("(Total off-target sites).+?(\d+)", all_line).group(2)
    return contaminated_reads, dsCir_reads, overlap_reads, broken_points, offtarget_sites


def extract_mapping_information(file):
    # extract mapping information from the file - mapping.statistic
    with open(file, 'r') as f:
        all_line = str(f.readlines())
        total_reads = re.search("(\d+)?\s\\+.+?in total", all_line).group(1)
        mapped_reads = re.search("duplicates.+?(\d+)?\s\\+.+?mapped\s\\((\d+\\.\d+%).+", all_line).group(1)
        mapping_rate = re.search("duplicates.+?(\d+)?\s\\+.+?mapped\s\\((\d+\\.\d+%).+", all_line).group(2)
    return total_reads, mapped_reads, mapping_rate


def compute_uniformity(file_1, file_2):
    # compute uniformity from the file in step4
    with open(file_1, 'r') as f1, open(file_2, 'r') as f2:
        all_line_1 = "".join(f1.readlines())
        target_reads = re.search("\\[Target\\]\sTarget\sReads\t(\d+)\n", all_line_1).group(1)
        average_depth = re.search("\\[Target\\]\s+Average\s+depth\s+(\d+.+?)\n", all_line_1).group(1)
        total_base, over_uniformity_base = 0, 0
        for line in f2:
            items = line.strip('\n').split()
            total_base += int(items[1])
            if int(items[0]) >= (float(average_depth) * 0.2):
                over_uniformity_base += int(items[1])
            else:
                pass
        uniformity = '%.2f%%' % (over_uniformity_base / total_base * 100)
        average_depth_format = average_depth + 'x'
    return target_reads, average_depth_format, uniformity


def import_sample(sample_filename):
    dic_sample = {}
    with open(sample_filename, "r") as f:
        for line in f:
            items = line.rstrip("\n").split("\t")
            sample_id = items[0]
            read_1 = items[1]
            read_2 = items[2]
            target = items[3]
            chromosome = items[4]
            start = items[5]
            end = items[6]
            dic_sample[sample_id] = (read_1, read_2, target, chromosome, start, end)
        return dic_sample


def reads_qc(dic_sample, adapter_filename):
    current_path = os.getcwd()
    for m, (j, k, l, c, s, e) in dic_sample.items():
        os.chdir(current_path)
        new_path = os.path.join(current_path, m)
        if not os.path.exists(new_path):
            os.mkdir(new_path)
        else:
            pass
        os.chdir(new_path)
        read_1_basename = os.path.basename(j)
        read_2_basename = os.path.basename(k)
        clean_read_1_name = 'clean_' + read_1_basename
        clean_read_2_name = 'clean_' + read_2_basename
        logging.info("Start QC for %s ..." % m)
        # check whether this step has been done
        if os.path.exists(clean_read_1_name) and os.path.exists(clean_read_2_name):
            logging.warning('The files of clean reads have been existed, skipping QC step!!! If not, please delete them'
                            ' and restart...')
        else:
            command_qc = 'fastp -i ' + j + ' -o ' + clean_read_1_name + ' -I ' + k + ' -O ' + clean_read_2_name + \
                         ' --adapter_fasta ' + adapter_filename + \
                         ' --thread 6 --length_required 30 -j report.json -h report.html -R ' + m + ' > log.txt 2>&1'
            logging.info("Command: %s" % command_qc)
            subprocess.check_call(command_qc, shell=True)
    logging.info("Finish reads QC.")


def align_bwa(dic_sample, bwa_index, home_path):
    current_path = os.getcwd()
    for m, (j, k, l, c, s, e) in dic_sample.items():
        os.chdir(current_path)
        new_path = os.path.join(current_path, m)
        if not os.path.exists(new_path):
            os.mkdir(new_path)
        else:
            pass
        os.chdir(new_path)
        read_1_basename = os.path.basename(j)
        read_2_basename = os.path.basename(k)
        clean_read_1_name = 'clean_' + read_1_basename
        clean_read_2_name = 'clean_' + read_2_basename
        clean_read_1_path = os.path.join(home_path, '01.reads_QC', m, clean_read_1_name)
        clean_read_2_path = os.path.join(home_path, '01.reads_QC', m, clean_read_2_name)

        logging.info("Start aligning %s to reference ..." % m)
        step_complete = True
        generate_file = ['bwa.align.log', 'align.pe.bwa.sorted.bam', 'mapping.statistic']
        # check whether this step has been done
        for i in generate_file:
            if not os.path.exists(os.path.join(new_path, i)):
                step_complete = False
                break
        if not step_complete:
            command_aln = 'bwa mem -t 6 ' + bwa_index + ' ' + clean_read_1_path + ' ' + clean_read_2_path + \
                          ' -o align.pe.bwa.sam > bwa.align.log 2>&1'
            logging.info("Command: %s" % command_aln)
            subprocess.check_call(command_aln, shell=True)
            samtool_command_1 = "samtools view -bS align.pe.bwa.sam | samtools sort > align.pe.bwa.sorted.bam"
            samtool_command_2 = "samtools flagstat align.pe.bwa.sorted.bam > mapping.statistic"
            subprocess.check_call(samtool_command_1, shell=True)
            subprocess.check_call(samtool_command_2, shell=True)
        else:
            logging.warning("The SAM and BAM files have been existed, skipping alignment!!! If not, please delete them"
                            " and restart...")
    logging.info("Finish reads alignment.")


def extract_cleavage(dic_sample, home_path, script_path, python_command, genome, gtf, search_radius=20, gap=5,
                     overlap_threshold=5, mismatch_threshold=6, top_number=50):
    current_path = os.getcwd()
    for m, (j, k, l, c, s, e) in dic_sample.items():
        os.chdir(current_path)
        new_path = os.path.join(current_path, m)
        if not os.path.exists(new_path):
            os.mkdir(new_path)
        else:
            pass
        os.chdir(new_path)
        sam_path = os.path.join(home_path, '02.mapping', m, 'align.pe.bwa.sam')
        analysis_script_path = os.path.join(script_path, 'analysis_paired.py')
        logging.info("Start extracting %s's cleavage information ..." % m)
        step_complete = True
        generate_file = ['coordinate_output.txt', 'offtarget_output.txt', 'unmatch_output.txt', 'log.txt']
        # check whether this step has been done
        for i in generate_file:
            if not os.path.exists(os.path.join(new_path, i)):
                step_complete = False
                break
        if not step_complete:
            command_extraction = python_command + ' ' + analysis_script_path + ' -r ' + genome + ' -g ' + gtf + ' -S ' + \
                                 sam_path + ' -t ' + l + ' --gap ' + str(gap) + ' --search_radius ' + str(search_radius) + \
                                 ' --overlap_threshold ' + str(overlap_threshold) + ' --mismatch_threshold ' + \
                                 str(mismatch_threshold) + ' -c coordinate_output.txt -o offtarget_output.txt ' \
                                                           '-u unmatch_output.txt > log.txt 2>&1'
            logging.info("Command: %s" % command_extraction)
            subprocess.check_call(command_extraction, shell=True)
            logging.info("Finish %s's cleavage analysis" % m)
        else:
            logging.warning("The extraction has been done in advance, skipping extracting off-target sites!!! If not, "
                            "please delete them and restart...")

        # visualization
        logging.info("Visualize %s's offtarget sites ..." % m)
        visualization_script_path = os.path.join(script_path, 'visualization.py')
        command_visualization = python_command + ' ' + visualization_script_path + ' -i offtarget_output.txt -s ' + l + \
                                ' --top_number ' + str(top_number) + ' -o visual -t ' + m
        logging.info("Command: %s" % command_visualization)
        subprocess.check_call(command_visualization, shell=True)
    logging.info("Finish cleavage extraction and off-target visualization.")


def bamdst(dic_sample, home_path):
    current_path = os.getcwd()
    for m, (j, k, l, c, s, e) in dic_sample.items():
        os.chdir(current_path)
        new_path = os.path.join(current_path, m)
        if not os.path.exists(new_path):
            os.mkdir(new_path)
        else:
            pass
        os.chdir(new_path)
        content = shlex.quote(c + '\t' + s + '\t' + e)
        command_bed = 'echo -e ' + content + ' > target.bed'
        subprocess.check_call(command_bed, shell=True)
        logging.info("Summarize sequencing information of %s ..." % m)

        step_complete = True
        generate_file = ['chromosomes.report', 'coverage.report', 'depth_distribution.plot', 'depth.tsv.gz',
                         'insertsize.plot', 'region.tsv.gz', 'uncover.bed']
        # check whether this step has been done
        for i in generate_file:
            if not os.path.exists(os.path.join(new_path, i)):
                step_complete = False
                break
        if not step_complete:
            bed_file_path = os.path.join(new_path, 'target.bed')
            bam_file_path = os.path.join(home_path, '02.mapping', m, 'align.pe.bwa.sorted.bam')
            command_bamdst = 'bamdst --cutoffdepth 10 --maxdepth 10000 -q 10 -p ' + bed_file_path + ' -o ./ ' + \
                             bam_file_path + ' > log.txt 2>&1'
            logging.info("Summarize %s's sequencing information ..." % m)
            logging.info("Command: %s" % command_bamdst)
            subprocess.check_call(command_bamdst, shell=True)
        else:
            logging.warning("The 'bamdst' step has been done, skipping this step!!! If not, please delete these files "
                            "and restart...")


def report(dic_sample, home_path):
    dic_info = dict()
    for m, (j, k, l, c, s, e) in dic_sample.items():
        mapping_file_path = os.path.join(home_path, '02.mapping', m, 'mapping.statistic')
        total_reads, mapped_reads, mapping_rate = extract_mapping_information(mapping_file_path)
        offtarget_log_path = os.path.join(home_path, '03.extract_cleavage', m, 'log.txt')
        linear_reads, circular_reads, overlap_reads, broken_points, offtarget_sites = summarize_offtarget(
            offtarget_log_path)
        uniformity_file1_path = os.path.join(home_path, '04.bamdst', m, 'coverage.report')
        uniformity_file2_path = os.path.join(home_path, '04.bamdst', m, 'depth_distribution.plot')
        target_reads, average_depth, uniformity = compute_uniformity(uniformity_file1_path, uniformity_file2_path)
        dic_info[m] = (total_reads, mapped_reads, mapping_rate, linear_reads, circular_reads, overlap_reads,
                       broken_points, offtarget_sites, target_reads, average_depth, uniformity)

        # remove SAM file
        sam_file_path = os.path.join(home_path, '02.mapping', m, 'align.pe.bwa.sam')
        if os.path.exists(sam_file_path):
            os.remove(sam_file_path)
        else:
            pass

    with open("Final_statistic_report.txt", 'a') as fo:
        print("#" * 100, file=fo)
        print("Note:", file=fo)
        print("01. Linear_reads represent the contaminated reads derived from linear DNA fragment, its pattern is as"
              " follows:", file=fo)
        print("\t\t---------->\t<----------", file=fo)
        print("02. Circular_reads represent reads which are stemmed from circular DNA fragment, its pattern is as"
              " follows:", file=fo)
        print("\t\t<----------\t---------->", file=fo)
        print("03. overlap_reads represent reads derived from circular DNA but there is a overlap between two reads"
              " after locating on the reference, its pattern is as follows:", file=fo)
        print("\t\t<----------", file=fo)
        print("\t\t       ---------->", file=fo)
        print("04: Uniformity is the ratio that the sum of bases, whose depth is bigger than the 20 percent of average,"
              " divides by the total bases of target sites", file=fo)
        print("#" * 100, file=fo)
        print("\n\n", file=fo)

        print("Sample\tTotal_reads\tMapped_reads\tMapping_rate\tLinear_reads\tCircular_reads\tOverlap_reads\t"
              "Broken_points\tOff-target_sites\tOn-target_reads\tAverage_depth\tUniformity", file=fo)
        for k, j in dic_info.items():
            print(k, "\t", str(j[0]), "\t", str(j[1]), "\t", j[2], "\t", str(j[3]), "\t", str(j[4]), "\t", str(j[5]),
                  "\t", str(j[6]), "\t", str(j[7]), "\t", str(j[8]), "\t", j[9], "\t", j[10], file=fo)
    logging.info("Finish report generation.")


def main():
    parse = argparse.ArgumentParser(description="This pipeline is for dsCir-Seq.")

    parse.add_argument("--sample_file", help="the file with requested format storing whole sample information",
                       required=True)
    parse.add_argument("--output_path", help="the output path", required=True)
    parse.add_argument("--python_command",
                       help="the python command when running python on Linux, such as python or python3", required=True)
    parse.add_argument("--adapter_file", help="the path of adapter file", required=True)
    parse.add_argument("--bwa_index", help="the path of bwa index", required=True)
    parse.add_argument("--reference", "-r", help="Reference Genome Fasta file", required=True)
    parse.add_argument("--GTF_file", "-g", help="GTF file", required=True)
    parse.add_argument("--search_radius", help="Search radius around the position window, default: 20", default=20,
                       type=int)
    parse.add_argument("--gap", help="Gap threshold, default: 5", default=5, type=int)
    parse.add_argument("--overlap_threshold", help="overlap threshold between pair reads, default: 5", default=5,
                       type=int)
    parse.add_argument("--mismatch_threshold", help="Maximum score threshold, default: 6", default=6, type=int)
    parse.add_argument("--top_number", help="output the top # of target sites for visualization, default: 50",
                       default=50, type=int)

    args = parse.parse_args()

    script_path = sys.path[0]
    output_path = os.path.abspath(args.output_path)
    adapter_path = os.path.abspath(args.adapter_file)
    sample_path = os.path.abspath(args.sample_file)
    log_file_path = os.path.join(output_path, 'my.log')

    # define the format of log
    LOG_FORMAT = "[%(asctime)s][%(levelname)s][%(module)s] %(message)s"
    DATE_FORMAT = "%m/%d/%Y %H:%M:%S %p"
    logging.basicConfig(filename=log_file_path, level=logging.DEBUG, format=LOG_FORMAT, datefmt=DATE_FORMAT)

    logging.info("Start analyzing the data of dsCir-Seq, please make sure all parameters are right...")
    parameters = ' '.join(sys.argv[1:])
    logging.info("Parameters: %s" % parameters)

    logging.info("Loading sample file...")
    try:
        dic_sample = import_sample(sample_path)
    except:
        logging.error("The sample file has problems, please make sure its format meets the requirement...")
        sys.exit()

    os.chdir(output_path)
    step1_path = os.path.join(output_path, '01.reads_QC')
    step2_path = os.path.join(output_path, '02.mapping')
    step3_path = os.path.join(output_path, '03.extract_cleavage')
    step4_path = os.path.join(output_path, '04.bamdst')
    step5_path = os.path.join(output_path, '05.report')
    paths = [step1_path, step2_path, step3_path, step4_path, step5_path]
    # check whether these folders have been built
    for path in paths:
        if not os.path.exists(path):
            os.mkdir(path)
        else:
            pass

    logging.info("Step1: reads QC...")
    os.chdir(step1_path)
    try:
        reads_qc(dic_sample, adapter_path)
    except:
        logging.error("Reads QC failed, please check the log file generated in this step.")
        sys.exit()

    logging.info("Step2: reads alignment....")
    os.chdir(step2_path)
    try:
        align_bwa(dic_sample, args.bwa_index, output_path)
    except:
        logging.error("Reads alignment failed, please check the log file generated in this step.")
        sys.exit()

    logging.info("Step3: extract cleavage and report offtarget...")
    os.chdir(step3_path)
    try:
        extract_cleavage(dic_sample, output_path, script_path, args.python_command, args.reference, args.GTF_file,
                         args.search_radius, args.gap, args.overlap_threshold, args.mismatch_threshold, args.top_number)
    except:
        logging.error("Cleavage extraction failed, please check the log file generated in this step.")
        sys.exit()

    logging.info("Step4: summarizing sequencing depth and uniformity...")
    os.chdir(step4_path)
    try:
        bamdst(dic_sample, output_path)
    except:
        logging.error("Computing sequencing information failed, please check the log file generated in this step.")
        sys.exit()

    logging.info("Step5: generate report...")
    os.chdir(step5_path)
    try:
        report(dic_sample, output_path)
    except:
        logging.error("report generation failed, please check the files produced in precedent steps are complete.")


if __name__ == "__main__":
    main()

