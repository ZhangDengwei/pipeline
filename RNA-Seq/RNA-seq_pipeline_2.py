#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="filter method, default:SOAPnuke", choices=["SOAPnuke", "TrimGalore"], default="SOAPnuke")
parser.add_argument("-Q", choices=[1, 2], default= 2, type=int, help="quality system 1:illumina, phred64. 2:sanger, phred33. default: 2")
parser.add_argument("-L", choices=["SE","PE"], required=True, help="single-end or paired-end")
parser.add_argument("-P", "--Pathfile", required=True, help="Sample_id:read1:read2")
parser.add_argument("-x", required=True, help="hisat2 index, argument of hisat2")
parser.add_argument("-G", required=True, help="Genome annotation, xxx.annotation.gtf, argument of stringtie")
parser.add_argument("-s", required=True, help="Genome sequence, xxx.genome.fa, argument of cuffcompare")
parser.add_argument("-O", required=True, help="output path")
args = parser.parse_args()
os.chdir(args.O)
dic_sample = {}
if args.L == "SE":
    with open(args.P,"r") as file_in:
        for line in file_in:
            dic_sample[line.rstrip("\n").split(":")[0]] = line.rstrip("\n").split(":")[1]
elif args.L == "PE":
    with open(args.P,"r") as file_in:
        for line in file_in:
            dic_sample[line.rstrip("\n").split(":")[0]] = [line.rstrip("\n").split(":")[1],line.rstrip("\n").split(":")[2]]

# 1.fastqc
os.chdir(args.O)
os.system("mkdir 1.fastqc && cd  1.fastqc")
for key,value in dic_sample.items():
    os.chdir(args.O + "/1.fastqc")
    if len(list(value)) == 1:
        cmd1 = "fastqc " + value + " -o ./ >log_" + key + ".txt 2>&1"
        os.system(cmd1)
    elif len(list(value)) == 2:
        for x in value:
            cmd2 = "fastqc " + x + " -o ./ >log_" + key + ".txt 2>&1"
            os.system(cmd2)
os.system("echo '1.fastqc: finished!'")

# 2.filter
os.chdir(args.O)
os.system("mkdir 2.filter && cd 2.filter")
for key,value in dic_sample.items():
    os.chdir(args.O + "/2.filter")
    os.system("mkdir " + key)
    os.chdir(args.O + "/2.filter/" + key)
    if args.L == "SE":
        if args.f == "SOAPnuke":
            cmd3 = "SOAPnuke filter -l 15 -q 0.2 -n 0.05 -Q " + args.Q + " -1 " + value + " > log.txt 2>&1"
            os.system(cmd3)
        elif args.f == "TrimGalore":
            if args.Q == 1:
                cmd3 = "trim_galore --phred64 --gzip " + value + " > log.txt 2>&1"
                os.system(cmd3)
            elif args.Q == 2:
                cmd3 = "trim_galore --phred33 --gzip " + value + " > log.txt 2>&1"
                os.system(cmd3)
    elif args.L == "PE":
        if args.f == "SOAPnuke":
            cmd3 = "SOAPnuke filter -l 15 -q 0.2 -n 0.05 -Q " + args.Q + " -1 " + value[0] + " -2 " + value[1] + " > log.txt 2>&1"
            os.system(cmd3)
        elif args.f == "TrimGalore":
            if args.Q == 1:
                cmd3 = "trim_galore --phred64 --gzip " + value[0] + " " + value[1] + " > log.txt 2>&1"
                os.system(cmd3)
            elif args.Q == 2:
                cmd3 = "trim_galore --phred33 --gzip " + value[0] + " " + value[1] + " > log.txt 2>&1"
                os.system(cmd3)
os.system("echo '2.filter: finished!'")

# 3.mapping
os.chdir(args.O)
os.system("mkdir 3.mapping && cd 3.mapping")
for key,value in dic_sample.items():
    os.chdir(args.O + "/3.mapping")
    os.system("mkdir " + key)
    os.chdir(args.O + "/3.mapping/" + key)
    if args.L == "SE":
        if args.Q == 1:
            cmd4 = "hisat2 -k 1 -t -p 12 --phred64 -x " + args.x + " -S accepted_hits.sam -q --no-unal --dta --un-gz unmapped.sam.gz -U " +\
                   "../../2.filter/" + key + "/*.fq.gz" + " > log.txt 2>&1"
            os.system(cmd4)
        elif args.Q == 2:
            cmd4 = "hisat2 -k 1 -t -p 12 --phred33 -x " + args.x + " -S accepted_hits.sam -q --no-unal --dta --un-gz unmapped.sam.gz -U " + \
                   "../../2.filter/" + key + "/*.fq.gz" + " > log.txt 2>&1"
            os.system(cmd4)
    elif args.L == "PE":
        if args.Q == 1:
            cmd4 = "hisat2 -k 1 -t -p 12 --phred64 -x " + args.x + " -S accepted_hits.sam -q --no-unal --dta --un-conc-gz unmapped.sam.gz -1 " +\
                   "../../2.filter/" + key + "/*_1.fq.gz" + " -2 ../../2.filter/" + key + "/*_2.fq.gz" + " > log.txt 2>&1"
            os.system(cmd4)
        elif args.Q == 2:
            cmd4 = "hisat2 -k 1 -t -p 12 --phred33 -x " + args.x + " -S accepted_hits.sam -q --no-unal --dta --un-conc-gz unmapped.sam.gz -1 " + \
                   "../../2.filter/" + key + "/*_1.fq.gz" + " -2 ../../2.filter/" + key + "/*_2.fq.gz" + " > log.txt 2>&1"
            os.system(cmd4)
    os.system("samtools view -Su accepted_hits.sam | samtools sort -o accepted_hits.sorted.bam -O bam")
    cmd5 = "bamCoverage -b accepted_hits.sorted.bam -bs 10 --normalizeUsingRPKM -o " + key + ".bw"
    os.system(cmd5)
    os.system("rm accepted_hits.sam")
os.system("echo '3.mapping: finished!'")

# 4.Gene_Expression
os.chdir(args.O)
os.system("mkdir 4.Gene_Expression && cd 4.Gene_Expression")
for key,value in dic_sample.items():
    os.chdir(args.O + "/4.Gene_Expression")
    os.system("mkdir " + key)
    os.chdir(args.O + "/4.Gene_Expression/" + key)
    cmd6 = "stringtie -t -C cov_refs.gtf.txt -e -B -A gene_abund.tab -p 12 -G " + \
           args.G + " -o Annotated.gtf " + "../../3.mapping/" + key + "/accepted_hits.sorted.bam"
    os.system(cmd6)
    cmd7 = "perl ../../scripts/sumifs.pl gene_abund.tab gene_abund.uniq.tab && sort gene_abund.uniq.tab > gene_abund.uniq.sorted.tab"
    os.system(cmd7)
os.chdir(args.O + "/4.Gene_Expression/" + list(dic_sample.keys())[0])
os.system("prepDE.py -i ../")
os.chdir(args.O + "/4.Gene_Expression/")
cmd8 = "paste "
for key in dic_sample.keys():
    cmd8 = cmd8 + "./" + key + "/gene_abund.uniq.sorted.tab "
os.system(cmd8 + "> pre_merged.gene_fpkm.txt")
for key in dic_sample.keys():
    os.system("echo " + key + " >> sample_list.txt")
os.system("perl ../scripts/merged.pl sample_list.txt pre_merged.gene_fpkm.txt merged.gene_fpkm.txt")
os.system("echo '4.Gene_Expression: finished!'")

# 5.StringTie_assembly
os.chdir(args.O)
os.system("mkdir 5.StringTie_assembly && cd 5.StringTie_assembly")
for key,value in dic_sample.items():
    os.chdir(args.O + "/5.StringTie_assembly")
    os.system("mkdir " + key)
    os.chdir(args.O + "/5.StringTie_assembly/" + key)
    cmd9 = "stringtie -p 12 -m 50 -t -o transcripts.gtf -l " + key + "../../3.mapping/" + key +"/accepted_hits.sorted.bam"
    os.system(cmd9)
os.chdir(args.O + "/5.StringTie_assembly")
for key in dic_sample.keys():
    string = key + "/transcripts.gtf"
    os.system("echo " + string + " >> assembly_GTF_list.merged.txt ")
cmd10 = "cuffcompare -r " + args.G + " -R -s " + args.s + " -i assembly_GTF_list.merged.txt"
os.system(cmd10)
os.system("grep -P -v \"class_code .s.\" cuffcmp.combined.gtf > merged.gtf")
os.system("cut -f9 merged.gtf | cut -f2 -d \";\" | cut -f2 -d \"\\\"\" | sort | uniq > merged.list")
os.system("perl ../scripts/0b_gene_locus.1.pl merged.gtf merged.locus")
cmd11 = "gffread -w isoforms.FASTA -g " + args.s + " merged.gtf"
os.system(cmd11)
os.system("perl ../scripts/FASTA_to_fa.1.pl isoforms.FASTA isoforms.fa")
os.system("perl ../scripts/fa_to_Tab.pl isoforms.fa isoforms.Tab")
os.system("perl ../scripts/1_TCONS_XLOC_seq.pl .merged.gtf isoforms.Tab isoforms.1.Tab")
os.system("echo '5.StringTie_assembly: finished!'")
# 6.Transcript_assembly_expression
os.chdir(args.O)
os.system("mkdir 6.Transcript_assembly_expression && cd 6.Transcript_assembly_expression")
for key,value in dic_sample.items():
    os.system("mkdir " + key)
    os.chdir(args.O + "/6.Transcript_assembly_expression/" + key)
    cmd12 = "stringtie -t -C cov_refs.gtf.txt -e -B -A gene_abund.tab -p 12 -G " + \
           "../../5.StringTie_assembly/merged.gtf" + " -o Annotated.gtf " + "../../3.mapping/" + key + "/accepted_hits.sorted.bam"
    os.system(cmd12)
    cmd13 = "perl ../../scripts/sumifs.pl gene_abund.tab gene_abund.uniq.tab && sort gene_abund.uniq.tab > gene_abund.uniq.sorted.tab"
    os.system(cmd13)
os.chdir(args.O + "/6.Transcript_assembly_expression/" + list(dic_sample.keys())[0])
os.system("prepDE.py -i ../")
os.chdir(args.O + "/6.Transcript_assembly_expression/")
cmd14 = "paste "
for key in dic_sample.keys():
    cmd14 = cmd14 + "./" + key + "/gene_abund.uniq.sorted.tab "
os.system(cmd14 + "> pre_merged.gene_fpkm.txt")
for key in dic_sample.keys():
    os.system("echo " + key + " >> sample_list.txt")
os.system("perl ../scripts/merged.pl sample_list.txt pre_merged.gene_fpkm.txt merged.gene_fpkm.txt")
os.system("echo '6.Transcript_assembly_expression: finished!'")
os.system("echo 'Successfully, all done!'")