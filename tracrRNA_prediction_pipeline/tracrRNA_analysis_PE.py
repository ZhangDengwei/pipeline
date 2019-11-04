#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
from sys import argv
import os

os.system("echo ==========start at : `date` ==========")
script, pathfile = argv
dic_path = {}
currentpath = os.getcwd()
with open(pathfile,"r") as path: # pathfile sample:read1:read2
    for line in path:
        pathinf = line.rstrip("\n").split(":")
        dic_path[pathinf[0]] = [pathinf[1],pathinf[2]]
# 1.TrimGalore
os.system("mkdir 1.TrimGalore && cd 1.TrimGalore")
for key,value in dic_path.items():
    os.chdir(currentpath + "/1.TrimGalore")
    os.system("mkdir " + key)
    os.chdir(currentpath + "/1.TrimGalore/" + key)
    # minimal insert fragment: 30 bp
    cmd = "trim_galore --phred33 --gzip -a AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -a2 AAGTCGGATCGTAGCCATGTCGTTC" +\
           " --length 30 --paired " + value[0] + " " + value[1] + " -o ./ > log.txt 2>&1"
    os.system(cmd)
os.system("echo '1.TrimGalore finished!'")
# 2.flash
os.chdir(currentpath)
os.system("mkdir 2.flash")
for key,value in dic_path.items():
    os.chdir(currentpath + "/2.flash")
    os.system("mkdir " + key)
    os.chdir(currentpath + "/2.flash/" + key)
    seqid1 = value[0].split("/")[-1].rstrip(".fq.gz")
    seqid2 = value[1].split("/")[-1].rstrip(".fq.gz")
    cmd = "flash ../../1.TrimGalore/" + key + "/" + seqid1 + "_val_1.fq.gz" + " ../../1.TrimGalore/" + key + "/" + seqid2 + "_val_2.fq.gz > log.txt 2>&1"
    os.system(cmd)
os.system("echo '2.flash finished!'")
# 3.bwa
os.chdir(currentpath)
os.system("mkdir 3.bwa")
os.chdir(currentpath + "/3.bwa")
os.system("mkdir index && mkdir alignment")
os.chdir(currentpath + "/3.bwa/index")
os.system("cp ../../CRISPR_20kb.fa ./")
os.system("bwa index CRISPR_20kb.fa -p CRISPR_flank")
for key,value in dic_path.items():
    os.chdir(currentpath + "/3.bwa/alignment/")
    os.system("mkdir " + key)
    os.chdir(currentpath + "/3.bwa/alignment/" + key)
    cmd1 = "bwa mem -t 12 ../../index/CRISPR_flank ../../../2.flash/" + key + "/out.extendedFrags.fastq > extendedFrags.sam"
    os.system(cmd1)
    cmd2 = "bwa mem -t 12 ../../index/CRISPR_flank ../../../2.flash/" + key + "/out.notCombined_1.fastq > notCombined_1.sam"
    os.system(cmd2)
    cmd3 = "bwa mem -t 12 ../../index/CRISPR_flank ../../../2.flash/" + key + "/out.notCombined_2.fastq > notCombined_2.sam"
    os.system(cmd3)
    os.system("cat extendedFrags.sam notCombined_1.sam notCombined_2.sam > merged.sam")
    cmd4 = "awk -F\"\t\" '{if ($3 == \"CRISPR_20kb\") print $0}' merged.sam > mapped.sam"
    os.system(cmd4)
    os.system("python3 ../../../scripts/uniq_2.py mapped.sam mapped_uniq.sam")
    os.system("python3 ../../../scripts/transfer_fasta.py mapped_uniq.sam mapped_uniq.fa")
    os.system("samtools faidx ../../../CRISPR_20kb.fa")
    os.system("samtools view -ht ../../../CRISPR_20kb.fa.fai mapped_uniq.sam > map_uniq_header.sam")
    os.system("samtools view -bS map_uniq_header.sam | samtools sort -O bam >  map_uniq_sorted.bam && samtools index map_uniq_sorted.bam")
os.system("echo '3.bwa finished!'")
# 4.blast
os.chdir(currentpath)
os.system("mkdir 4.blast")
os.chdir(currentpath + "/4.blast")
os.system("mkdir index && mkdir alignment")
os.chdir(currentpath + "/4.blast/index")
os.system("cp ../../CRISPR_Array.fa ./")
os.system("makeblastdb -in CRISPR_Array.fa -dbtype nucl")
for key in dic_path.keys():
    os.chdir(currentpath + "/4.blast/alignment")
    cmd5 = "blastn -task blastn-short -query ../../3.bwa/alignment/"+ key + "/mapped_uniq.fa" +\
          " -db ../index/CRISPR_Array.fa -out outfmt6_" + key +\
            " -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand\" -num_threads 12 -max_target_seqs 10 -dust no"
    os.system(cmd5)
    cmd6 = "python3 ../../scripts/blast_seq.py outfmt6_" + key + " ../../3.bwa/alignment/" + key + "/mapped_uniq.fa blast_read_" + key + ".fa"
    os.system(cmd6)
os.system("echo '4.blast finished'")
# 5.cdhit
os.chdir(currentpath)
os.system("mkdir 5.cdhit")
os.chdir(currentpath + "/5.cdhit")
for key in dic_path.keys():
    os.chdir(currentpath + "/5.cdhit")
    os.system("mkdir " + key)
    os.chdir(currentpath + "/5.cdhit/" + key)
    cmd7 = "cd-hit -i ../../4.blast/alignment/blast_read_" + key + ".fa -o cluster -c 0.8 -aS 0.8 -d 0 -M 6000 > log.txt 2>&1"
    os.system(cmd7)
    os.system("python3 ../../scripts/clust.py cluster.clstr cluster cluster_seq.txt")
    os.system("python3 ../../scripts/sort_cluster.py cluster_seq.txt cluster_seq_sort.txt")
os.system("echo '5.CDHIT: finished!'")
os.system("echo 'Successfully, all done!'")
os.system("echo ==========end at : `date` ==========")

