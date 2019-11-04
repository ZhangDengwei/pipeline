#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
from sys import argv
import os

os.system("echo ==========start at : `date` ==========")
script, pathfile = argv
dic_path = {}
currentpath = os.getcwd()
with open(pathfile,"r") as path:
    for line in path:
        pathinf = line.rstrip("\n").split(":")
        dic_path[pathinf[0]] = pathinf[1]
# 1.SOAPnuke filtersRNA
os.system("mkdir 1.SOAPnuke && cd 1.SOAPnuke")
for key,value in dic_path.items():
    os.chdir(currentpath + "/1.SOAPnuke")
    cmd2 = "mkdir " + key
    os.system(cmd2)
    os.chdir(currentpath + "/1.SOAPnuke/" + key)
    # minimal insert fragment: 30 bp
    cmd4 = "SOAPnuke filtersRNA -f" + value + \
          " -3 AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -5 GAACGACATGGCTACGATCCGACTT -Q 2 -q -z 30 -x " + \
          "Clean_" + key + " > log.txt 2>&1"
    os.system(cmd4)
os.system("echo '1.SOAPnuke: finished!'")
# 2.bwa
os.chdir(currentpath)
os.system("mkdir 2.bwa")
os.chdir(currentpath + "/2.bwa")
os.system("mkdir index && mkdir alignment")
os.chdir(currentpath + "/2.bwa/index")
os.system("cp ../../CRISPR_20kb.fa ./")
os.system("bwa index CRISPR_20kb.fa -p CRISPR_flank")
os.system("cd ../alignment/")
for key in dic_path.keys():
    os.chdir(currentpath + "/2.bwa/alignment")
    os.system("mkdir " + key)
    os.chdir(currentpath + "/2.bwa/alignment/" + key)
    cmd5 = "bwa mem -t 12 ../../index/CRISPR_flank ../../../1.SOAPnuke/" + key + "/Clean_* > output.sam"
    os.system(cmd5)
    cmd6 = "awk -F\"\t\" '{if ($2 != \"4\") print $0}' output.sam > mapped.sam"
    os.system(cmd6)
    os.system("samtools view -bS mapped.sam | samtools sort -O bam > mapped_sorted.bam")
    os.system("samtools index mapped_sorted.bam")
    os.system("python3 ../../../scripts/transfer_fasta.py mapped.sam mapped.fa")
os.system("echo '2.bwa: finished!'")
# 3.blast
os.chdir(currentpath)
os.system("mkdir 3.blast")
os.chdir(currentpath + "/3.blast")
os.system("mkdir index && mkdir alignment")
os.chdir(currentpath + "/3.blast/index")
os.system("cp ../../CRISPR_Array.fa ./")
os.system("makeblastdb -in CRISPR_Array.fa -dbtype nucl && cd ../alignment")
for key in dic_path.keys():
    os.chdir(currentpath + "/3.blast/alignment")
    cmd9 = "blastn -task blastn-short -query ../../2.bwa/alignment/"+ key + "/mapped.fa" +\
          " -db ../index/CRISPR_Array.fa -out outfmt6_" + key +\
            " -outfmt \"6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand\" -num_threads 12 -max_target_seqs 10 -dust no"
    os.system(cmd9)
    cmd10 = "python3 ../../scripts/blast_seq.py outfmt6_" + key + " ../../2.bwa/alignment/" + key + "/mapped.fa blast_read_" + key + ".fa"
    os.system(cmd10)
os.system("echo '3.blast: finished!'")
# 4.cd-hit
os.chdir(currentpath)
os.system("mkdir 4.cdhit")
os.chdir(currentpath + "/4.cdhit")
for key in dic_path.keys():
    os.chdir(currentpath + "/4.cdhit")
    os.system("mkdir " + key)
    os.chdir(currentpath + "/4.cdhit/" + key)
    cmd13 = "cd-hit -i ../../3.blast/alignment/blast_read_" + key + ".fa -o cluster -c 0.8 -aS 0.8 -d 0 -M 6000 > log.txt 2>&1"
    os.system(cmd13)
    os.system("python3 ../../scripts/clust.py cluster.clstr cluster cluster_seq.txt")
    os.system("python3 ../../scripts/sort_cluster.py cluster_seq.txt cluster_seq_sort.txt")
os.system("echo '4.CDHIT: finished!'")
os.system("echo 'Successfully, all done!'")
os.system("echo ==========end at : `date` ==========")

