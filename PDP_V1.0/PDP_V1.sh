#!/bin/bash
#****************************information****************************#
echo "Description: PAM Discovery Pipeline(PDP) V1.0"
echo "Anotation:\n\tThis pipeline is only suitable for discovering potential PAM sequence which is result from PAM determination screen!"
#****************************FASTQC****************************#
path=`pwd`
echo "First step: Test sequencing quality using FASTQC"
mkdir Fastqc
cd ./Fastqc
cat ../raw_reads_list | while read line;do
	fastqc ${line} -o ./
done
echo "Fist step: Done!"
#****************************SOAPnuke****************************#
echo "Second step: QC to remove low-quality reads using SOAPnuke"
cd ${path}
mkdir SOAPnuke
cd ./SOAPnuke
row_num=`sed -n '$=' ../raw_reads_list`	# total row numbers
i=0
while read line
do
	arr[i]=$line
	let "i++"
done < ../raw_reads_list
for ((i=0;i<$row_num;i+=2))
do
	SOAPnuke filter -l 15 -q 0.2 -n 0.05 -Q 2 -1 ${arr[$i]} -2 ${arr[$i+1]}
done
echo "Second step: Done!"
#****************************FLASh****************************#
echo "Third step: Merge paired-end reads into one read using FLASh"
cd ${path}
mkdir FLASh
cd ./FLASh
m=0
for var in ${arr[@]}
do
	arr_2[$m]=${var##*/}
	arr_3[$m]=${arr_2[$m]%%.*}
	let "m++"
done
for ((i=0;i<$row_num;i+=2))
do
	cd ${path}/FLASh
	mkdir ${arr_3[$i]} ; cd ${arr_3[$i]}
	flash ../../SOAPnuke/Clean_${arr_2[$i]} ../../SOAPnuke/Clean_${arr_2[$i+1]}
done
echo "Third step: Done!"
#****************************Bowtie2****************************#
echo "Fourth step: Map merged reads to reference sequence using Bowtie2"
cd ${path}
mkdir -p mapping/index
cd ./mapping/index
bowtie2-build -f ../../reference.fa reference
cd ..
mkdir alignment ; cd alignment
for ((i=0;i<$row_num;i+=2))
do
	cd ${path}/mapping/alignment
	mkdir ${arr_3[$i]} ; cd ${arr_3[$i]}
	bowtie2 -k 1 -p 12 -x ../../index/reference --un-gz unmapped.sam.gz -q -U ${path}/FLASh/${arr_3[$i]}/out.extendedFrags.fastq -S output.sam
done
echo "Fourth step: Done!"
#****************************python3****************************#
echo "Fifth step: Extracting 7nt inserted sequence using python3"
cd ${path} 
mkdir PAM_seq ; cd PAM_seq/
cat << EOF > count_PAM.py
#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
from sys import argv
import os
script, file_in, file_out = argv
path = os.getcwd()
filename = path.split("/")[-1]
cont_in = open(file_in,"r",encoding="utf-8") # PAM_sequence.txt
cont_out = open(file_out,"w",encoding="utf-8") # PAM_count.txt
dic_num = {}
for line in cont_in:
    line_cont = line.rstrip().split("\t")
    if line_cont[1] not in dic_num.keys():
        dic_num[line_cont[1]] = 1
    else:
        dic_num[line_cont[1]] += 1
for key,value in dic_num.items():
    cont_out.writelines([key,"\t",str(value),"\t",filename,"\n"])
cont_in.close()
cont_out.close()
EOF
cat << EOF > single_insertion.py
#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
from sys import argv
import re 
script, file_in_1, file_in_2, file_out_1, file_out_2, file_out_3, file_out_4 = argv
cont_in_1 = open(file_in_1,"r",encoding="utf-8") # output.sam
cont_in_2 = open(file_in_2,"r",encoding="utf-8") # ref_seq.txt
cont_out_1 = open(file_out_1,"w",encoding="utf-8") # PAM sequence
cont_out_2 = open(file_out_2,"w",encoding="utf-8") # calling reads
cont_out_3 = open(file_out_3,"w",encoding="utf-8") # uncalling reads
cont_out_4 = open(file_out_4,"w",encoding="utf-8") # site_frequency.txt
cont = cont_in_2.readline()
left_seq = cont.rstrip("\n").split(",")[0]
rigth_seq = cont.rstrip("\n").split(",")[1]
dic_site_frequency = [i for i in range(27)]
for x in range(27):
	dic_site_frequency[x] =  {"A":0,"T":0,"G":0,"C":0}
row_num = 0
for line in cont_in_1:
	row_num += 1
	if row_num <= 3:
		continue
	else:
		line_cont = line.rstrip().split("\t")
		query = line_cont[9]
		if line_cont[2] == "*" or line_cont[5] == "*":
			cont_out_3.write(line)
		else:
			CIGAR = line_cont[5]
			if re.search(r'I',CIGAR): 
				if left_seq[-10:] in query and rigth_seq[:10] in query:
					pattern = re.compile(r""+left_seq[-10:]+"(\w*)"+rigth_seq[:10]+"")
					try:
						PAM_seq = pattern.search(query).group(1)
					except:
						continue
					else:
						if all(["N" not in PAM_seq, len(PAM_seq) == 7]):
							cont_out_1.writelines([line_cont[0],"\t",PAM_seq,"\n"])
							cont_out_2.write(line)
							partsequence = left_seq[-10:] + PAM_seq + rigth_seq[:10]
							i = 0
							for x in partsequence:
								dic_site_frequency[i][x] += 1
								i += 1
						else:
							cont_out_3.write(line)
				else:
					cont_out_3.write(line)
			else:
				cont_out_3.write(line)
cont_out_4.write("position\tA\tG\tC\tT\n")
num = 1
for site in dic_site_frequency:
	cont_out_4.writelines([str(num),"\t"])
	count_total = sum([x for x in site.values()])
	A_frequency = round(site["A"]/count_total,4)
	G_frequency = round(site["G"]/count_total,4)
	C_frequency = round(site["C"]/count_total,4)
	T_frequency = round(site["T"]/count_total,4)
	cont_out_4.writelines([str(A_frequency),"\t",str(G_frequency),"\t",str(C_frequency),"\t",str(T_frequency),"\n"])
	num += 1
cont_in_1.close()
cont_in_2.close()
cont_out_1.close()
cont_out_2.close()
cont_out_3.close()
cont_out_4.close()
EOF
cat << EOF > double_insertion.py
#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
from sys import argv
import re 
script, file_in_1, file_in_2, file_out_1, file_out_2, file_out_3, file_out_4 = argv
cont_in_1 = open(file_in_1,"r",encoding="utf-8") # output.sam
cont_in_2 = open(file_in_2,"r",encoding="utf-8") # ref_seq.txt
cont_out_1 = open(file_out_1,"w",encoding="utf-8") # PAM_sequence.txt
cont_out_2 = open(file_out_2,"w",encoding="utf-8") # calling_reads.txt
cont_out_3 = open(file_out_3,"w",encoding="utf-8") # uncalling_reads.txt
cont_out_4 = open(file_out_4,"w",encoding="utf-8") # site_frequency.txt
cont_1 = cont_in_2.readline()
left_seq_1 = cont_1.rstrip("\n").split(",")[0]
rigth_seq_1= cont_1.rstrip("\n").split(",")[1]
cont_2 = cont_in_2.readline()
left_seq_2 = cont_2.rstrip("\n").split(",")[0]
rigth_seq_2= cont_2.rstrip("\n").split(",")[1]
dic_site_frequency = [i for i in range(27)]
for x in range(27):
	dic_site_frequency[x] =  {"A":0,"T":0,"G":0,"C":0}
def extractPAM(leftseq,rightseq):
	completeseq= leftseq + rightseq
	pattern = re.compile(r""+leftseq[-10:]+"(\w*)"+rightseq[:10]+"")
	try:
		PAM_seq = pattern.search(query).group(1)
	except:
		PAM_seq = ""
		return PAM_seq
	else:
		return PAM_seq
row_num = 0
for line in cont_in_1:
	row_num += 1
	if row_num <= 3:
		continue
	else:
		line_cont = line.rstrip().split("\t")
		query = line_cont[9]
		if line_cont[2] == "*" or line_cont[5] == "*":
			cont_out_3.write(line)
		else:
			CIGAR = line_cont[5]
			if re.search(r'I',CIGAR): 
				if left_seq_1[-10:] in query and rigth_seq_1[:10] in query:
					PAM_1 = extractPAM(left_seq_1,rigth_seq_1)
					if "N" not in PAM_1 and len(PAM_1) == 7:
						cont_out_1.writelines([line_cont[0],"\t",PAM_1,"\n"])
						cont_out_2.write(line)
						partsequence = left_seq_1[-10:] + PAM_1 + rigth_seq_1[:10]
						i = 0
						for x in partsequence:
							dic_site_frequency[i][x] += 1
							i += 1
					elif left_seq_2[-10:] in query and rigth_seq_2[:10] in query:
						PAM_2 = extractPAM(left_seq_2,rigth_seq_2)
						if "N" not in PAM_2 and len(PAM_2) == 7:
							cont_out_1.writelines([line_cont[0],"\t",PAM_2,"\n"])
							cont_out_2.write(line)
							partsequence = left_seq_2[-10:] + PAM_2 + rigth_seq_2[:10]
							i = 0
							for x in partsequence:
								dic_site_frequency[i][x] += 1
								i += 1
						else:
							cont_out_3.write(line)
					else:
						cont_out_3.write(line)
				elif left_seq_2[-10:] in query and rigth_seq_2[:10] in query:
					PAM = extractPAM(left_seq_2,rigth_seq_2)
					if "N" not in PAM and len(PAM) == 7:
						cont_out_1.writelines([line_cont[0],"\t",PAM,"\n"])
						cont_out_2.write(line)
						partsequence = left_seq_2[-10:] + PAM + rigth_seq_2[:10]
						i = 0
						for x in partsequence:
							dic_site_frequency[i][x] += 1
							i += 1
					else:
						cont_out_3.write(line)
					
				else:
					cont_out_3.write(line)
			else:
				cont_out_3.write(line)
num = 1
cont_out_4.write("position\tA\tG\tC\tT\n")
for site in dic_site_frequency:
	cont_out_4.writelines([str(num),"\t"])
	count_total = sum([x for x in site.values()])
	A_frequency = round(site["A"]/count_total,4)
	G_frequency = round(site["G"]/count_total,4)
	C_frequency = round(site["C"]/count_total,4)
	T_frequency = round(site["T"]/count_total,4)
	cont_out_4.writelines([str(A_frequency),"\t",str(G_frequency),"\t",str(C_frequency),"\t",str(T_frequency),"\n"])
	num += 1
cont_in_1.close()
cont_in_2.close()
cont_out_1.close()
cont_out_2.close()
cont_out_3.close()
cont_out_4.close()
EOF
num=`sed -n '$=' ../ref_seq.txt`
if [[ $num -eq 1 ]]
then
	for ((i=0;i<$row_num;i+=2))
	do
		cd ${path}/PAM_seq
		mkdir ${arr_3[$i]} ; cd ${arr_3[$i]}
		python3 ../single_insertion.py ${path}/mapping/alignment/${arr_3[$i]}/output.sam ${path}/ref_seq.txt PAM_sequence.txt calling_reads.txt uncalling_reads.txt site_frequency.txt
		python3 ../count_PAM.py PAM_sequence.txt PAM_count.txt
	done
elif [[ $num -eq 2 ]]
then
	
	for ((i=0;i<$row_num;i+=2))
	do
		cd ${path}/PAM_seq
		mkdir ${arr_3[$i]} ; cd ${arr_3[$i]}
		python3 ../double_insertion.py ${path}/mapping/alignment/${arr_3[$i]}/output.sam ${path}/ref_seq.txt PAM_sequence.txt calling_reads.txt uncalling_reads.txt site_frequency.txt
		python3 ../count_PAM.py PAM_sequence.txt PAM_count.txt
	done
else
	echo "ref_seq.txt : Error!"
fi
echo "Fifth step: Done!"
cd ${path}/PAM_seq
touch merger.txt merger_site_frequency.txt
for ((i=0;i<$row_num;i+=2))
do
	cat merger.txt ./${arr_3[$i]}/PAM_count.txt >> merger.txt
	echo "${arr_3[$i]}" >> merger_site_frequency.txt
	cat merger_site_frequency.txt ./${arr_3[$i]}/site_frequency.txt >> merger_site_frequency.txt
done
cat << EOF > merger.pl
#!/usr/bin/perl -w
die("Argument error!")if(@ARGV!=2);
open(R_1,"< ./\$ARGV[0]")||die("Cannot open \$ARGV[0]: $!\n");
open(W_1,"> ./\$ARGV[1]");
while (\$_=<R_1>) {
	chomp(\$_);
	@array_0=split(/\t/,\$_);
	(\$miRNA,\$Reads,\$Sample)=(\$array_0[0],\$array_0[1],\$array_0[2]);
	\$\$Sample{\$miRNA}=\$Reads;
	\$id{\$miRNA}=\$i;
	unless(exists \$hash_0{\$Sample}){
		\$array_1[++\$#array_1]=\$Sample;
	}
	\$hash_0{\$Sample}=\$i;
}
print W_1 "ID\t";
for(\$i=0;\$i<=\$#array_1;\$i++){
	if(\$i<\$#array_1){
		print W_1 "\$array_1[\$i]\t";
	}
	else{
	 print W_1 "\$array_1[\$i]\n";
	}
}
foreach \$keys(keys %id){
	print W_1 "\$keys\t";
	for(\$i=0;\$i<=\$#array_1;\$i++){
	\$temp=\$array_1[\$i];
	if(exists \$\$temp{\$keys}){
		print W_1 \$\$temp{\$keys};
	}
	else{
	print W_1 "0";
	}
	if(\$i<\$#array_1){
	print W_1 "\t";
	}
	else{
	print W_1 "\n";
	}
	}
}
close R_1;
close W_1; 
EOF
cd ${path}
perl ./PAM_seq/merger.pl ./PAM_seq/merger.txt merger_states.txt
mv ${path}/PAM_seq/merger_site_frequency.txt ./
echo "All done!"

