#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
# -*- coding: utf-8 -*-
# @Time : 2018/10/26 21:09
# @Author : DW Zhang
# @File : transfer_fasta.py
# @Software: PyCharm

from sys import argv

script, file1, file2 = argv
with open(file1,"r") as file_in: # mapped.sam
    with open(file2,"w") as file_out: # mapped.fa
        for line in file_in:
            if line[0] != "@":
                line_cont = line.rstrip("\n").split("\t")
                position = line_cont[3] + ":" + str(int(line_cont[3])+len(line_cont[9]))
                ID = ">" + line_cont[0] + "_region_" + position
                seq = line_cont[9]
                file_out.writelines([ID,"\n",seq,"\n"])
            else:
                pass
