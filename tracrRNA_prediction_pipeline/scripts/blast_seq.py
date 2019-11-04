#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
# -*- coding: utf-8 -*-
# @Time : 2018/10/26 21:42
# @Author : DW Zhang
# @File : blast_seq.py
# @Software: PyCharm
from sys import argv
import re
script, file1, file2, file3 = argv
with open(file1,"r") as file_in_1: # outfmt_xx.txt
    with open(file2,"r") as file_in_2: # mapped_xx.fa
        with open(file3,"w") as file_out: # blast_read.txt
            dic_list, ID_list = {}, []
            all_mapped = file_in_2.readlines()
            for x in range(len(all_mapped)):
                if x % 2 == 0:
                    dic_list[all_mapped[x].rstrip("\n")[1:]] = all_mapped[x+1].rstrip("\n")
                else:
                    pass
            for line in file_in_1:
                line_cont = line.rstrip("\n").split("\t")
                if int(line_cont[3]) >= 8:
                    ID_list.append(line_cont[0])
                else:
                    pass
            for x in set(ID_list):
                file_out.writelines([">",x,"\n",dic_list[x],"\n"])


