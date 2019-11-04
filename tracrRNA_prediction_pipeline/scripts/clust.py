#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
# -*- coding: utf-8 -*-
# @Time : 2018/10/27 14:48
# @Author : DW Zhang
# @File : clust.py
# @Software: PyCharm
from sys import argv
import re

script, file1, file2, file3 = argv
with open(file1,"r") as file_in_1: # cluster_25.clstr
    with open(file2,"r") as file_in_2: # cluster_25
        with open(file3,"w") as file_out: # cluster_seq.txt
            dic_clust, dic_seq, clust_list = {}, {}, []
            for line in file_in_1:
                if line[0] == ">":
                    dic_clust[line.rstrip("\n")] = [0,"None"]
                    clust_list.append(line.rstrip("\n"))
                else:
                    dic_clust[clust_list[-1]][0] += 1
                    if line.rstrip("\n")[-1] == "*":
                        dic_clust[clust_list[-1]][1] = re.search(r'>(.+)\.{3}',line).group(1)
            all_clust = file_in_2.readlines()
            for x in range(len(all_clust)):
                if x % 2 == 0:
                    dic_seq[all_clust[x].rstrip("\n")[1:]] = all_clust[x+1].rstrip("\n")
                else:
                    pass
            for x in dic_clust.keys():
                if dic_clust[x][0] >= 20:
                    ID = x + "| " + dic_clust[x][1] + "| cluster:" + str(dic_clust[x][0])
                    seq = dic_seq[dic_clust[x][1]]
                    file_out.writelines([ID,"\n",seq,"\n"])
                else:
                    pass

