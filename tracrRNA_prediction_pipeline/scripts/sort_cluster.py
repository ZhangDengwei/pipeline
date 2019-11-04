#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
from sys import argv

script, file1, file2 = argv
with open(file1,"r") as file_in:
    with open(file2,"w") as file_out:
        dic_num = {}
        allcont = file_in.readlines()
        for x in range(len(allcont)):
            if x % 2 == 0:
                name = allcont[x].split("|")[1].strip()
                num = int(allcont[x].split("|")[2].split(":")[1].rstrip("\n"))
                seq = allcont[x+1].rstrip("\n")
                dic_num[name] = [num,seq]
            else:
                pass
        file_out.write("ID\tCluster_number\tSequence\n")
        for x in sorted(dic_num.items(),key=lambda item:item[1][0],reverse=True):
            file_out.writelines([x[0],"\t",str(x[1][0]),"\t",x[1][1],"\n"])

