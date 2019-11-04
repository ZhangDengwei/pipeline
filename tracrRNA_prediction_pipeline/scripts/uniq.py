#!/hwfssz5/CNGB_WRITE/USER/zhangdengwei/toolkit/python3.7.0/bin/ python3
from sys import argv

script, file1, file2 = argv
with open(file1,"r") as file_in:
    with open(file2,"w") as file_out:
        dicseq = {}
        for line in file_in:
            line_cont = line.rstrip("\n").split("\t")
            if "/" in line_cont[0]:
                readname = line_cont[0].split("/")[0]
                seq = line_cont[9]
                if dicseq.get(readname,0) == 0:
                    dicseq[readname] = line
                else:
                    if dicseq[readname].rstrip("\n").split("\t")[9] in seq or seq in dicseq[readname].rstrip("\n").split("\t")[9]:
                        dicseq[readname] = line if len(seq) > len(dicseq[readname].rstrip("\n").split("\t")[9]) else dicseq[readname]
                    else:
                        file_out.write(line)
            else:
                file_out.write(line)
        for key,values in dicseq.items():
            file_out.write(values)

