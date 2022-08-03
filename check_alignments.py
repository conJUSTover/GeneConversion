#!/usr/bin/python

import sys

def check_differences(line1, line2, header):
    differences = [line1[1]]
    for i in range(len(line1) - 5):
        if line1[i+1] != line2[i+1]:
            differences.append(header[i])
    return differences
#    if len(differences) == 2 and differences[1] == "D10":
#        return([line1[header.index("D10") + 1], line2[header.index("D10") + 1]])
#    else:
#        return



header = []
with open(sys.argv[2], "r") as handle:
    line = handle.read()
    header = line.strip().split()

ref = header.index(sys.argv[3]) + 1
alt = header.index(sys.argv[4]) + 1

bad_pos = []
bed = []
last_pos = ""
good_group = ""
previous_A = []

with open(sys.argv[1], "r") as handle:
    for line in handle:
        line = line.strip().split()
        if line[0] == "1":
            if previous_A:
                if line[ref] == previous_A[ref] and line[alt] == previous_A[alt] and line[ref] != ".;." and line[alt] != ".;.":
                    dif = check_differences(line, previous_A, header)
                    if dif:
                        bad_pos.append(dif)
            previous_A = line

for i in bad_pos:
    for j in range(1,len(i)):
#        print(i[0] + "\t" + i[j])
        print(i[j])
#    if len(i) < 4:
#        print(*i[1:], sep="\t")
#for i in bad_pos:
#    if len(i) == 2:
#        print(i[1])
#print(len(bad_pos))
