#!/usr/bin/python 

import sys

def process_info(line):
    info = line[7].strip().split(';')
    for i in range(len(info)):
        info[i] = info[i].split('=')[1]
        info[i] = info[i].split(',')
    chroms = info[0]
    pos = subtract_ref(line[1], info[1], info[2])
    temp_out = []
    for i in range(len(pos)):
        temp_out.append(";".join([str(chroms[i]), pos[i]]))
    return temp_out

def subtract_ref(pos, asm_end, ref_end):
    temp_out = []
    pos = int(pos)
    for i in range(len(asm_end)):
        asm = int(asm_end[i])
        ref = int(ref_end[i])
        dif = ref - pos
        temp_out.append(str(asm-dif))
    return temp_out


output = []

with open(sys.argv[1], "r") as handle:
    for line in handle:
        if line[0] != "#":
            line = line.strip().split('\t')
            temp_out = []
            REF_ALT = [line[3], line[4]]
            #print(REF_ALT)
            GTs = []
            INFO = process_info(line)
            for i in line[9:]:
                if i == ".":
                    temp_out.append(".;.")
                    GTs.append(".")
                else:
                    temp_out.append(INFO.pop(0))
                    #print(type(int(i)))
                    GTs.append(REF_ALT[int(i)])
#            output.append(temp_out)
            temp_out += temp_out[0].split(';')
            temp_out += temp_out[1].split(';')
            temp_out += temp_out[2].split(';')
            temp_out += GTs
            print(*temp_out, sep="\t")
        elif line[0:6] == "#CHROM":
            line = line.strip().split('\t')
            print(*line[9:], sep="\t")


for i in output:
    print(*i, sep="\t")



