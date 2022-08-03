#!/usr/bin/python

import sys


def process_info(position, line):
    info = line.strip().split('Q')
    for i in range(len(info)):
        info[i] = info[i].split('=')[1]
        info[i] = info[i].split(',')
    chroms = info[0]
    pos = subtract_ref(position, info[1], info[2])
    print(*pos, sep='\t')
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
        print(dif)
        temp_out.append(str(asm-dif))
    return temp_out


results = process_info(sys.argv[1], sys.argv[2])
print(*results, sep="\t")
