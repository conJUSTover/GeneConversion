#!/usr/bin/python

import sys

infile = sys.argv[1]

with open(sys.argv[1], "r") as handle:
    for line in handle:
        if line[0] == '#':
            print(line.strip())
        else:
            #print("not a header lien")
            #print(line)
            line = line.strip().split('\t')
            if "END" in line[7]:
                INFO_string = line[7].split(';')
                END_string = INFO_string[-1].split('=')[-1]
                temp_string = "REF_End=" + str(END_string)
                INFO_string.insert(4, temp_string)
                line[7] = ";".join(INFO_string)
                print(*line, sep="\t")

            else:
                line[7] += ";" + "REF_End=" + str(line[1])
                print(*line, sep="\t")
