#!/usr/bin/python

#TODO create randomized positions
#TODO divide vcf by user-inputted blocks and run GC_check on each block
#TODO Allow for multiple species in d1, d2, p1, p2 upon input

import sys
import argparse

def read_file(d1, p1, p2, d2, vcf):
    #Read file and filter out non-informative SNPs 
    #read header
    GC_SNPs = []
    SNPheader = False


    with open(vcf, "r") as infile:
        #read rest of file 
        for line in infile:
            line = line.strip().split('\t')
            if not SNPheader:
                if line[0] == "#CHROM":
                    SNPheader = process_header(line, d1, p1, p2, d2)
                    print(*SNPheader, sep = '\t')
            else:
                SNP_value = check_SNP(line, SNPheader)
                if SNP_value != 0:
                    GC_SNPs.append([line[0], line[1], str(line[0]) + "_" + str(line[1]), SNP_value])
    return GC_SNPs, SNPheader



def GC_check(SNPs, d1, d2):
    running_count = 0
    max_begin = 0
    min_begin = 0
    last_pos = 0
    d1_count = 0
    d2_count = 0
    all_GC_tally = []
    for line in SNPs:
        SNP_value = line[3]
        if SNP_value > 0:
            if running_count == 0:
                #edit stats
                max_begin = last_pos
                min_begin = line[1]
            #update all remaining stats
            running_count += 1
            last_pos = line[1]
            if SNP_value == 1:
                d1_count += 1
            if SNP_value == 2:
                d2_count += 1
        elif SNP_value == -1:
            if running_count > 1:
                #Process GC site 
                max_end = line[1]
                min_end = last_pos
                if d2_count == 0:
                    donor_dip = d1 + "_donor"
                elif d1_count == 0:
                    donor_dip = d2 + "_donor"
                else: donor_dip = "mixed_donors"
                if max_begin == 0: 
                    note = "terminated_beginning"
                    max_begin = min_begin
                else: note = ""
                all_GC_tally.append([line[0], max_begin, min_begin, min_end, max_end, running_count, d1_count, d2_count, donor_dip, note])
                max_begin = 0
                min_begin = 0
                donor_dip = 0
            d1_count = 0
            d2_count = 0
            max_begin = 0
            min_begin = 0
            running_count = 0
            last_pos = line[1] 
    if running_count > 0:
        #process final string of SNPs
        max_end = last_pos
        min_end = last_pos
        if d2_count == 0:
            donor_dip = d1 + "_donor"
        elif d1_count == 0:
            donor_dip = d2 + "_donor"
        else: donor_dip = "mixed_donors"
        note = "terminated_ending"
        all_GC_tally.append([line[0], max_begin, min_begin, min_end, max_end, running_count, d1_count, d2_count, donor_dip, note])
    return all_GC_tally



def check_SNP(vcf_line, header_list):
    d1_SNPs = [["1", "1", "1", "0"], ["0", "0", "0", "1"]]
    d2_SNPs = [["1", "0", "0", "0"], ["0", "1", "1", "1"]]
    bad_SNPs = [["0", "0", "1", "1"], ["1", "1", "0", "0"]]
    check = [vcf_line[i] for i in header_list]
    if check in d1_SNPs:
        return 1
    if check in d2_SNPs:
        return 2
    if check in bad_SNPs:
        return -1
    else:
        return 0



def process_header(header, d1, p1, p2, d2):
    dip1 = header.index(d1)
    pol1 = header.index(p1)
    pol2 = header.index(p2)
    dip2 = header.index(d2)
    pos = [dip1, pol1, pol2, dip2]
    return pos



def parse_args():
    #Parse Arguments from input
    parser = argparse.ArgumentParser(usage="GC [-h] -c CHROM -d1 DIP_PREFIX -p1 POLY_PREFIX -p2 POLY_PREFIX -d2 DIP_PREFIX -vcf VCF_FILE -o OUTFILE")
    parser.add_argument("d1", help="VCF header column name for Diploid 1.")
    parser.add_argument("p1", help="VCF header column name for Polyploid Subgenome 1.")
    parser.add_argument("p2", help="VCF header column name for Polyploid Subgenome 2.")
    parser.add_argument("d2", help="VCF header column name for Diploid 2.")
    parser.add_argument("vcf", help="Input VCF file. Must contain only one chromosome, and be sorted by position number, and contain header line starting with \"#CHROM\"")
    parser.add_argument("o", help="Outfile to print resulting areas of potential Homoeologous Gene Conversion.")
    args = parser.parse_args()
    print(args.d1)
    print(args.p1)
    print(args.p2)
    print(args.d2)
    print(args.o)
    filtered_vcf, species = read_file(args.d1, args.p1, args.p2, args.d2, args.vcf)
    print(str(len(filtered_vcf)))
    GC_sites = GC_check(filtered_vcf, args.d1, args.d2)
    print_output(GC_sites, args.o)







def print_output(output, outfile):
    with open(outfile, "w+") as handle:
        for item in output:
            items = [str(i) for i in item]
            handle.write("\t".join(items) + "\n")



parse_args()
