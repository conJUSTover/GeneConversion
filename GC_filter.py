#!/usr/bin/python

#TODO create randomized positions
#TODO Allow for multiple species in d1, d2, p1, p2 upon input

import sys
import argparse

def read_file(d1, p1, p2, d2, vcf, complete):
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
                SNP_value = check_SNP(line, SNPheader, complete)
                if SNP_value != 0:
                    GC_SNPs.append([str(line[0]) + "," + str(line[1]), SNP_value])
    return GC_SNPs, SNPheader



def GC_check(SNPs, d1, d2, homoeo, indels):
    running_count = 0
    max_begin = 0
    min_begin = 0
    last_pos = 0
    d1_count = 0
    d2_count = 0
    nodir_count = 0
    previous = 0
    all_GC_tally = []
    for line in SNPs:
        SNP_value = line[0]
        if SNP_value > 0:
            if running_count == 0:
                #edit stats
                max_begin = last_pos
                min_begin = line[2]
            #update all remaining stats
            running_count += 1
            last_pos = line[2]
            if SNP_value == 1:
                d1_count += 1
            if SNP_value == 2:
                d2_count += 1
            if SNP_value ==3:
                nodir_count += 1
        elif SNP_value == -1:
            if running_count > 0:
                #Process GC site 
                max_end = line[2]
                min_end = last_pos
                if d1_count + d2_count == 0:
                    donor_dip = "Unknown_donor"
                elif d2_count == 0:
                    donor_dip = d1 + "_donor"
                elif d1_count == 0:
                    donor_dip = d2 + "_donor"
                else: donor_dip = "mixed_donors"
                if max_begin == 0: 
                    note = "terminated_beginning"
                    max_begin = min_begin
                else: note = "complete"
                all_GC_tally.append([line[1], max_begin, min_begin, min_end, max_end, running_count, d1_count, d2_count, nodir_count, donor_dip, note])
                max_begin = max_end
#                min_begin = 0
                donor_dip = 0
            elif previous == -1 and homoeo:
                max_end = line[2]
                donor_dip = "homoeoSNP"
                note = "complete"
                max_begin = last_pos
                all_GC_tally.append([line[1], max_begin, min_begin, min_end, max_end, running_count, d1_count, d2_count, nodir_count, donor_dip, note])
            d1_count = 0
            d2_count = 0
            max_begin = 0
            min_begin = 0
            running_count = 0
            last_pos = line[2] 
            nodir_count = 0
            min_end = 0 
        previous = SNP_value
    if running_count > 0:
        #process final string of SNPs
        max_end = last_pos
        min_end = last_pos
        if d1_count + d2_count == 0:
            donor_dip = "Unknown_donor"
        elif d2_count == 0:
            donor_dip = d1 + "_donor"
        elif d1_count == 0:
            donor_dip = d2 + "_donor"
        else: donor_dip = "mixed_donors"
        note = "terminated_ending"
        all_GC_tally.append([line[1], max_begin, min_begin, min_end, max_end, running_count, d1_count, d2_count, nodir_count, donor_dip, note])
    return all_GC_tally



def check_SNP(vcf_line, header_list, complete):
    d1_SNPs = [[1.0, 1.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    d2_SNPs = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 1.0, 1.0]]
    nodirection_SNPs = [[0.0, 1.0, 1.0, 0.0], [1.0, 0.0, 0.0, 1.0]]
    bad_SNPs = [[0.0, 0.0, 1.0, 1.0], [1.0, 1.0, 0.0, 0.0]]
    check = [SNP_freq(vcf_line, i, complete) for i in header_list]
    if check in d1_SNPs:
        return 1
    if check in d2_SNPs:
        return 2
    if check in bad_SNPs:
        return -1
    if check in nodirection_SNPs:
        return 3
    else:
        return 0


def SNP_freq(SNPs, species, complete):
    positions = 0.0
    derived = 0.0
    check = [SNPs[i] for i in species]
    for i in check:
        if i == ".":
            if complete: return None
        elif i == "1":
            positions += 1
            derived += 1
        elif i == "0":
            positions += 1
        else: 
            sys.exit("Unrecognized SNP value: " + i)
    if positions == 0.0: return None #If no positions remain for species group
    return derived / positions

def process_header(header, d1, p1, p2, d2):
    dip1 = [header.index(i) for i in d1.split(',')]
    pol1 = [header.index(i) for i in p1.split(',')]
    pol2 = [header.index(i) for i in p2.split(',')]
    dip2 = [header.index(i) for i in d2.split(',')]
    pos = [dip1, pol1, pol2, dip2]
    return pos


def read_indel(indel_file):
    indel_return = []
    with open(indel_file, "r") as handle:
        for pos in handle:
            temp_pos = pos.strip().split('\t')
            temp_pos[-1] = int(temp_pos[-1])
            indel_return.append(temp_pos)
    return indel_return


def read_bed(bed_file):
    regions = []
    chroms = []
    with open(bed_file, "r") as handle:
        for line in handle:
            line = line.strip().split('\t')
            line[1] = int(line[1])
            line[2] = int(line[2])
            regions.append(line)
            if line[0] not in chroms:
                chroms.append(line[0])
    return regions, chroms


def score_vcf(vcf, columnID, chroms, bed, d1, d2, homoeo, indels):
    final_groups = []
    vcf_short = [[i[1]] + i[columnID].split(',') for i in vcf]
    print(vcf_short[0])
    for i in chroms:
        chrom_vcf = [j for j in vcf_short if j[1] == i]
        temp_bed = [j for j in bed if j[0] == i]
        temp_indel = [j for j in indels if j[0] == i]
        for k in temp_bed:
            temp_vcf = [l for l in chrom_vcf if int(l[2]) >= k[1] and int(l[2]) <= k[2]]
            final_groups += GC_check(temp_vcf, d1, d2, homoeo, temp_indel)
    return final_groups


def parse_args():
    #Parse Arguments from input
    parser = argparse.ArgumentParser(usage="GC [-h] [--no-complete] -d1 DIP_PREFIX -p1 POLY_PREFIX -p2 POLY_PREFIX -d2 DIP_PREFIX -vcf VCF_FILE -bed BEDFILE -o OUTFILE")
    parser.add_argument("-d1", default = None, help="VCF header column name for Diploid 1. Can allow multiple species in comma-delimited list.")
    parser.add_argument("-p1", default = None, help="VCF header column name for Polyploid Subgenome 1. Can allow multiple species in comma-delimited list.")
    parser.add_argument("-p2", default = None, help="VCF header column name for Polyploid Subgenome 2. Can allow multiple species in comma-delimited list.")
    parser.add_argument("-d2", default = None, help="VCF header column name for Diploid 2. Can allow multiple species in comma-delimited list.")
    parser.add_argument("-vcf", default = None, help="Input VCF file. Must contain only one chromosome, and be sorted by position number, and contain header line starting with \"#CHROM\"")
    parser.add_argument("-o", default = None, help="Outfile to print resulting areas of potential Homoeologous Gene Conversion.")
    parser.add_argument("-bed", default = None, help="bed file with regions of the genome to analyze. Regions should be non-overlapping and contain no inversions/translocations within them.")
    parser.add_argument("-indel", default = None, help="Tab-delimited file with chromsome, position, and size of indel in the reference genome relative to ancestral state. Importantly, sign of indel length is size change relative to ancestral state (i.e. alternative allele - reference allele).")
    parser.add_argument("--no-complete", dest='complete', action='store_false', help="When inputting multiple species on a branch tip, this flag allows for missing data in all but one of the species on that tip.")
    parser.add_argument("--homoeoSNPs", dest='homoeoSNPs', action='store_true', help="only identifies homoeoSNPs with zero SNPs in between them. Null distribution of sizes between homoeoSNPs.")
    parser.set_defaults(complete=True)
    parser.set_defaults(homoeoSNPs=False)
    args = parser.parse_args()
    print(str(args.complete))
    print(args.d1)
    print(args.p1)
    print(args.p2)
    print(args.d2)
    print(args.o)
    inbed,chroms = read_bed(args.bed)
    indels = read_indel(args.indel)
#    print(*chroms, sep='\t')
#    print(*inbed, sep='\t')
    filtered_vcf, species = read_file(args.d1, args.p1, args.p2, args.d2, args.vcf, args.complete)
    print(str(len(filtered_vcf)))
#    GC_sites = GC_check(filtered_vcf, args.d1, args.d2)
    GC_sites = score_vcf(filtered_vcf, 0, chroms, inbed, args.d1, args.d2, args.homoeoSNPs, indels)
    print_output(GC_sites, args.o)


def print_output(output, outfile):
    with open(outfile, "w+") as handle:
        for item in output:
            items = [str(i) for i in item]
            handle.write("\t".join(items) + "\n")



parse_args()
