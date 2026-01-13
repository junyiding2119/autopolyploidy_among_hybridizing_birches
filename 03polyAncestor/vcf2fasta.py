import vcf
import re
import os
import sys
from collections import Counter

iupac_code = {'G': 'G', 'C': 'C', 'T': 'T', 'A': 'A', 'N': '.',
              'R': ['A', 'G'], 'Y': ['C', 'T'], 'M': ['A', 'C'],
              'K': ['G', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
              'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
              'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
            }

def convert_ambigous_base(gt_bases):
    if gt_bases == '*/*' or gt_bases == '*|*' or gt_bases == '*/*/*/*'  or gt_bases == 'N/N' or gt_bases == 'N/N/N/N':
        return 'N'
    elif gt_bases == None:
        return 'N'
    else:
        gt_bases = gt_bases.replace('*', '')
        bases = sorted(re.split(r'/|\|', gt_bases))
        bases = [base for base in bases if base]
        bases_set = set(Counter(bases).keys())
        for convert_item in iupac_code.items():
            if bases_set == set(convert_item[1]):
                return convert_item[0]

    return sys.exit(f"{gt_bases} ERROR\n")

def init(ouf_dir, chrom, pri_chrom, pos, pri_pos, loci, seq_dic):
    end(ouf_dir, pri_chrom, pri_pos, loci, seq_dic)
    # init
    seq_dic = {sample.sample: convert_ambigous_base(sample.gt_bases) for sample in loci.samples}
    pri_chrom = chrom
    return seq_dic, pri_pos, pri_chrom

def end(ouf_dir, chrom, pos, loci, seq_dic):

    if len(seq_dic[ind_for_count_len]) >= 1000: # fasta length
        #print("fasta length",len(seq_dic[ind_for_count_len]))
        t_ouf = open(os.path.join(ouf_dir, f"{chrom}.{pos}.fa"), "w")
        for sample in seq_dic:
            t_ouf.write(">"+sample+"\n"+seq_dic[sample]+"\n")
        t_ouf.close()

def vcf_to_fasta(vcf_reader, ouf_dir):
    start = 0
    seq_dic = {sam_:"" for sam_ in vcf_reader.samples}
    for loci in vcf_reader:

        ref = loci.REF
        alts = loci.ALT
        chrom = loci.CHROM.split('_')[1]
        pos = loci.POS

        if start == 0:
            pri_pos = pos
            pri_chrom = chrom
        #print("start new loc, chr is %s, previous loci is %d, current loci is %d" %(chrom,pri_pos,pos) )

        if chrom != pri_chrom: # split chrom
            seq_dic, pri_pos, pri_chrom = init(ouf_dir, chrom, pri_chrom, pos, pri_pos, loci, seq_dic)
        elif (pos - pri_pos) >= 500: # gap
            seq_dic, pri_pos, pri_chrom = init(ouf_dir, chrom, pri_chrom, pos, pri_pos, loci, seq_dic)
        else:
            for sample in loci.samples:
                seq_dic[sample.sample] += convert_ambigous_base(sample.gt_bases)

        pri_pos = pos
        start = 1

    end(ouf_dir, chrom, pos, loci, seq_dic)

if __name__ == '__main__':
    vcf_filename = sys.argv[1]
    vcf_reader = vcf.Reader(filename=vcf_filename)
    ouf_dir = sys.argv[2]
    ind_for_count_len = vcf_reader.samples[1]

    if not os.path.exists(ouf_dir):
        os.mkdir(ouf_dir)

    vcf_to_fasta(vcf_reader, ouf_dir)