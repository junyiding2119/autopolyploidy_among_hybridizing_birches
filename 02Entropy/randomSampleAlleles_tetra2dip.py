import sys
import random
from pathlib import Path
import gzip

vcf_name="07all.bisnp.QC.nocall0.5.maf0.01.LD.vcf.gz"

def sample_alleles(GT_String):
    # 0/0/1/1 subsample as 0/1
    alleles_sampled = sorted(random.sample(GT_String.split("/"),2))
    GT_sampled = '/'.join(alleles_sampled)
    return GT_sampled

## use gzip? 
read_mode = 'r'
if vcf_name.endswith(".gz"):
    read_mode = 'rt'
    ofunc = gzip.open
else:  
    ofunc = open

infile = ofunc(vcf_name, read_mode)
lines = infile.readlines()
infile.close()

header = []
locis = []

for line in lines:
    if line.startswith('#'):
        header.append(line.strip())
        if line.startswith('#CHROM'):
            row = line.strip().split()
            startcol = row.index('FORMAT')
            indnames = [x for x in row[startcol+1:]]
    else:
        locis.append(line)

loci_convert = ['\n'.join(header)]

for line in locis:
    loci = line.strip().split()
    loci[startcol] = "GT"
    for Individual in range(startcol+1,len(loci)):
        GT_String = loci[Individual].split(":")[0]
        if len(GT_String) > 3:
            loci[Individual] = sample_alleles(GT_String)
        elif len(GT_String) == 3:
            loci[Individual] = GT_String
    loci_convert.append('\t'.join(loci))


Path('./08sampled.bisnp.QC.nocall0.5.maf0.01.LD.vcf').write_text('\n'.join(loci_convert))

