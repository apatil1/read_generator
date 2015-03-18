#!/usr/bin/env python
import sys

fin ="out_bowtie.sam_concordant.txt"
fout = "new_out.txt"
fopen = file (fin, "r")
f = file(fout, "w")

f.write("\t".join(['QNAME', 'RNAME', 'POS', 'XN', 'XM', 'XO', 'XG', 'NM']))

for line in fopen:
    line = line.strip()
    vals = line.split()
    extravals = {}
    for val in vals[11:]:
        k, _, v = val.split(":")
        extravals[k] = v

    output_vals = [
        vals[0], # QNAME
        vals[2], # RNAME
        vals[3], # POS
        extravals["XN"],
        extravals["XM"],
        extravals["XO"],
        extravals["XG"],
        extravals["NM"],
    ]
    
    f.write("\t".join(output_vals))
    f.write("\n")
fopen.close
f.close()
