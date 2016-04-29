import sys
import os
import re
import logging
import pysam

vcffile = open(sys.argv[1],'r')

for line in vcffile:
    line = line.rstrip()
    record = line.split('\t')

    chr_org = record[0]
    start_org = record[1]
    end_org = record[1]
    rsid_freq = record[2]
    ref_org = record[3]
    alt_org = record[4]

    pattern = re.compile("[A-Z_]+")
    fh = {'AF':'', 'JPT_AF':'', 'JPT_NS':'', 'ASN_AF':'', 'EUR_AF':'', 'AMR_AF':'', 'AFR_AF':'', 'RSID':''}

    if rsid_freq != '.':
        for freq_item in rsid_freq.split(','):
            match = pattern.search(freq_item)
            if match != None: 
                freq_id = match.group(0)
                fh[freq_id] = freq_item.replace(freq_id, "")
            else:
                fh['RSID'] = freq_item

    alt_list = alt_org.split(',')
    for alt in alt_list:
        start = int(start_org)
        ref = ref_org
        ref_len = len(ref)
        alt_len = len(alt)
      
        # SNV
        if alt_len == 1 and ref_len == 1:
            end = int(end_org)

        # deletion
        elif alt_len < ref_len and ref.startswith(alt):
            ref = ref[alt_len:]
            start = start + ref_len - len(ref)
            end = start + len(ref) - 1
            alt = '-'

        # deletion
        elif alt_len > ref_len and alt.startswith(ref):
            alt = alt[ref_len:]
            start = start + alt_len - len(alt)
            end = start
            ref = '-'

        # block substitution 1
        else:
            continue

        print chr_org +"\t"+ str(start-1) +"\t"+ str(end) +"\t"+ ref +"\t"+ alt +"\t"+ fh['RSID'] +';;'+ fh['AF'] +";;"+ fh['JPT_AF'] +';;'+ fh['JPT_NS'] +';;'+ fh['ASN_AF'] +';;'+ fh['EUR_AF'] +';;'+ fh['AMR_AF'] +';;'+ fh['AFR_AF'] + ';;' + ';;'.join(map(str,record[5:12]))
                
####
vcffile.close()


