import os
import re
gffs = os.listdir('./annotation')
data = open('data.csv','w')
for gff in gffs:
    with open('./annotation/'+gff, 'r') as f:
        seqid = gff.split('_')[1]
        lines = f.readlines()
        for line in lines:
            if line.startswith('#'):
                continue
            else:
                line = line.strip()
                tmp = line.split('\t')
                REF = tmp[8].split(';')[2][4:]
                ALT = tmp[8].split(';')[3][4:]
                mutaition_type = tmp[8].split(';')[4].split(',')[0][4:]
                pattern = re.compile(r'QHD\d+.\d+')
                try:
                    pro_id = re.findall(pattern, tmp[8])[0]
                except IndexError:
                    pro_id = ""
                data.write(",".join((seqid, tmp[1], tmp[3], tmp[4], REF, ALT, mutaition_type, pro_id))+'\n')
data.close()