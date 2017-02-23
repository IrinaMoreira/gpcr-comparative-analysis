#This uses a fasta file and converts it to aligned an aligned table format.
#It should be used from the command line as follows:
###python fasta_toalign.py {fasta.file}

import sys

fasta = sys.argv[1]

fasta_name = fasta.split('.')[-2]
o = open(fasta,'r')
lines = o.readlines()
o.close()
lines = [line.strip() for line in lines]
to_write = ''
pre_to_write = ''
for line in lines[1:]:
	pre_to_write += ','.join(line)

to_write += fasta_name + ',,' + pre_to_write

o = open(fasta_name + '.csv','w')
o.write(to_write)
o.close()
	
