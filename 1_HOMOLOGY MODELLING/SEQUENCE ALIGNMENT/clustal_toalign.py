#This script finds .clustal files in its path and converts them to aligned
#tables. 


from glob import glob

all_clustal = glob('*.clustal')

for clustal in all_clustal:
	clustal_name = clustal.split('.')[-2]
	o = open(clustal,'r')
	lines = o.readlines()
	o.close()
	lines = [line.strip() for line in lines]
	lines = [line for line in lines if '*' not in line]
	lines = [line for line in lines if line != '']
	protein_dict = dict()
	protein_list = list()
	for line in lines[1:]:
		data = line.split()
		if data[0] not in protein_dict:
			protein_list.append(data[0])
			protein_dict[data[0]] = data[1]
		else:
			protein_dict[data[0]] += data[1]
	protein_list = sorted(protein_list)
	to_write = ''
	for protein in protein_list:
		pre_write = ',' + protein_dict[protein]
		pre_to_write = ','.join(pre_write)
		to_write += protein + pre_to_write + '\n'
	
	o = open(clustal_name + '.csv','w')
	o.write(to_write)
	o.close()
		
