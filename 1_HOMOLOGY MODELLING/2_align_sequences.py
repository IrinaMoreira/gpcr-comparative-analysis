from modeller import *

env = environ()
aln = alignment(env)
#We are now creating a model with the model() function and creating a sequence
#alignment which will be used for the modelling.
mdl = model(env, file='3sn6_A', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='3sn6', atom_files='3sn6_A.pdb')
aln.append(file='d1r.ali', align_codes='d1r')
aln.align2d()
aln.write(file='d1r-3sn6.ali', alignment_format='PIR')
aln.write(file='d1r-3sn6.pap', alignment_format='PAP')

