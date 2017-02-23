#importing all functions from modeller
from modeller import * 

#creating an environment to which all of our files will be temporarily added
log.verbose()
env = environ() 

###File input
#The sequence database being utilized for our example is the same as the one 
#used for the basic example in the MODELLER website. It contains all
#non-redundant PDB sequences at a 95% threshold. 

#This creates a sequence database with the function sequence_db within our 
#environment - env.
sdb = sequence_db(env)
sdb.read(seq_database_file='pdb_95.pir', seq_database_format='PIR',
         chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)

#Converts our file into a binary file and reads it. This is useful since the
#binary format is processed much faster than the PIR format.
sdb.write(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
          chains_list='ALL')
sdb.read(seq_database_file='pdb_95.bin', seq_database_format='BINARY',
         chains_list='ALL')

#This will be reading our example sequence alignment - d1r.ali
#MODELLER only works with sequences in the .ali format.
#anl.to_profile() will create a new profile for our sequence
#alignment.
aln = alignment(env)
aln.append(file='d1r.ali', alignment_format='PIR', align_codes='ALL')
prf = aln.to_profile()

#We will now be using our sequence database to look for similar sequences
#and writes them in a build_profile.prf output
prf.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
          gap_penalties_1d=(-500, -50), n_prof_iterations=1,
          check_profile=False, max_aln_evalue=0.01)
prf.write(file='build_profile.prf', profile_format='TEXT')

#Converts the profile back to an alignment and builds the profile into 
#build_profile.ali
aln = prf.to_alignment()
aln.write(file='build_profile.ali', alignment_format='PIR')
