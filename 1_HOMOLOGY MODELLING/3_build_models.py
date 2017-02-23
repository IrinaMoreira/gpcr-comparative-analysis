from modeller import *
#Here we import an aditional function, automodel, which takes as inputs
#all the necessary parameters for the creation of models.
from modeller.automodel import *

env = environ()
a = automodel(env, alnfile='d1r-3sn6.ali',
              knowns='3sn6', sequence='d1r',
              assess_methods=(assess.DOPE,
                              assess.GA341))
#The following lines define the starting model and the ending model which 
#roughly translates into the number of models we will be producing
#make() is the function that will start the model generating process
a.starting_model = 1
a.ending_model = 100
a.make()
