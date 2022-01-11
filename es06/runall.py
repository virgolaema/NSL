import sys
import os
import subprocess
import shutil

#os.system('touch u.out')
#os.system('touch c.out')
#os.system('touch m.out')
#os.system('touch x.out')

simulations = ["G","G_m","M_m", "M" ] 
for sim in simulations:              
    print ('')
    print ("Execution of ", sim)
    
    os.system ('touch out') #just to have an out file to delete if there's not any
    os.system('make %s' % sim)
    os.system('make')
    os.system('make sposta')