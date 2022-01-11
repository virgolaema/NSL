import sys
import os
import subprocess
import shutil

simulations = ["U1","G1","U2","G2"]
for sim in simulations:              
    print ('')
    print ("Execution of ", sim)
        
    os.system('make %s' % sim)
    os.system('make')