import sys
import os
import subprocess
import shutil

simulations = ["gas", "liquid", "solid"]

for sim in simulations:
    print ("----------------------------")
    print('EXECUTION OF %s' % sim)
    os.system('make %s' % sim)
    os.system('make fast')
    os.system('mv *output* 1_2_%s' % sim)