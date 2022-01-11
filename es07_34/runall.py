import sys
import os
import subprocess
import shutil

simulations = ["gas"] #, "liquid", "solid"

for sim in simulations:
    print ("----------------------------")
    print('EXECUTION OF %s' % sim)
    #os.system('cp configs/config.fcc config.0')
    os.system('cp configs/config0.%s config.0')
    os.system('cp configs/config00.%s config.00')
    os.system('cp inputs/input.%s input.dat' % sim)
    os.system('rm *.out')
    os.system('make fast')
    # os.system('mv config.final configs/config0.%s' % sim)
    # os.system('mv old.final configs/config00.%s' % sim)
    os.system('mv *.out 4_%s' % sim)