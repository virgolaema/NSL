import sys
import os
import subprocess
import shutil

simulations = ['gas'] #,"solid", 'liquid']

for sim in simulations:              
    print ('')
    print ("Execution of ", sim)
    os.system('make prepare')
    os.system('make %s' % sim)
    os.system('cp configs/config.0_%s config.0' % sim) #equilibrate
    os.system('cp configs/config.00_%s config.00' % sim) #equilibrate
    os.system('make') #compila
    os.system('make esegui')
    os.system('cp *out* 3_%s' % sim)