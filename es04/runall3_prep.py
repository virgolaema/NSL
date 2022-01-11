import sys
import os
import subprocess
import shutil

simulations = [ 'gas', 'solid', 'liquid'] 
for sim in simulations:              
    print ('')
    print ("Execution of ", sim)
    os.system('make prepare')
    os.system('make %s' % sim)
    os.system('make fcc') #equilibrate
    os.system('make') #compila
    os.system('make esegui')
    os.system('cp config.final configs/config.0_%s' % sim)
    os.system('cp config.00 configs/config.00_%s' % sim)