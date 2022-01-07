import sys
import os
import subprocess
import shutil

for sim in ["s","c"]: 
    os.system('cp input.%s input.dat' %sim)
    os.system('make')
    os.system('cp *.out %s' %sim)