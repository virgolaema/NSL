import sys
import os
import subprocess
import shutil

for sim in ["s1","s2", "c1", "c2"]: 
    os.system('cp inputs/input.%s input.dat' %sim)
    os.system('make')
    os.system('mv *.out %s' %sim)