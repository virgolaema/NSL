import sys
import os
import subprocess
import shutil

for sim in ["s1"]: #c1 and c2 done s2 
    os.system('cp input.%s input.dat' %sim)
    os.system('make')
    os.system('cp *.out %s' %sim)