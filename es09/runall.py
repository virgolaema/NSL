import sys
import os
import subprocess
import shutil

for sim in ["s1","s2"]: #c1 and c2 done
    os.system('cp input.%s input.dat' %sim)
    os.system('make')
    os.system('cp *.out %s' %sim)