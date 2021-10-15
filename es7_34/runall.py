import sys
import os
import subprocess
import shutil

for sim in ["liquid", "gas", "solid"]:
    os.system('make %s' %sim)
    os.system('make')
    os.system('cp output* %s' %sim)

    