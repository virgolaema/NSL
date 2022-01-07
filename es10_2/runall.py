import os
import shutil

for cores in (1,2,3,4):
    os.system('make compila')
    print('mpirun -np %i ./main %i' %(cores,cores))
    os.system('mpirun -np %s ./main %s'  %(cores, cores))
    os.system('mv final.out results/final_%s.out' %cores)
    os.system('mv meanL.out results/meanL_%s.out' %cores)
    os.system('mv fitness.out results/fitness_%s.out' %cores)
    os.system('mv cities.out results/cities_%s.out' %cores)