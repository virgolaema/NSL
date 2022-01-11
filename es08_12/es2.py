import sys
import os
import subprocess
import shutil
import numpy as np

#performing a grid search to find the best parameters mu, sigma

min_H = 1e5
min_sigma = 0 
min_mu = 0

for sigma in np.linspace(0.5,0.7,10) :  
    for mu in np.linspace(0.7,0.9,10):
        
        print ('')
        print ('Submitting with =====>>>> mu = ', mu, 'sigma =  ', sigma)
        print ('')
    
        input_file = "input.dat"

        template = [
            '10000', 
            '20',
            '0.1',
            '3.2',
            '0',
            '{mu}',
            '{sigma}',
            'grid_search',
            ' ',
            'in >> M',
            'in >> N',
            'in >> x; //starting point',
            'in >> delta; //3.2 gives 50 acceptance',
            'in >> T; (0 for uniform, 1 for gaussian)',
            'in >> mu;',
            'in >> sigma;',
            'in >> filename; (for output)'
        ]
        
        template = '\n'.join(template)
        template = template.format(mu = mu,
                                sigma= sigma
        )
        
        with open( input_file , 'w') as f:
            f.write(template)

        # execute     
        os.system('make fast')

        try_H = np.loadtxt ("grids_final.out", usecols=(0), delimiter=' ', unpack='true')
        if (try_H < min_H):
            min_H = try_H
            print ("new minimum is ", min_H)
            min_sigma = sigma
            min_mu = mu

print ("Optimized valued for mu and sigma after grid search are mu = ", min_mu, " sigma = ", min_sigma, "H min = ", min_H)

