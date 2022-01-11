import sys
import os
import subprocess
import shutil
import numpy as np


for T in [1.25, 2,3,5,10,20,50]:  
    
    print ('')
    print ('Submitting with =====>>>> T = ', T)
    print ('')

    input_file = "input.dat"

    template = [
        'timeslices				30',
        'temperature				{T}',
        'imaginaryTimePropagation		8',
        ' ',                                                                            
        'brownianMotionReconstructions           0',
        'delta_translation			1.8',
        'brownianBridgeReconstructions		20',
        'brownianBridgeAttempts			4',
        ' ',                                                                                                        
        'MCSTEPS					4000',
        'equilibration				2000',
        'blocks					20',
         ' ',                                                                                                               
        'histogram_bins				400',
        'histogram_start				-5',
        'histogram_end				5',
        'timeslices_interval_for_averages	1 29',

        'which_wave      0',
        'mu              0.79',
        'sigma           0.61',

    ]
    
    template = '\n'.join(template)
    template = template.format(T = T
    )
    
    with open( input_file , 'w') as f:
        f.write(template)

    # execute     
    os.system('make esegui')
    os.system('cp probability.out ./es2/probability%s.out' %T)