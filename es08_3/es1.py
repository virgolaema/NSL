import sys
import os
import subprocess
import shutil
import numpy as np


for imT in [4,3,2,1]:  
    print ('')
    print ('Submitting with =====>>>> imT = ', imT)
    print ('')

    input_file = "input.dat"

    template = [
        'timeslices				300',
        'temperature				0.0',
        'imaginaryTimePropagation		{imT}',
        ' ',                                                                            
        'brownianMotionReconstructions           100',
        'delta_translation			1.1',
        'brownianBridgeReconstructions		160',
        'brownianBridgeAttempts			4',
        ' ',                                                                                                        
        'MCSTEPS					4000',
        'equilibration				2000',
        'blocks					20',
         ' ',                                                                                                               
        'histogram_bins				400',
        'histogram_start				-5'
        'histogram_end				5'
        'timeslices_interval_for_averages	120 180'

        'which_wave      0'
        'mu              0.79'
        'sigma           0.61'
    ]
    
    template = '\n'.join(template)
    template = template.format(imT = imT
    )
    
    with open( input_file , 'w') as f:
        f.write(template)

    # execute     
    os.system('make esegui')
    os.system('cp probability.out ./es1/probability%s.out' %imT)