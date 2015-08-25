'''
Created on 24.8.2015

@author: tohekorh
'''

def get_lammps_params():
    
    parameters = {'pair_style':'rebo',
              'pair_coeff':['* * CH.airebo C H'],
              'mass'      :['1 12.0', '2 1.0'],
              'units'     :'metal', 
              'boundary'  :'f f f'}
    
    return parameters

def get_simulParams(edge):
    
    T       =   10      # K
    dt      =   1       #[fs]
    fric    =   .002
    
    if edge == 'zz':    thres_Z =   4       # Angst
    elif edge == 'ac':  thres_Z =   5.0     # Angst
    else: raise
    
    v       =   .001    # Angst/fs
    bond    =   1.39695 # Angst
    h       =   3.3705  # Angst 
    
    return T, dt, fric, thres_Z, v, bond, h