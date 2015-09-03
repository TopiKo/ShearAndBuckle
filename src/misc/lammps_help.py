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
    
    if edge == 'zz':    thres_Z =   5.       # Angst
    elif edge == 'ac':  thres_Z =   7.     # Angst
    else: raise
    
    v       =   .0005   # Angst/fs
    # v has to be such that the elastic wave resulting from the bend 
    # can travel the ribbon length twice from end to end before 'too'
    # large deformations occur. If bend that is 5%of the length of the
    # ribbon is not too large .001Angst/fs bending velocity should be
    # slow enough. Sound speend for in-plane waves in graphene > 10km/s 
    # = .1Angst/fs.
    bond    =   1.39695 # Angst
    h       =   3.3705  # Angst 
    
    
    params_dic  =   {}
    params_dic['T'] =   T
    params_dic['dt'] =   dt
    params_dic['fric'] =   fric
    params_dic['edge'] =   edge
    params_dic['thresZ'] =   thres_Z
    params_dic['bond'] =   bond
    params_dic['vmax'] =   v
    params_dic['h'] =   h
    params_dic['tau']   =   10./fric/5
    
    
    
    return T, dt, fric, thres_Z, v, bond, h, params_dic