'''
Created on 17.8.2015

@author: tohekorh
'''
import numpy as np
from scipy.optimize import fmin

def get_R(L, deltaY):
    
    g = lambda R: (deltaY + R * np.cos(L / R) - R)**2 
    
    if deltaY != 0:
        return fmin(g, L**2 / 2 * deltaY, disp = 0)[0]

    else: return 0

def get_dir(atoms, rend_b, rend_t):
    
    vec     =   atoms.positions[rend_t,:][0] - atoms.positions[rend_b,:][0]
    vec[2]  =   0.
    return vec/np.linalg.norm(vec)

def int_toAngst(thatInt, edge, key = 'width', C_C = 1.42):
    
    if edge == 'zz':
        if key == 'width':
            return thatInt * 3. / 2 * C_C - C_C
        if key == 'length':
            return (thatInt - 1) * np.sqrt(3) / 2 * C_C 
    
    if edge == 'ac':
        if key == 'width':
            return (thatInt -1) * np.sqrt(3) * C_C / 2.  
        if key == 'length':
            return thatInt * 3. / 2 * C_C - C_C 
        


def get_Ld_quess(width, ratio, edge):
    
    bond        =   1.4
    Theta       =   3.14/(6*ratio)
    
    if edge == 'ac': W       =   np.sqrt(3)/2*(width - 1)*bond
    elif edge == 'zz': W    =   3./2*bond*width - bond  
    
    Ld_rat = 44400*Theta**3/W
    return Ld_rat

    
    