'''
Created on 12.10.2015

@author: tohekorh
'''



from ase.structure import graphene_nanoribbon
from ase.visualize import view
import numpy as np
from scipy.optimize import fmin
from ase.calculators.lammpsrun import LAMMPS
from scipy.optimize import minimize

bond    =   1.4
#atoms   =   graphene_nanoribbon(1, 1, C_C=bond, sheet=True)
atoms   =   graphene_nanoribbon(1,1, type = 'armchair', C_C=bond, sheet = True)
atoms.rotate([1,0,0], -np.pi/2, rotate_cell=True)
atoms.rotate([0,0,1], -np.pi/2, rotate_cell=True)

del atoms[[0,3]]

atoms.set_cell([[3./2*bond, np.sqrt(3)/2*bond, 0],
                [3./2*bond, -np.sqrt(3)/2*bond, 0],
                [0,0,9]])
atoms.center()
atoms.set_pbc   =   [True, True, False]




def setCell(bond):
    
    cell    =   np.array([[3./2*bond, np.sqrt(3)/2.*bond, 0],
                          [3./2*bond, -np.sqrt(3)/2.*bond, 0],
                          [0,0,9.]])
    atoms.set_cell(cell, scale_atoms = True)
    
def energy(bond):
    setCell(bond)
    e   =   atoms.get_potential_energy()
    print bond, e
    return e


parameters = {'pair_style':'rebo',
              'pair_coeff':['* * CH.airebo C'],
              'mass'      :['1 12.0'],
              'units'     :'metal', 
              'boundary'  :'p p f'}


# CALCULATOR
calc          =   LAMMPS(parameters=parameters) 
atoms.set_calculator(calc)

#(bond*.97, bond*1.03)
#opm_bond    =   fmin(energy, 1.4, xtol = 1e-15, ftol = 1e-15)[0]
res     =   minimize(energy, [bond], args=(), method='L-BFGS-B', 
                     jac=None, bounds=None, tol=None, 
                     callback=None,  
                     options={'disp': None, 'iprint': -1, 
                              'gtol': 1e-015, 'eps': 1e-05, 
                              'maxiter': 150,
                              'ftol':   2.220446049250313e-15, 
                              'maxcor': 10, 
                              'maxfun': 100}) #['x'][0] 'factr': 10000
opm_bond=   res['x'][0]
print 'number of iterations = %i' %res['nit']
            

setCell(opm_bond)
view(atoms)
print opm_bond