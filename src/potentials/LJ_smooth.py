'''
Created on 17.8.2015

@author: tohekorh
'''
import numpy as np

class LJ_potential_smooth:
    
    
    def __init__(self, bond):
        
        self.ecc        =   0.002843732471143
        self.sigmacc    =   3.4
        self.cPerA      =   4. / (3 * np.sqrt(3) * bond**2)
        
    def adjust_positions(self, oldpositions, newpositions):
        pass
    
    
    def adjust_forces(self, posits, forces):
        
        forces[:,2] +=  8*np.pi * self.cPerA * self.ecc * \
                        (self.sigmacc**12 / posits[:,2]**11 \
                       - self.sigmacc**6 / posits[:,2]**5)
   
     
    def adjust_potential_energy(self, posits, energy):
        
        e           =  2. / 5 * np.pi * self.cPerA * self.ecc * \
                        (2 * (self.sigmacc**6 / posits[:,2]**5)**2 \
                       - 5 * (self.sigmacc**3 / posits[:,2]**2)**2) 
   
        return np.sum(e)
    
# TEST
'''
from ase import Atoms
from ase.calculators.lammpsrun import LAMMPS
import matplotlib.pyplot as plt


atoms       =   Atoms('C', positions = [[0.,0.,3.4]])
atoms.cell  =   (5,5,10)
add_LJ      =   LJ_potential_smooth(1.39695)
n           =   1000
# CALCULATOR LAMMPS 
parameters = {'pair_style':'rebo',
              'pair_coeff':['* * CH.airebo C'],
              'mass'      :['1 12.0'],
              'units'     :'metal', 
              'boundary'  :'f f f'}

calc    =   LAMMPS(parameters=parameters) 
atoms.set_calculator(calc)
# END CALCULATOR

atoms.set_constraint(add_LJ)

LJ_e    =   np.zeros(n)
LJ_f    =   np.zeros((n,3))
h       =   np.linspace(2.8, 10, n)

for i in range(n):
    atoms.positions[:,2]    =   h[i]
    LJ_e[i]                 =   atoms.get_potential_energy()
    LJ_f[i,:]               =   atoms.get_forces()


plt.plot(h, LJ_e)
plt.show()
LJ_e1   =   LJ_e.copy()
LJ_e2   =   LJ_e.copy()
LJ_e1   =   np.delete(LJ_e1, -1, 0)
LJ_e2   =   np.delete(LJ_e2, 0, 0)
dh      =   h[1] - h[0]
LJ_ed   =   (LJ_e2 - LJ_e1)/dh



plt.plot(h, LJ_f[:,2])
plt.plot(h[:-1], -LJ_ed)

plt.show()
'''
    
    