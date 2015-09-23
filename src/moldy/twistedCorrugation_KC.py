'''
Created on 22.9.2015

@author: tohekorh
'''
import numpy as np
from misc.lammps_help import get_simulParams, get_lammps_params
from potentials.KC_imagBottom import KC_potential_p
from structure import trans_atomsKC
from analyze.plot import plot_posits
from potentials.tests import get_optimal_h
from ase.constraints import FixedLine
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS
from structure import create_stucture
from ase.visualize import view
import matplotlib.pyplot as plt 
import os

path        =   '/space/tohekorh/ShearSlide/files/KC_corrugation/'
    
def corr_KC(width, edge):
    #width   =   5
    #edge    =   'ac'
    params0 =   get_simulParams(edge)[-1]
    bond    =   params0['bond']
    
    
    atoms   =   create_stucture(1, width, edge, key = 'top', a = bond)[0]
    
    atoms.set_cell([40, 40, 20])
    atoms.center()
    atoms.positions[:,2]    =   3.4
    #atoms   =   graphene_nanoribbon(5, 3, type= 'armchair', C_C=bond, saturated = False)
    #atoms.rotate([1,0,0], np.pi/2, rotate_cell = True)
    #atoms.rotate([0,0,1], -np.pi/2, rotate_cell = True)
    #atoms.set_cell([60, 60, 10])
    #atoms.center()
    #del atoms[[48, 51, 52, 55, 56, 59]]
    h_t =   []
    for i in range(len(atoms)):
        if atoms[i].number == 1:
            h_t.append(i)
    
    del atoms[h_t]
    
    
    params  =   {}
    params['positions']         =   atoms.positions
    params['chemical_symbols']  =   atoms.get_chemical_symbols()   
    params['ia_dist']           =   10
    params['edge']              =   edge
    params['bond']              =   bond    
    params['ncores']            =   2
    add_KC                      =   KC_potential_p(params)
    
    
    constraints =   []
    for i in range(len(atoms)):
        fix_l   =   FixedLine(i, [0., 0., 1.])
        constraints.append(fix_l)
    
    constraints.append(add_KC)
    
    lamp_parameters =   get_lammps_params(H=False)
    calc            =   LAMMPS(parameters = lamp_parameters) #, files=['lammps.data'])
    atoms.set_calculator(calc)
    atoms.set_constraint(constraints)
    
    #dyn     =   BFGS(atoms, trajectory = 'test.traj')
    #dyn.run(fmax=0.05)
    
    trans_vec   =   trans_atomsKC(atoms.positions[0], edge, bond)
    atoms.translate(trans_vec)
    
    init_pos    =   atoms.positions.copy()
    middle      =   [np.average(init_pos[:,0]), 
                     np.average(init_pos[:,1]), 
                     np.average(init_pos[:,2])]
     
    
    thetas      =   np.linspace(np.pi/2, np.pi, 91)
    L           =   4*bond
    ds          =   .05
    n           =   int(L/ds)
    
    for i, theta in enumerate(thetas):
        fname   =   path + 'corr_w=%02d_%s_theta=%.2f.data' %(width, edge, theta/(2*np.pi)*360)
        if not os.path.isfile(fname): 
            print 'Calculating w=%i, theta = %.2f' %(width, theta/(2*np.pi)*360)
            
            atoms.positions =   init_pos
            atoms.rotate([0,0,1], theta, center = middle)
            trans_vec       =   np.array([-np.sin(theta), np.cos(theta), 0])
            data            =   np.zeros((n, 3))
            
            for j in range(n):
                atoms.translate(ds*trans_vec)
                emin, hmin  =   get_optimal_h(atoms, len(atoms), dyn = False)
                data[j, :]  =   [j*ds, emin, hmin]   
                #plot_posits(atoms, edge, bond)
            
            header  =   '%s runs along x-dir, angle measured from x-axis, natoms = %i. x (Angs), e (eV/atom), hmin' %(edge, len(atoms))
            np.savetxt(fname, data, header = header)
    

def study_files(width, edge):
    
    theta_dat   =   []
    for filen in os.listdir(path):
        if filen[:12] == 'corr_w=%02d_%s' %(width, edge):
            print filen.split('theta=')[-1].split('.data')[0]
            theta   =   float(filen.split('theta=')[-1].split('.data')[0])
            data    =   np.loadtxt(path + filen)
            data[:,1] -=    np.min(data[:,1]) 
            
            plt.plot(data[:,0],data[:,1])
            plt.show()
            emax    =   np.max(data[:,1])
            theta_dat.append([theta, emax])
    
    # TOPI SORT
    theta_dat   =   np.array(theta_dat) 
    tdat_new    =   np.zeros(theta_dat.shape)
    for i, theta in enumerate(np.sort(theta_dat[:,0])):
        idx = np.where(theta_dat[:,0] == theta)[0][0]
        tdat_new[i,:]   =   [theta, theta_dat[idx, 1]]
        
    return tdat_new

def plot_corr(wedges):
    
    for width, edge in wedges:
        theta_dat   =   study_files(width, edge)
        
        plt.plot(theta_dat[:,0], 1000*theta_dat[:,1], label = 'w=%i, %s' %(width, edge))
        
    plt.xlabel('theta, deg')
    plt.ylabel('corrugation wall, meV/atom')
    plt.legend(frameon = False)
    plt.show()
#corr_KC(5, 'ac')
#plot_corr([[4,'zz'], [5, 'ac']])

#widths  =   [5,7,9,11]
#for width in widths:
#    corr_KC(width, 'ac') 

widths  =   [4,6,8,10]
for width in widths:
    corr_KC(width, 'zz') 
