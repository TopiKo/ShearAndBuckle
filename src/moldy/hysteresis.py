'''
Created on 24.8.2015

@author: tohekorh
'''
from misc.lammps_help import get_simulParams, get_lammps_params
from ase.io.trajectory import PickleTrajectory
from moldy.structure import get_constraints
from ase.calculators.lammpsrun import LAMMPS
from ase.md.langevin import Langevin
from misc.ui import write_line_own, read_simul_params_file
from write_structure import saveAndPrint
from misc.solvers import get_R, get_dir
from ase import units
import numpy as np
from subprocess import call
import os


taito   =   False
bond    =   1.39695

def hyster_study(edge):
    
    folder  =   os.getcwd()
    print folder
    for fileC in os.listdir(folder):
        if fileC[-6:] == '.simul':
            _, length, _, _, v, T, dt, fric, _, \
            thresZ, interval, deltaY, M, edge   =   read_simul_params_file(fileC)
            mdfile      =   fileC[:-6] + '.traj'
            traj        =   PickleTrajectory(mdfile, 'r')
            atoms_init  =   traj[0]
            
            constraints, _, rend_b, rend_t     =   get_constraints(atoms_init, edge, bond)
            atoms       =   traj[-1]
            atoms.set_constraint(constraints)
            
            
            
            vels        =   (traj[-2].positions - traj[-1].positions) / (interval * dt)
            atoms.set_velocities(vels)
            
            calc        =   LAMMPS(parameters=get_lammps_params()) 
            atoms.set_calculator(calc)
            
            dyn         =   Langevin(atoms, dt*units.fs, T*units.kB, fric)
            dyn.run(10 * interval)
            
            
            traj_new=   PickleTrajectory(fileC[:-6] + '_hyst.traj', 'w', atoms)
            mdlogf  =   fileC[:-6] + '_hyst.log'
            do_dynamics(mdlogf, atoms, dyn, rend_b, rend_t, v, dt, deltaY, length, thresZ, \
                        interval, traj_new, M)
            
            mdhystf =   fileC[:-6] + '_hyst.traj'
            logfile =   fileC[:-6] + '.log'
            
            call(['cat', logfile, mdlogf, '>', logfile])
            call(['ase-gui', mdfile, mdhystf, '-o', logfile[:-4] + '_comb.log'])            
            
            
def do_dynamics(mdlogf, atoms, dyn, rend_b, rend_t, v, dt, deltaY, \
                L, thres_Z, interval, traj, M, save = True):
    
    write_line_own(mdlogf, '# Hysteresis simul begins \n', 'w')
    
    dy          =   v * dt
    kink_vanis  =   False
    dir_vec     =   np.array([0.,-1.,0.])
    i,m,n       =   0, 0, 0
    while m <= int(M / 30):
        
        deltaY  +=  dy*dir_vec[1]
        for ind in rend_b:
            atoms[ind].position[:] += dy*dir_vec 
            
            
        dyn.run(1)
        
        if i%interval == 0:
            
            epot, ekin  =   saveAndPrint(atoms, traj, False)[:2]
            
            R           =   get_R(L, deltaY)
            dir_vec     =   -get_dir(atoms, rend_b, rend_t)
            data        =   [i*dt, deltaY, R, epot, ekin, epot + ekin]
            
            if save:
                stringi =   ''
                for k,d in enumerate(data):
                    if k == 0:           
                        stringi += '%.2f ' %d
                    elif k == 1 or k == 2:
                        stringi += '%.4f ' %d
                    else:
                        stringi += '%.12f ' %d
                write_line_own(mdlogf, stringi +  '\n', 'a')
                  

            n += 1
            
            if np.max(atoms.positions[:,2]) < thres_Z and m == 0:
                write_line_own(mdlogf, '# Kink vanish! ' +  '\n', 'a')
                print ' kink vanished! '
                kink_vanis     =   True
                
                
        if i % int(M/100) == 0:    deltaY
        if kink_vanis:     m   +=  1
        i       +=  1
            
            
           
hyster_study('zz')
