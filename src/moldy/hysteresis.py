'''
Created on 24.8.2015

@author: tohekorh
'''
from misc.lammps_help import  get_lammps_params
from ase.io.trajectory import PickleTrajectory
from structure import get_constraints
from ase.calculators.lammpsrun import LAMMPS
from ase.md.langevin import Langevin
from misc.ui import write_line_own, read_simul_params_file, append_files
from write_structure import saveAndPrint
from misc.solvers import get_R
from ase import units
from ase.visualize import view
import numpy as np
from subprocess import call
import os


taito   =   False
bond    =   1.39695

def hyster_study(edge, folder = None):

    if folder == None:  folder  =   os.getcwd()
        
    print folder
    for fileC in os.listdir(folder):
        if fileC[-6:] == '.simul':
            fileC       =   folder + fileC
            _, length, _, _, v, T, dt, fric, dtheta, \
            thresZ, interval, deltaY, theta, M, edge   =   read_simul_params_file(fileC)
            mdfile      =   fileC[:-6] + '.traj'
            traj        =   PickleTrajectory(mdfile, 'r')
            atoms_init  =   traj[0]
            
            
            constraints, _, twist, rend_b, rend_t  =   get_constraints(atoms_init, edge, \
                                                                   bond, None, key = 'twist_p')
            
            #constraints, _, rend_b, rend_t     =   get_constraints(atoms_init, edge, \
            #                                                       bond, None, key = 'twist_p')
            atoms       =   traj[-1]
            atoms.set_constraint(constraints)
            
            
            
            vels        =   (traj[-2].positions - traj[-1].positions) / (interval * dt)
            atoms.set_velocities(vels)
            
            calc        =   LAMMPS(parameters=get_lammps_params()) 
            atoms.set_calculator(calc)
            
            view(atoms)
            dyn         =   Langevin(atoms, dt*units.fs, T*units.kB, fric)
            twist.set_angle(theta)
            dyn.run(10 * interval)
            view(atoms)
            
            traj_new=   PickleTrajectory(fileC[:-6] + '_hyst.traj', 'w', atoms)
            mdlogf  =   fileC[:-6] + '_hyst.log'
            
            
            
            do_dynamics(mdlogf, atoms, dyn, rend_b, rend_t, v, dt, deltaY, \
                        theta, dtheta, length, thresZ, \
                        interval, traj_new, M, twist)
            
            mdhystf =   fileC[:-6] + '_hyst.traj'
            logfile =   fileC[:-6] + '.log'
            
            append_files(logfile, mdlogf, logfile[:-4] + '_comp.log')
            call(['ase-gui', mdfile, mdhystf, '-o', logfile[:-4] + '_comp.traj'])            

            
            
def do_dynamics(mdlogf, atoms, dyn, rend_b, rend_t, v, dt, deltaY_old, theta, dtheta, \
                L, thres_Z, interval, traj, M, twist, save = True):
    
    write_line_own(mdlogf, '# Hysteresis simul begins \n', 'w')
    
    #kink_vanis  =   False
    
    y0          =   atoms.positions[rend_b, 1]      
    i,m         =   0, 0
    time        =   0.

    while 0 < theta:
        
        theta  -=   dtheta
        twist.set_angle(theta)
        dyn.run(1)
        
        if i % interval == 0:
            
            epot, ekin  =   saveAndPrint(atoms, traj, False)[:2]
            
            deltaY      =   deltaY_old + atoms.positions[rend_b, 1] - y0
            R           =   get_R(L, deltaY)
            hmax        =   np.max(atoms.positions[:,2]) - 3.4 #subsract the height of the plane
            data        =   [time, deltaY, R, theta, hmax, epot, ekin, epot + ekin]
            #data        =   [i * dt, deltaY, R, theta, epot, ekin, epot + ekin]
            
            
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
                  
            print i/interval, theta / (2*np.pi) * 360, R 
            
            if np.max(atoms.positions[:,2]) < thres_Z and m == 0:
                write_line_own(mdlogf, '# Kink vanish! ' +  '\n', 'a')
                print ' kink vanished! '
                m += 1
                #kink_vanis     =   True
                
                
        if i % int(M/100) == 0:    theta / (2 * np.pi) *360
        #if kink_vanis:  m   +=  1
        
        i      +=  1
        time   +=   dt       
            
           
hyster_study('zz')
