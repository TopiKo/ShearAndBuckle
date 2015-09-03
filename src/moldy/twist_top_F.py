'''
Created on 27.8.2015

@author: tohekorh
'''


import numpy as np
from potentials.LJ_smooth import LJ_potential_smooth
from misc.solvers import get_R, get_dir
from misc.ui import write_line_own, make_simul_param_file, read_simul_params_file
from misc.lammps_help import get_lammps_params, get_simulParams
from filenames import get_fileName
from write_structure import saveAndPrint
from structure import create_stucture, get_posInds, get_constraints
from ase.calculators.lammpsrun import LAMMPS
from ase.visualize import view
from ase.io.trajectory import PickleTrajectory
from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as mbd
from ase import units
import sys

#width, edge, rmin, rmax =  int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])


taito       =   False
width, edge, rmin, rmax =   7, 'ac', 7, 11


T, dt, fric, thres_Z, vmax, bond, h =   get_simulParams(edge)

#bond        =   1.39695
#h           =   3.3705 
#ncores      =   2   
#T           =   10 #K
#dt, fric    =   2, 0.002
#v           =   .001 #Angst/fs


tau         =   10./fric/5
#thres_Z     =   4.

def shearDyn(width, ratio, edge, save = False):
    
    atoms, L, W, length_int, b_idxs     =   create_stucture(ratio, width, edge, key = 'top')
    
    view(atoms)
    # FIXES
    constraints, add_LJ, twist, rend_b, rend_t =   get_constraints(atoms, edge, bond, b_idxs, key = 'twist')
    # END FIXES
    
    
    # CALCULATOR LAMMPS 
    calc    =   LAMMPS(parameters=get_lammps_params()) 
    atoms.set_calculator(calc)
    # END CALCULATOR
    
    
    # TRAJECTORY
    mdfile, mdlogfile, mdrelax, simulfile   =   get_fileName('LJ', edge + '_twist', width, \
                                                 length_int, int(T), taito)
    
    if save:    traj    =   PickleTrajectory(mdfile, 'w', atoms)
    else:       traj    =   None
    
    #data    =   np.zeros((M/interval, 5))
    
    # RELAX
    atoms.set_constraint(add_LJ)
    dyn     =   BFGS(atoms, trajectory = mdrelax)
    dyn.run(fmax=0.05)
    
    # FIX AFTER RELAXATION
    atoms.set_constraint(constraints)
    
    # DYNAMICS
    dyn     =   Langevin(atoms, dt*units.fs, T*units.kB, fric)
    header  =   '#t [fs], shift y [Angstrom], Rad, epot_tot [eV], \
                ekin_tot [eV], etot_tot [eV], F [eV/angst] \n'
    write_line_own(mdlogfile, header, 'w')

    if T != 0:
        # put initial MaxwellBoltzmann velocity distribution
        mbd(atoms, T*units.kB)
    
    y0          =   atoms.positions[rend_b, 1][0]
    
    kink_formed =   False
    dir_vec     =   np.array([0.,1.,0.])
    i,m         =   0, 0
    interval    =   20
    print 'width = %i, length = %i' %(width, length_int)
    
    
    M_therm     =   int(tau / dt)
    
    dyn.run(M_therm)
    F               =   0.
    
    while not kink_formed:
        
        pos_0       =   atoms.positions[rend_b][0]
        F_vec       =   F*np.array([dir_vec[1], -dir_vec[0], 0.])
        twist.set_F(F_vec)
        dyn.run(interval)
        
        i          +=   interval
        epot, ekin  =   saveAndPrint(atoms, traj, False)[:2]
        pos_1       =   atoms.positions[rend_b][0]
        deltaY      =   pos_1[1] - y0
        
        
        v           =   np.sign(pos_1[1] - pos_0[1])*np.linalg.norm(( pos_1 - pos_0 ) / ( interval * dt ))
        R           =   get_R(L, deltaY)
        dir_vec     =   get_dir(atoms, rend_b, rend_t)
        data        =   [i*dt, deltaY, R, epot, ekin, epot + ekin, F_vec[0], F_vec[1], F_vec[2]]
        
        #if      vmax/5 < v:   F   -=  .01
        if    v < vmax/5 and not kink_formed:   F   +=  .01
        
        
        print deltaY, v, F
        
        if save:
            stringi =   ''
            for k,d in enumerate(data):
                if k == 0:           
                    stringi += '%.2f ' %d
                elif k == 1 or k == 2:
                    stringi += '%.4f ' %d
                else:
                    stringi += '%.12f ' %d
            write_line_own(mdlogfile, stringi +  '\n', 'a')
              

        if thres_Z  <   np.max(atoms.positions[:,2]) and m == 0:
            idxs    =   np.where(thres_Z < atoms.positions[:,2])
            write_line_own(mdlogfile, '# Kink! at idxs %s' %str(idxs) +  '\n', 'a')
            print ' kink formed! '
            kink_formed     =   True
            F       =   0
            
                
        if save and T != 0 and i*dt == tau:
            write_line_own(mdlogfile, '# Thermalization complete. ' +  '\n', 'a')
            
        if i % 100 == 0:    print i
        if kink_formed:     m   +=  1
        
    
    make_simul_param_file(simulfile, W, L, width, length_int, v, deltaY/i, T, \
                          dt, fric, thres_Z, interval, deltaY, i, edge)
    

for rat in range(rmin, rmax):
    shearDyn(width, rat, edge, True)
    #shearDyn(7, rat, 'ac', True)