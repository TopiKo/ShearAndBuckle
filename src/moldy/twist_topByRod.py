'''
Created on 27.8.2015

@author: tohekorh
'''


import numpy as np
from misc.solvers import get_R
from misc.ui import write_line_own, make_simul_param_file
from misc.lammps_help import get_lammps_params
from filenames import get_fileName
from write_structure import saveAndPrint
from structure import create_stucture, get_constraints
from ase.calculators.lammpsrun import LAMMPS
#from ase.visualize import view
from ase.io.trajectory import PickleTrajectory
from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as mbd
from ase import units


def shearDyn(params_set, save = False):
    
    bond    =   params_set['bond']
    T       =   params_set['T']
    taito   =   params_set['taito']
    dt, fric=   params_set['dt'], params_set['fric']
    tau     =   params_set['tau']
    vmax    =   params_set['vmax']
    vMAX    =   params_set['vMAX']
    thres_Z =   params_set['thresZ']
    width   =   params_set['width']
    ratio   =   params_set['ratio']
    edge    =   params_set['edge']
    
    
    
    atoms, L, W, length_int, b_idxs     =   \
            create_stucture(ratio, width, edge, key = 'top')
    
    mdfile, mdlogfile, mdrelax, simulfile, folder, relaxed \
                        =   get_fileName('LJ', edge + '_twistRod', width, \
                                        length_int, vmax * 1000, int(T), taito)
    
    if relaxed:
        atoms   =   PickleTrajectory(mdrelax, 'r')[-1]
    
    # FIXES
    constraints, add_LJ, twist, rend_b, rend_t =   \
            get_constraints(atoms, edge, bond, b_idxs, key = 'twist_p')
    # END FIXES
    
    
    # CALCULATOR LAMMPS 
    calc    =   LAMMPS(parameters=get_lammps_params()) 
    atoms.set_calculator(calc)
    # END CALCULATOR
    
    # TRAJECTORY
    
    if save:    traj    =   PickleTrajectory(mdfile, 'w', atoms)
    else:       traj    =   None
    
    #data    =   np.zeros((M/interval, 5))
    
    # RELAX
    atoms.set_constraint(add_LJ)
    dyn     =   BFGS(atoms, trajectory = mdrelax)
    dyn.run(fmax=0.05)
    
    dist    =   np.linalg.norm(atoms.positions[rend_b] - atoms.positions[rend_t])
    twist.set_dist(dist)
    # FIX AFTER RELAXATION
    atoms.set_constraint(constraints)
    
    # DYNAMICS
    dyn     =   Langevin(atoms, dt*units.fs, T*units.kB, fric)
    header  =   '#t [fs], shift y [Angstrom], Rad, theta [rad], hmax [A], epot_tot [eV], ekin_tot [eV], etot_tot [eV], F [eV/angst] \n'
    write_line_own(mdlogfile, header, 'w')

    if T != 0:
        # put initial MaxwellBoltzmann velocity distribution
        mbd(atoms, T*units.kB)
    
    y0          =   atoms.positions[rend_b, 1]
    
    kink_formed =   False
    kink_vanished   =   False
    i           =   0
    print 'width = %i, length = %i, v=%.6f' %(width, length_int, vmax)
    
    
    M_therm     =   int(tau / dt)
    dyn.run(M_therm)
    
    M               =   int(2 * L / (np.pi * dt * vmax))
    M_min           =   int(2 * L / (np.pi * dt * vMAX))
    dtheta          =   np.pi / 2 / M
    dtheta_max      =   np.pi / 2 / M_min
    
    interval        =   int( M / 1000 ) 
    theta, time, m  =   0., 0., 0
    i_kink          =   0
    
    
    while 0. <= theta:
        
        if not kink_formed:
            if theta < np.pi/4:
                theta      +=   dtheta_max
            else:
                theta      +=   dtheta
            twist.set_angle(theta)
        else:
            if i_kink / 10 < m: 
                theta      -=   dtheta
                twist.set_angle(theta)
            
        dyn.run(1)
        
        if i % interval == 0:
            epot, ekin  =   saveAndPrint(atoms, traj, False)[:2]
            deltaY      =   atoms.positions[rend_b, 1] - y0
            
            hmax        =   np.max(atoms.positions[:,2]) #substract the height of the plane?
            R           =   get_R(L, deltaY)
            data        =   [time, deltaY, R, theta, hmax, epot, ekin, epot + ekin]
            
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
        
            if thres_Z  <   hmax and not kink_formed:
                idxs    =   np.where(thres_Z < atoms.positions[:,2])
                write_line_own(mdlogfile, '# Kink! at idxs %s' %str(idxs) +  '\n', 'a')
                print ' kink formed! ' + str(i / interval)
                kink_formed     =   True
                i_kink  =   i
                
            if hmax < 3.6 and kink_formed and not kink_vanished:
                write_line_own(mdlogfile, '# Kink vanished! \n', 'a')
                print ' kink vanished '
                kink_vanished   =   True


            print i/interval, theta / (2*np.pi) * 360, R
        
        if kink_formed: m += 1
        
        i      +=   1
        time   +=   dt     

    
    make_simul_param_file(simulfile, W, L, width, length_int, dtheta/dt, dtheta, T, \
                          dt, fric, thres_Z, interval, deltaY, theta, i, edge)
    
    return folder

