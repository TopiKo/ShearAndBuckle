'''
Created on 29.9.2015

@author: tohekorh
'''
from structure import create_bend_stucture, strip_Hend, get_idxOfEnds
from filenames import get_fileName
from ase.io.trajectory import PickleTrajectory
from ase.constraints import FixAtoms
import numpy as np
from misc.lammps_help import get_lammps_params
from ase.calculators.lammpsrun import LAMMPS
from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as mbd
from potentials.LJ_smooth import LJ_potential_smooth
from ase.visualize import view
from analyze.plot import plot_posits
from structure import trans_atomsKC
from ase import units
from misc.ui import write_line_own
from write_structure import saveAndPrint
from potentials.KC_imagBottom import KC_potential_p
from misc.ui import make_stick_simul_param_file

def relaxBend(bend, left_idxs, right_idxs, edge, bond, mdrelax):
    
    constraints =   []
    constraints.append(FixAtoms(indices = left_idxs))        
    constraints.append(FixAtoms(indices = right_idxs))        
    #twist       =   twistConst_Rod(bend, 1, edge, bond ,F = 20)
    #twist.set_angle(np.pi/3 + 2./180*np.pi)
    #constraints.append(twist)
    
    add_pot =   LJ_potential_smooth(bend, bond)
    constraints.append(add_pot)
    
    calc    =   LAMMPS(parameters=get_lammps_params()) 
    bend.set_calculator(calc)
    # END CALCULATOR
    
    # RELAX
    bend.set_constraint(constraints)
    dyn     =   BFGS(bend, trajectory = mdrelax)
    dyn.run(fmax=0.05)
    
    

def runAndStudy(params_set, pot_key, save = False):
    
    bond    =   params_set['bond']
    T       =   params_set['T']
    taito   =   params_set['taito']
    dt, fric=   params_set['dt'], params_set['fric']
    tau     =   params_set['tau']
    width   =   params_set['width']
    ratio   =   params_set['ratio']
    edge    =   params_set['edge']
    ncores  =   params_set['ncores']
    LdL_r   =   params_set['LdildeL_ratio']
    
    bend, straight, [matchL_idx, matchR_idx, vec], [L_bend, L_straight], [left_idxs, right_idxs]\
            =   create_bend_stucture(width, ratio, LdL_r, edge, bond)
    
    mdfile, mdlogfile, mdrelax, simulfile, folder, relaxed \
            =   get_fileName(pot_key, edge + '_corrStick', width, \
                             L_bend, L_straight, int(T), taito, key = 'corrStick')
            
    if relaxed:
        bend    =   PickleTrajectory(mdrelax, 'r')[-1]
    else:
        relaxBend(bend, left_idxs, right_idxs, edge, bond, mdrelax)
        bend.set_constraint([])
        
    shift_v =   -straight.positions[matchR_idx] + (bend.positions[matchL_idx] + vec) 
    straight.translate(shift_v)
    
    atoms   =   bend + straight
    cell    =   [1.5*(L_bend + L_straight), L_bend + L_straight, 20]   
    atoms.set_cell(cell)
    atoms.positions[:,2]    =   3.4
    
    trans_vec   =   trans_atomsKC(straight.positions[matchR_idx], edge, bond)
    atoms.translate(trans_vec)
    
    #plot_posits(atoms, edge, bond)
    
    if edge == 'ac':
        nx  =   int((cell[0]/5 - np.min(atoms.positions[:,0]))/(3*bond))  
        ny  =   int((cell[1]/5 - np.min(atoms.positions[:,1]))/(np.sqrt(3)*bond))  
        atoms.translate([nx*3.*bond, ny*np.sqrt(3)*bond, 0])
        width_f =   np.sqrt(3)/2.*bond*(width - 1)
    elif edge == 'zz':
        nx  =   int((cell[0]/5 - np.min(atoms.positions[:,0]))/(np.sqrt(3)*bond))  
        ny  =   int((cell[1]/5 - np.min(atoms.positions[:,1]))/(3*bond))  
        atoms.translate([nx*np.sqrt(3)*bond, ny*3*bond, 0])
        width_f =   (3./2.*width - 1)*bond
    
    cminx, cmaxx    =   strip_Hend(atoms, 'right')
    left_b          =   get_idxOfEnds(atoms, cminx, cmaxx)[0]
    
    # CONSTRAINTS
    constraints =   []
    constraints.append(FixAtoms(indices = left_b))
    
    params      =   {}
    params['positions']         =   atoms.positions
    params['chemical_symbols']  =   atoms.get_chemical_symbols()   
    params['ia_dist']           =   10
    params['edge']              =   edge
    params['bond']              =   bond    
    params['ncores']            =   ncores
    add_pot     =   KC_potential_p(params)
    constraints.append(add_pot)
    atoms.set_constraint(constraints)
    ##
    
    # CALCULATOR
    calc    =   LAMMPS(parameters=get_lammps_params()) 
    atoms.set_calculator(calc)
    ##
    
    # DYNAMICS
    dyn     =   Langevin(atoms, dt*units.fs, T*units.kB, fric)
    header  =   '#t [fs], shift y [Angstrom], Rad, theta [rad], hmax [A], epot_tot [eV], ekin_tot [eV], etot_tot [eV], F [eV/angst] \n'
    write_line_own(mdlogfile, header, 'w')
    traj    =   PickleTrajectory(mdfile, 'w', atoms)
    
    if T != 0:
        # put initial MaxwellBoltzmann velocity distribution
        mbd(atoms, T*units.kB)
    ####
    
    # SIMULATION PARAMS 
    nframes     =   1000
    M           =   int(10*tau/dt)
    interval    =   int(M/nframes)
    thres_cancel=   3*bond
    stick       =   'True'
    xmax_idx    =   np.where(atoms.positions[:,0] == np.max(atoms.positions[:,0]))[0][0]
    r_init      =   atoms.positions[xmax_idx].copy()  
    
    R   =   L_bend/np.pi*3.
    print 'R, width, length bend, theta'
    print R, width_f, L_bend, width_f/(2*R)
    # SIMULATION LOOP
    for i in range(nframes):
        
        print float(i)/nframes*100.
        dyn.run(interval)
        
        epot, ekin  =   saveAndPrint(atoms, traj, False)[:2]
        data        =   [i*interval*dt, epot, ekin, epot + ekin]
        
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
        
        #print np.linalg.norm(atoms.positions[xmax_idx] - r_init)
        if thres_cancel <   np.linalg.norm(atoms.positions[xmax_idx] - r_init):
            stick   =   'false'
            break 
    
    make_stick_simul_param_file(simulfile, width, L_bend, L_straight, T, \
                                dt, fric, interval, M, edge, stick)
    
        
    #plot_posits(atoms, edge, bond)
    #view(atoms)