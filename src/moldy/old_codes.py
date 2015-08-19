'''
Created on 17.8.2015

@author: tohekorh
'''
import numpy as np
from kc_parallel import KC_potential_p
from filenames import get_fileName
from write_structure import saveAndPrint
from structure import create_stucture, get_topInds, get_rightInds
from ase.calculators.lammpsrun import LAMMPS
from ase.visualize import view
from ase.constraints import FixAtoms, FixedLine
from ase.io.trajectory import PickleTrajectory
from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as mbd
from ase import units


M           =   10000
bond        =   1.39695
h           =   3.3705 
ncores      =   2   
T           =   10 #K
dt, fric    =   2, 0.002
v           =   1.
interval    =   10
dy          =   dt*v/1000.
tau         =   10./fric/5

def shearDyn(width, edge, save = False):
    
    ratio   =   8
    atoms   =   create_stucture(ratio, width, edge, key = 'top')
    
    
    # FIXES
    constraints     =   []
    top, bot        =   get_topInds(atoms)
    rend            =   get_rightInds(atoms, top)
    fix_bot         =   FixAtoms(indices = bot)
    
    view(atoms)
    constraints.append(fix_bot)
    for i in rend:
        constraints.append(FixedLine(i, (1,0,0)))
    
    # KC
    params          =   get_params(atoms)
    params['top_inds']  \
                    =   top
    add_kc          =   KC_potential_p(params)
    constraints.append(add_kc)
    # END FIXES
    
    
    
    # CALCULATOR LAMMPS 
    parameters = {'pair_style':'rebo',
                  'pair_coeff':['* * CH.airebo C H'],
                  'mass'      :['1 12.0', '2 1.0'],
                  'units'     :'metal', 
                  'boundary'  :'p p f'}
    
    calc    =   LAMMPS(parameters=parameters) 
    atoms.set_calculator(calc)
    # END CALCULATOR
    
    
    # TRAJECTORY
    mdfile, mdlogfile, mdrelax  =   get_fileName(edge, width, ratio, v, taito = False)
    
    if save:    traj    =   PickleTrajectory(mdfile, 'w', atoms)
    else:       traj    =   None
    
    #data    =   np.zeros((M/interval, 5))
    
    # RELAX
    atoms.set_constraint(add_kc)
    dyn     =   BFGS(atoms, trajectory = mdrelax)
    dyn.run(fmax=0.05)
    
    # FIX AFTER RELAXATION
    atoms.set_constraint(constraints)
    
    # DYNAMICS
    dyn     =   Langevin(atoms, dt*units.fs, T*units.kB, fric)
    n       =   0
    header  =   '#t [fs], d [Angstrom], epot_tot [eV], ekin_tot [eV], etot_tot [eV] \n'
    log_f   =   open(mdlogfile, 'w')
    log_f.write(header)            
    log_f.close()

    if T != 0:
        # put initial MaxwellBoltzmann velocity distribution
        mbd(atoms, T*units.kB)
    
    
    for i in range(0, M):
        
        if tau < i*dt:
            hw   =   i*dy
            for ind in rend:
                atoms[ind].position[1] += dy 
        
        dyn.run(1)
        
        if i%interval == 0:

            epot, ekin  =   saveAndPrint(atoms, traj, False)[:2]
            
            if T != 0:
                if tau < i*dt:  hw   =   i*dy - tau*v
                else: hw =   0
            else:   hw = i*dy
                
            data        =   [i*dt, hw, epot, ekin, epot + ekin]
            
            if save:
                log_f   =   open(mdlogfile, 'a')
                stringi =   ''
                for k,d in enumerate(data):
                    if k == 0:           
                        stringi += '%.2f ' %d
                    elif k == 1:
                        stringi += '%.6f ' %d
                    else:
                        stringi += '%.12f ' %d
                
                log_f.write(stringi +  '\n')
                log_f.close()
                  

            n += 1
        
        
        if save and T != 0 and i*dt == tau:
            log_f   =   open(mdlogfile, 'a')
            log_f.write('# Thermalization complete. ' +  '\n')
            log_f.close()
            
            
        if 1e2 <= M:    
            if i%(int(M/100)) == 0: print 'ready = %.1f' %(i/(int(M/100))) + '%' 
    

def get_params(atoms):

    params              =   {'bond':bond, 'a':np.sqrt(3)*bond, 'h':h}
    params['ncores']    =   ncores
    params['positions'] =   atoms.positions.copy() 
    params['pbc']       =   atoms.get_pbc()
    params['cell']      =   atoms.get_cell().diagonal()
    params['ia_dist']   =   10
    params['chemical_symbols']  \
                        =   atoms.get_chemical_symbols()
    return params