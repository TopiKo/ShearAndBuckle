'''
Created on 12.8.2015

@author: tohekorh
'''

import numpy as np
from potentials.LJ_smooth import LJ_potential_smooth
from misc.solvers import get_R
from filenames import get_fileName
from write_structure import saveAndPrint
from structure import create_stucture, get_posInds
from ase.calculators.lammpsrun import LAMMPS
from ase.visualize import view
from ase.constraints import FixedPlane, FixedLine, FixAtoms
from ase.io.trajectory import PickleTrajectory
from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as mbd
from ase import units


M           =   20000
bond        =   1.39695
h           =   3.3705 
ncores      =   2   
T           =   10 #K
dt, fric    =   2, 0.002
v           =   .001 #Angst/fs
# v has to be such that the elastic wave resulting from the bend 
# can travel the ribbon length twice from end to end before 'too'
# large deformations occur. If bend that is 5%of the length of the
# ribbon is not too large .001Angst/fs bending velocity should be
# slow enough.
tau         =   10./fric/5
thres_Z     =   4.

def shearDyn(width, ratio, edge, save = False):
    
    atoms, L, W, length_int  =   create_stucture(ratio, width, edge, key = 'top')
    
    # FIXES
    constraints     =   []
    if edge == 'zz':    fixL    =   np.sqrt(3) * bond * 2.05
    if edge == 'ac':    fixL    =   5 * bond

    rend_b          =   get_posInds(atoms, 'redge')[1]
    lend_s, lend_h  =   get_posInds(atoms, 'ledge', fixL)
    
    #view(atoms)
    for i in rend_b:
        constraints.append(FixedPlane(i, (0,1,0)))
    
    for i in lend_s:
        constraints.append(FixedLine(i, (0,0,1)))
        
    constraints.append(FixAtoms(indices = lend_h))
    
    # KC
    add_LJ          =   LJ_potential_smooth(bond)
    constraints.append(add_LJ)
    # END FIXES
    
    
    
    # CALCULATOR LAMMPS 
    parameters = {'pair_style':'rebo',
                  'pair_coeff':['* * CH.airebo C H'],
                  'mass'      :['1 12.0', '2 1.0'],
                  'units'     :'metal', 
                  'boundary'  :'f f f'}
    
    calc    =   LAMMPS(parameters=parameters) 
    atoms.set_calculator(calc)
    # END CALCULATOR
    
    
    # TRAJECTORY
    mdfile, mdlogfile, mdrelax  =   get_fileName('LJ', edge, width, \
                                                 length_int, \
                                                 taito = False)
    
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
    n       =   0
    header  =   '#t [fs], shift y [Angstrom], Rad, epot_tot [eV], ekin_tot [eV], etot_tot [eV] \n'
    write_line_own(mdlogfile, header, 'w')

    if T != 0:
        # put initial MaxwellBoltzmann velocity distribution
        mbd(atoms, T*units.kB)
    
    deltaY      =   0.
    eta         =   1.1 # ratio between the lengths of the streched and compressed side of the ribbon.
    r           =   L / W # length to width ratio.
    deltaYMax   =   W / 2. * (eta + 1)/(eta - 1) * (1 - np.cos(2 * r * (eta - 1) / (eta + 1))) 
    # Estimate for the required shift
    dy          =   v * dt
    M           =   deltaYMax / dy
    interval    =   int(M / 1000)
    info_line   =   '# width = %.4f [Angstrom], length = %.4f [Angstrom], vel = %.6f [Angstrom/fs]\n' \
                    %(W, L, dy/dt)
    write_line_own(mdlogfile, info_line, 'a')

    kink_formed =   False
    i,m         =   0, 0
    print 'width = %i, length = %i' %(width, length_int)
    while m <= int(M / 50):
        
        if tau < i*dt:
            deltaY  +=  dy
            for ind in rend_b:
                atoms[ind].position[1] += dy 
            
            
        dyn.run(1)
        
        if i%interval == 0:
            
            epot, ekin  =   saveAndPrint(atoms, traj, False)[:2]
            
            R           =   get_R(L, deltaY)
            
                
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
                write_line_own(mdlogfile, stringi +  '\n', 'a')
                  

            n += 1
            
            if thres_Z  <   np.max(atoms.positions[:,2]) and m == 0:
                write_line_own(mdlogfile, '# Kink! ' +  '\n', 'a')
                print ' kink formed! '
                kink_formed     =   True
                
                
                
        if save and T != 0 and i*dt == tau:
            write_line_own(mdlogfile, '# Thermalization complete. ' +  '\n', 'a')
            
        if i % int(M/100) == 0:    print str(i/int(M/100)) + ' ~ % done'
        if kink_formed:     m   +=  1
        i       +=  1

def write_line_own(wfile, line, key):
    log_f       =   open(wfile, key)
    log_f.write(line)            
    log_f.close()
    
    
for rat in range(10,20):
    shearDyn(8, rat, 'zz', True)
    #shearDyn(7, rat, 'ac', True)