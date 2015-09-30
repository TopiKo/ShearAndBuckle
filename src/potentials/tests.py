'''
Created on 21.9.2015

@author: tohekorh
'''

from ase.constraints import FixedLine
from moldy.structure import create_stucture
from LJ_smooth import LJ_potential_smooth
from ase.visualize import view
from ase.calculators.lammpsrun import LAMMPS
from KC_imagBottom import KC_potential_p
import matplotlib.pyplot as plt
import numpy as np
from ase.optimize import BFGS
from moldy.twist_topByRod import trans_atomsKC
from moldy.structure import get_constraints
from analyze.plot import plot_posits

def get_adhesion_energy(atoms, hmax, acc, natoms):
    
    h           =   3.4
    
    def get_epot(z):
        
        new_pos         =   atoms.positions.copy()
        new_pos[:,2]    =   z
        atoms.positions =   new_pos
        
        e   = atoms.get_potential_energy()/natoms
        
        
        print z, e
        return e  

    
    
    # Start to move the top layer in z direction
    zrange  =   np.linspace(h - .8, h + hmax, acc)
    adh_pot =   np.zeros((len(zrange), 2))

    
    
    # Here we lift the top layer:
    for i, z in enumerate(zrange):
        adh_pot[i]  =   [z, get_epot(z)]
        
    return adh_pot

def get_optimal_h(atoms, natoms, dyn = False, show = False):
    
    # This find the optimal h - when the top is sliding:
    
    if not dyn:
        from scipy.optimize import fmin
        pos_init    =   atoms.positions.copy()
        
        def get_epot(z):
            
            new_pos         =   pos_init
            new_pos[:,2]    =   z 
            atoms.positions =   new_pos
            
            e   = atoms.get_potential_energy()/natoms
            #print   z, e
            return  e  
        
        hmin    =   fmin(get_epot, 3.4, disp = 0)
        emin    =   get_epot(hmin)
        
        atoms.positions = pos_init
        
        if show:    print 'optimal height= %.2f and e=%.2f' %(hmin, emin) 
        return emin, hmin
    else:
        dyn         =   BFGS(atoms)
        dyn.run(fmax=0.03)
        e           =   atoms.get_potential_energy()/natoms
        hmin        =   np.average(atoms.positions[:,2])
        return e, hmin

def get_corrugation_energy(atoms, bond, natoms, edge, acc):
    
    if edge == 'zz':
        xset        =   np.linspace(0., np.sqrt(3)*bond, acc) 
    elif edge == 'ac':
        xset        =   np.linspace(0., 4*bond, acc) 
    atoms_init  =   atoms.copy()
    
    
    def get_epot(x):
        
        new_pos         =   atoms_init.positions.copy()
        new_pos[:,0]   +=   x 
        
        atoms.positions =   new_pos
        e, hmin         =   get_optimal_h(atoms, natoms, False)
        print  x,e        
        
        return e, hmin  

    corr_pot            =   np.zeros((acc, 3))
     
    for i, x in enumerate(xset):
        e, hmin         =   get_epot(x)
        corr_pot[i]     =   [x, hmin, e]
    
    corr_pot[:,2]   =   corr_pot[:,2] - np.min(corr_pot[:,2])
    
    return corr_pot


def test_pot():
    
    edge    =   'ac'
    width   =   5
    ratio   =   1
    bond    =   1.39695
    acc     =   100
    
    # Initial atoms graphite:
    #atoms_init              =   make_graphene_slab(a,h,width,length,N, (True, True, False))[3]
    atoms_init                  =   create_stucture(ratio, width, edge, key = 'top', a = bond)[0]

    _, _, _, rend_b, _ =   \
                get_constraints(atoms_init, edge, bond, None, \
                            key = 'twist_p', pot = 'KC')
    trans   =   trans_atomsKC(atoms_init.positions[rend_b], edge, bond)
    atoms_init.translate(trans) 
    
    plot_posits(atoms_init, edge, bond)
    
    natoms  =   0
    for atom in atoms_init:
        if atom.number == 6:
            natoms   +=  1
    
    
    view(atoms_init)
    params  =   {}
    params['positions']         =   atoms_init.positions
    params['chemical_symbols']  =   atoms_init.get_chemical_symbols()   
    params['ia_dist']           =   10
    params['edge']              =   edge
    params['bond']              =   bond    
    params['ncores']            =   2
    
    
    # This controls how far the upper layer is pulled:
    hmax                        =   5
    
    # FIX
    # Certain fixes are imposed. Note the KC and LJ - potential are as constraints
    # added to the atom system. The bottom layer is held fixed:
    
    constraints =   []
    for i in range(len(atoms_init)):
        fix_l   =   FixedLine(i, [0., 0., 1.])
        constraints.append(fix_l)
    
    # The constraints for LJ nad KC - potentials:
    add_KC_p        =   KC_potential_p(params)
    add_LJ          =   LJ_potential_smooth(atoms_init, params['bond'])
    
    
    # Different constraint sets for different calculations:
    constraint_KC_p =   [constraints, add_KC_p]
    constraint_LJ   =   [constraints, add_LJ]
    # END FIX
    
    
    # DEF CALC AND RELAX
    # Define proper calculator rebo when LJ or KC and airebe else:
    # Note! The potential files (CH.airebo etc.) have to be set as Enviromental variables
    # in order for lammps to find them!
    params_rebo         =   {'pair_style':'rebo',
                             'pair_coeff':['* * CH.airebo C H'],
                             'mass'      :['1 12.0', '2 1.0'],
                             'units'     :'metal', 
                             'boundary'  :'f f f'}   
   

    
    constraint_param_sets   = [['rebo_KC_p', constraint_KC_p, params_rebo],
                               ['rebo_LJ', constraint_LJ, params_rebo]] 
    
    data_adh    =   {}
    data_corr   =   {}
    
    indents     =   []
    
    for const_params in constraint_param_sets: 
        
        # Name for the method:
        indent          =   const_params[0]
        indents.append(indent)
        print indent
        
        if indent != 'LJ':
            constraints =   const_params[1]
            parameters  =   const_params[2]
        
            atoms       =   atoms_init.copy()
            init_posits =   atoms.positions.copy()
            
            # Caculator:
            calc        =   LAMMPS(parameters=parameters) #, files=['lammps.data'])
            atoms.set_calculator(calc)
            atoms.set_constraint(constraints)
            

            # SLAB IS RELAXED
            #
            
            print 'Adhesion'
            
            adh                 =   get_adhesion_energy(atoms, hmax, acc, natoms)
            adh[:,1]            =   adh[:,1] - np.min(adh[:,1])
            data_adh[indent]    =   adh - np.min(adh)   
            
            print 'run_moldy'
            atoms.positions     =   init_posits
            
            print 'Corrugation'
            
            data_corr[indent]   =   get_corrugation_energy(atoms, bond, natoms, edge, acc)
            
            '''
            atoms.positions     =   init_posits
            
            data_corr_s[indent] =   get_corrugation_energy_surf(atoms, constraints, \
                                                           bond, bottom, top, indent, acc)
            '''
        #else:
        #    data_adh['LJ']      =   get_adhesion_LJ(data_adh['rebo_KC'][:,0], CperArea) 
    
    #PLOT and save
    
    CperArea    =   4. / (3 * np.sqrt(3) * bond**2)
    ecc         =   0.002843732471143
    sigmacc     =   3.4
    zset        =   np.linspace(3, 10, 100)
    LJ          =   lambda z: 2./5*np.pi*CperArea*ecc*(2*(sigmacc**6/z**5)**2 - 5*(sigmacc**3/z**2)**2)
    LJset       =   LJ(zset) - np.min(LJ(zset))
     
    plt.plot(data_adh['rebo_KC_p'][:,0], data_adh['rebo_KC_p'][:,1])
    plt.plot(data_adh['rebo_LJ'][:,0], data_adh['rebo_LJ'][:,1])
    plt.plot(zset, LJset)
    
    plt.show() 

    plt.plot(data_corr['rebo_KC_p'][:,0], data_corr['rebo_KC_p'][:,2])
    plt.plot(data_corr['rebo_LJ'][:,0], data_corr['rebo_LJ'][:,2])
    
    plt.show() 

#test_pot()