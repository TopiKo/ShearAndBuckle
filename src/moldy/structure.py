'''
Created on 13.8.2015

@author: tohekorh
'''

from ase import Atoms, Atom
from ase.visualize import view
from ase.constraints import FixedPlane, FixedLine, FixAtoms
from ase.calculators.neighborlist import NeighborList
from potentials.LJ_smooth import LJ_potential_smooth
from potentials.KC_imagBottom import KC_potential_p

from potentials.twist_const import twist_const_F, twistConst_Rod

import numpy as np

def create_stucture(ratio, width, edge, key = 'top', a = 1.42):
    
    vacuum  =   10
    length  =   get_length(ratio, width, edge)
    length_b=   int(length * .8)*2
    width_b =   width * 4
    
    if edge == 'zz':    
        orig    =   [np.sqrt(3)*a, 2*a] # -> AB-stacking
        b_l, b_w=   4, 3*width
        if width % 2 == 1: raise
    if edge == 'ac':    
        orig    =   [2*a, np.sqrt(3)*a] # -> AB-stacking
        b_l, b_w=   2, 3*width 
        if width % 2 == 0: raise        
    
    bottom  =   graphene_nanoribbon2(length_b, width_b, edge_type=edge, 
                                    saturated=False, 
                                    C_H=1.09,
                                    C_C=a, 
                                    vacuum=2.5, 
                                    sheet=False, 
                                    main_element='C', 
                                    saturate_element='H')
    
    top     =   graphene_nanoribbon2(length, width, edge_type=edge, 
                                    saturated=True, 
                                    C_H=1.09,
                                    C_C=a, 
                                    vacuum=2.5, 
                                    sheet=False, 
                                    main_element='C', 
                                    saturate_element='H')
    
    base    =   graphene_nanoribbon2(b_l, b_w, edge_type=edge, 
                                    saturated=True, 
                                    C_H=1.09,
                                    C_C=a, 
                                    vacuum=2.5, 
                                    sheet=False, 
                                    main_element='C', 
                                    saturate_element='H')
    
    if key == 'top_bottom':
        top.translate([orig[0],orig[1],3.4])
        atoms       =   bottom 
        atoms.extend(top)
        atoms.cell  =   [bottom.cell[0,0], bottom.cell[1,1], 2 * vacuum + 3.4]
        atoms.center()
        atoms.pbc   =   [True, True, False]
        
        #view(atoms)
        return atoms

    if key == 'top':
        atoms       =   top
        L, W        =   top.cell[0,0], top.cell[1,1] 
        atoms.cell  =   [top.cell[0,0]*1.5, 2*top.cell[0,0] + top.cell[1,1], 2 * vacuum]
        atoms.center()
        atoms.positions[:,2] = 3.4
        atoms.pbc   =   [False, False, False]
        return atoms, L, W, length, None
    
    if key == 'rib+base':

        if edge == 'zz':    base.translate([-base.cell[0,0], -base.cell[1,1]/3 ,0])
        if edge == 'ac':    base.translate([-base.cell[0,0], -base.cell[1,1]/3 + np.sqrt(3) / 2. * a,0])

        atoms       =   top + base
        
        rads    =   np.ones(len(atoms))*.3
        nl = NeighborList(rads, self_interaction=False,
                 bothways=True)
        nl.update(atoms)
        
        coll        =   []
        for i, atom in enumerate(atoms):
            if atom.number == 1:
                if 1 < len(nl.get_neighbors(i)[0]):
                    for j in nl.get_neighbors(i)[0]:
                        if atoms.numbers[j] == 1:
                            coll.append(j)
        del atoms[coll]

        atoms_base  =   []
        for i, atom in enumerate(atoms):
            if atom.position[0] < 0:
                atoms_base.append(i)
                
        rads    =   np.ones(len(atoms))*.7
        nl = NeighborList(rads, self_interaction=False,
                 bothways=True)
        nl.update(atoms)

        
        new_pos =   []
        for i, atom in enumerate(atoms):
            if len(nl.get_neighbors(i)[0]) == 2 and atom.number == 6 and i not in atoms_base:
                pos =   np.zeros(3)
                for j in nl.get_neighbors(i)[0]:
                    pos +=  atoms.positions[j] - atoms.positions[i] 
                    if atoms.numbers[j] != 6: raise
                new_pos.append(atoms.positions[i] - pos)
        
        for pos in new_pos:
            atoms   +=  Atom('C', position = pos)
        
        rads    =   np.ones(len(atoms))*.7
        nl = NeighborList(rads, self_interaction=False,
                 bothways=True)
        nl.update(atoms)
        
        h_pos =   []
        for i, atom in enumerate(atoms):
            if len(nl.get_neighbors(i)[0]) == 2 and atom.number == 6:
                pos =   np.zeros(3)
                for j in nl.get_neighbors(i)[0]:
                    pos +=  atoms.positions[j] - atoms.positions[i] 
                    if atoms.numbers[j] != 6: raise
                pos =   pos/np.linalg.norm(pos)
                
                h_pos.append(atoms.positions[i] - pos)
        
        for pos in h_pos:
            atoms   +=  Atom('H', position = pos)
        
        
        
        L, W        =   top.cell[0,0], top.cell[1,1] 
        
        
        
        atoms.cell  =   [top.cell[0,0]*1.5, 2*top.cell[0,0] + top.cell[1,1], 2 * vacuum]
        atoms.center()
        atoms.positions[:,2] = 3.4
        atoms.pbc   =   [False, False, False]
        return atoms, L, W, length, atoms_base

def strip_Hend(atoms, end):
    
    cmaxx   =   0.
    cminx   =   1000.
    for atom in atoms:
        if atom.number == 6:
            if cmaxx < atom.position[0]:
                cmaxx   =   atom.position[0]    
            if atom.position[0] < cminx:
                cminx   =   atom.position[0]
    
    if end == 'right':
        coll_right_h    =   []
        for i, atom in enumerate(atoms):
            if atom.number == 1 and cmaxx < atom.position[0]:
                coll_right_h.append(i)

        del atoms[coll_right_h]

    elif end == 'left':           
        coll_left_h    =   []
        for i, atom in enumerate(atoms):
            if atom.number == 1 and atom.position[0] < cminx:
                coll_left_h.append(i)

        del atoms[coll_left_h] 
        
    return cminx, cmaxx

def get_idxOfEnds(atoms, cminx, cmaxx):
    
    left, right =   [],[]
    
    for i, atom in enumerate(atoms):
        if atom.position[0] < cminx + .02 and atom.number == 6:
            left.append(i)
        if cmaxx - .02 < atom.position[0] and atom.number == 6:
            right.append(i)

    return left, right


def create_bend_stucture(width, lw_ratio, LdildeL_ratio, edge, bond):
    
    
    L_int       =   get_length(lw_ratio, width, edge)
    
    Ldilde_i    =   int(L_int * LdildeL_ratio)
    
    bend        =   graphene_nanoribbon2(L_int, width, edge_type=edge, 
                                         saturated=True, 
                                         C_H=1.09,
                                         C_C=bond, 
                                         vacuum=2.5, 
                                         sheet=False, 
                                         main_element='C', 
                                         saturate_element='H')
    

    cminx, cmaxx    =   strip_Hend(bend, 'right')
    L_bend          =   cmaxx - cminx 
    left_b, right_b =   get_idxOfEnds(bend, cminx, cmaxx)
    matchL_idx      =   right_b[np.where(np.min(bend.positions[right_b, 1]) 
                                       == bend.positions[right_b, 1])[0][0]] 
    

    yav         =   np.average(bend.positions[:,1])    
    xarea       =   L_bend
    R           =   xarea/(np.pi/3)
    for atom in bend:
        old_pos =   atom.position
        angle   =   -(old_pos[0] - cminx)/(xarea)*np.pi/3 + np.pi/2
        rad     =   R + old_pos[1] - yav
        
        new_pos =   [rad*np.cos(angle), rad*np.sin(angle), old_pos[2]]
        atom.position  =   new_pos
    
    straight=   graphene_nanoribbon2(Ldilde_i, width, edge_type=edge, 
                                     saturated=True, 
                                     C_H=1.09,
                                     C_C=bond, 
                                     vacuum=2.5, 
                                     sheet=False, 
                                     main_element='C', 
                                     saturate_element='H')

    
    cminx, cmaxx    =   strip_Hend(straight, 'left')
    L_straight      =   cmaxx - cminx
    left_st         =   get_idxOfEnds(straight, cminx, cmaxx)[0]
    
    vec     =   np.array([bond/2, -np.sqrt(3)/2*bond, 0.])
    if edge == 'zz':
        del straight[left_st]
        cminx, cmaxx=   strip_Hend(straight, 'left')
        left_st     =   get_idxOfEnds(straight, cminx, cmaxx)[0]
        L_straight  =   cmaxx - cminx
        vec         =   np.array([0., -bond, 0.])
    
    matchR_idx      =   left_st[np.where(np.min(straight.positions[left_st, 1]) 
                                       == straight.positions[left_st, 1])[0][0]] 
    
    straight.rotate([0.,0.,1.], -np.pi/3)
    
    
    
    #cellxt  =   1.5*(L_straight + L_bend)
    #cellyt  =   2*(np.max(tot_str.positions[:,1]) - np.min(tot_str.positions[:,1]))
    #cellzt  =   20

    cellxb  =   1.5*L_bend
    cellyb  =   2*(np.max(bend.positions[:,1]) - np.min(bend.positions[:,1]))
    cellzb  =   20

    bend.set_cell([cellxb, cellyb, cellzb])
    bend.center()
    bend.positions[:,2] =   3.4

    
    #tot_str.set_cell([cellxt, cellyt, cellzt])
    #tot_str.center()
    
    #tot_str.positions[:,2] =   3.4
    
    return bend, straight, [matchL_idx, matchR_idx, vec], [L_bend, L_straight], [left_b, right_b]
    

def graphene_nanoribbon2(length, width, edge_type='zz', saturated=False, C_H=1.09,
                        C_C=1.42, vacuum=2.5, magnetic=None, initial_mag=1.12,
                        sheet=False, main_element='C', saturate_element='H'):
    """Create a graphene nanoribbon.

    Creates a graphene nanoribbon in the x-z plane, with the nanoribbon
    running along the z axis.

    Parameters:

    n: int
        The width of the nanoribbon.
    m: int
        The length of the nanoribbon.
    type: str
        The orientation of the ribbon.  Must be either 'zigzag'
        or 'armchair'.
    saturated: bool
        If true, hydrogen atoms are placed along the edge.
    C_H: float
        Carbon-hydrogen bond length.  Default: 1.09 Angstrom.
    C_C: float
        Carbon-carbon bond length.  Default: 1.42 Angstrom.
    vacuum: float
        Amount of vacuum added to both sides.  Default 2.5 Angstrom.
    magnetic: bool
        Make the edges magnetic.
    initial_mag: float
        Magnitude of magnetic moment if magnetic=True.
    sheet: bool
        If true, make an infinite sheet instead of a ribbon.
    """

    
    b = np.sqrt(3) * C_C
    ac_unit    =   Atoms(main_element + '2',
                          pbc=(1, 1, 0),
                          cell=[3 * C_C, b, 2 * vacuum])
    
    ac_unit.positions  = [[0, 0, 0],
                          [C_C, 0., 0.]]
    
    ac_sat_unit=   Atoms(saturate_element + '2',
                          pbc=(1, 1, 0),
                          cell=[3 * C_C, b, 2 * vacuum])
    
    ac_sat_unit.positions   =   [[0, 0, 0],
                                 [C_C + C_H, 0, 0]]
    
    
    def slab(m, n):
        atoms = Atoms()
            
        for i in range(n):
            if i % 2 == 0:
                layer   =   ac_unit.repeat((m, 1, 1))
            if i % 2 == 1:
                layer   =   ac_unit.repeat((m, 1, 1))
                layer.positions[:, 0] += 3./2 * C_C
                
            if saturated:
                if i == 0:  
                    sat =   ac_sat_unit.repeat((m,1,1))
                    sat.positions[:,1]  -= np.sqrt(3.) / 2 * C_H
                    sat.positions[:,0]  -= 1. / 2 * C_H
                    layer +=    sat
                elif i == n - 1:    
                    sat =   ac_sat_unit.repeat((m,1,1))
                    sat.positions[:,1]  += np.sqrt(3.) / 2 * C_H 
                    sat.positions[:,0]  -= 1. / 2 * C_H - 3./2 * C_C * (i % 2)
                    layer +=    sat
            
            layer.positions[:, 1] += b / 2 * i
            atoms += layer
        
        xmax    =   np.max(atoms.positions[:,0])
        for atom in atoms:
            if xmax - C_C/2. < atom.position[0]:
                atom.position[0] -=     m*3*C_C 
        
        if saturated:
            xmax    =   np.max(atoms.positions[:,0])
            xmin    =   np.min(atoms.positions[:,0])
            for atom in atoms:
                posit   =   atom.position
                if xmax - C_C/6. < posit[0] and atom.number == 6:
                    h_posit =  [posit[0] + C_H, posit[1], posit[2]]
                    atoms  +=  Atom('H', position = h_posit)  
                if  posit[0] < xmin + C_C/6. and atom.number == 6:
                    h_posit =  [posit[0] - C_H, posit[1], posit[2]]
                    atoms  +=  Atom('H', position = h_posit)  
        
        atoms.cell = [m * 3 * C_C, n * b/2, 2 * vacuum]
        atoms.center()
        return atoms
        
    if edge_type == 'ac':
        return slab(length, width)
    
    elif edge_type == 'zz':
        
        atoms       =   slab(width / 2 + width % 2, length)
        atoms.rotate('-z', np.pi / 2, rotate_cell = False)
        atoms.cell  =   [length* b/2 , (width / 2 + width % 2) * 3 * C_C , 2 * vacuum]
        atoms.center()
        
        if width % 2 == 1:
            atoms_new   =   Atoms()
            cell    =   [length* b/2 , (width / 2 + .5) * 3 * C_C , 2 * vacuum]
            ymax    =   atoms.get_cell()[1,1]
            for atom in atoms:
                if atom.position[1] < ymax - C_C * 3. / 2:
                    atoms_new += atom
                else:   
                    if atom.number == 1 and cell[1] < atom.position[1] - C_C * 3. / 2: 
                        
                        if cell[0]  <   atom.position[0] + b / 2:
                            atom.position[0]    +=   b / 2 - cell[0]
                        else:
                            atom.position[0]    +=   b / 2
                        
                        atom.position[1]        -=   C_C * 3. / 2
                        atoms_new               +=   atom
                    
                    
            atoms_new.cell  =   cell   
        
            return atoms_new
    
        return atoms

def get_length(ratio, width, edge):
    
    if edge == 'zz':
        length  =   int(ratio * ( 3. * width - 2) / np.sqrt(3) + 1)
        if length % 2 == 0: return length + 1
        else: return length
    
    elif edge == 'ac':
        length  =   int(ratio * np.sqrt(3) * (width - 1) / 6. + 2./3 )
        #length  =   int(ratio * np.sqrt(3) * (width - 1) / 3. + 2./3 )
        
        return length

def get_topInds(atoms):
    
    z_av    =   0.
    top, bot    =   [],[]
    for atom in atoms:
        z_av    +=  atom.position[2]/len(atoms)
    
    for i, atom in enumerate(atoms):
        if z_av < atom.position[2]: top.append(i) 
        elif atom.position[2] < z_av: bot.append(i) 
        else: raise
    return top, bot

def get_rightInds(atoms, topInds):

    xmax    =   0
    for i in topInds:
        if atoms[i].number == 6 and xmax < atoms[i].position[0]: xmax = atoms[i].position[0]
    
    rEdge   =   []
    for i in topInds:
        if atoms[i].number == 6 and xmax - .5 < atoms[i].position[0]:   rEdge.append(i)    
        
    return rEdge
        
def get_posInds(atoms, key = 'redge', *args):
    
    if key == 'redge':
        xmax    =   -10000
        for i in range(len(atoms)):
            if atoms[i].number == 6 and xmax < atoms[i].position[0]: xmax = atoms[i].position[0]
        
        redge       =   []
        for i in range(len(atoms)):
            if atoms[i].number == 6 and xmax - .5 < atoms[i].position[0]:   redge.append(i)    
        
        redge_bot   =   []
        for i in redge:
            if atoms.positions[i,1] <= np.min(atoms.positions[redge, 1]) + .1:
                redge_bot.append(i)

        if len(redge_bot) != 1: raise
        
        redge_top   =   []
        for i in redge:
            if np.max(atoms.positions[redge, 1]) - .1 < atoms.positions[i,1]:
                redge_top.append(i)

        if len(redge_top) != 1: raise


        
        return redge, redge_bot, redge_top
    
    if key == 'ledge':
        xmin    =   10000
        fixL    =   args[0]
        for i in range(len(atoms)):
            if atoms[i].number == 6 and atoms[i].position[0] < xmin: xmin = atoms[i].position[0]
        
        ledge_h   =   []
        for i in range(len(atoms)):
            if atoms[i].position[0] < xmin + fixL :  ledge_h.append(i)    

        ledge_s   =   []
        for i in range(len(atoms)):
            if xmin + fixL < atoms[i].position[0] < xmin + 2 * fixL:  ledge_s.append(i)    

            
        return ledge_s, ledge_h
    

def get_constraints(atoms, edge, bond, idxs, key = 'shear_p', pot = 'LJ'):
    
    if key == 'shear':
        # FIXES
        constraints     =   []
        if edge == 'zz':    fixL    =   np.sqrt(3) * bond * .51 #1.05
        if edge == 'ac':    fixL    =   1.01 * bond
    
        rend_b, rend_t  =   get_posInds(atoms, 'redge')[1:]
        lend_s, lend_h  =   get_posInds(atoms, 'ledge', fixL)
        
        #view(atoms)
        if idxs == None:
            for i in lend_s:
                constraints.append(FixedLine(i, (0,0,1)))
            constraints.append(FixAtoms(indices = lend_h))
        else:
            for i in idxs:
                constraints.append(FixedLine(i, (0,0,1)))
          
        for i in rend_b:
            constraints.append(FixedPlane(i, (0,1,0)))
        
        
        # KC
        add_LJ          =   LJ_potential_smooth(bond)
        constraints.append(add_LJ)
        # END FIXES
        
        return constraints, add_LJ, rend_b, rend_t
    
    elif key == 'twist_F':
        
        # FIXES
        constraints     =   []
        if edge == 'zz':    fixL    =   np.sqrt(3) * bond * 2.05
        if edge == 'ac':    fixL    =   5 * bond
    
        rend_b, rend_t  =   get_posInds(atoms, 'redge')[1:]
        lend_s, lend_h  =   get_posInds(atoms, 'ledge', fixL)
        
        #view(atoms)
        if idxs == None:
            for i in lend_s:
                constraints.append(FixedLine(i, (0,0,1)))
            constraints.append(FixAtoms(indices = lend_h))
        else:
            for i in idxs:
                constraints.append(FixedLine(i, (0,0,1)))
          
        # KC
        add_LJ          =   LJ_potential_smooth(bond)
        if len(rend_b) != len(rend_t) != 1: raise
        twist           =   twist_const_F(rend_b[0], rend_t[0], np.zeros(3))
        
        constraints.append(FixedPlane(rend_b[0], (0,0,1)))
        constraints.append(FixedPlane(rend_t[0], (0,0,1)))
        
        
        constraints.append(add_LJ)
        constraints.append(twist)
        # END FIXES
        
        return constraints, add_LJ, twist, rend_b, rend_t
    
    elif key == 'twist_p':
        
        # FIXES
        constraints     =   []
        if edge == 'zz':    fixL    =   np.sqrt(3) * bond * .51
        if edge == 'ac':    fixL    =   1.1 * bond
    
        rend_b, rend_t  =   get_posInds(atoms, 'redge')[1:]
        lend_s, lend_h  =   get_posInds(atoms, 'ledge', fixL)
        
        #view(atoms)
        if idxs == None:
            for i in lend_s:
                constraints.append(FixedLine(i, (0,0,1)))
            constraints.append(FixAtoms(indices = lend_h))
        else:
            for i in idxs:
                constraints.append(FixedLine(i, (0,0,1)))
          
        # KC
        if pot == 'KC':
            params  =   {}
            params['positions']         =   atoms.positions
            params['chemical_symbols']  =   atoms.get_chemical_symbols()   
            params['ia_dist']           =   10
            params['edge']              =   edge
            params['bond']              =   bond    
            params['ncores']            =   1
            add_pot     =   KC_potential_p(params)
            
        elif pot == 'LJ':
            add_pot =   LJ_potential_smooth(atoms, bond)
        
        
        if len(rend_b) != len(rend_t) != 1: raise
        
        #dist            =   np.linalg.norm(atoms.positions[rend_b[0]] - atoms.positions[rend_t[0]])
        #twist           =   twistConst_Rod(rend_b[0], rend_t[0], dist)
        twist          =   twistConst_Rod(atoms, 4, edge, bond)
        
        constraints.append(FixedPlane(rend_b[0], (0,0,1)))
        
        
        constraints.append(add_pot)
        constraints.append(twist)
        # END FIXES
        
        return constraints, add_pot, twist, rend_b[0], rend_t[0]

def trans_atomsKC(ru, edge, bond):
        
    
    if edge == 'ac':
        nx  =   int(ru[0]/(3*bond))
        ny  =   int(ru[1]/(np.sqrt(3)*bond))
        
        trans   =   np.array([3*bond*nx, np.sqrt(3)*bond*ny, 0]) \
                  - np.array([ru[0] - bond, ru[1], 0.])  
        
    elif edge == 'zz':
    
        nx  =   int(ru[0]/(np.sqrt(3)*bond))
        ny  =   int(ru[1]/(3*bond))
        
        trans   =   np.array([np.sqrt(3)*bond*nx, 3*bond*ny, 0]) \
                -   np.array([ru[0], ru[1] - bond, 0.])
        
    return trans
    

#create_bend_stucture(5, 7, .5, 'ac', 1.39)
    
    