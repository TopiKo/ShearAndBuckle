'''
Created on 13.8.2015

@author: tohekorh
'''

from ase import Atoms, Atom
from ase.visualize import view
from ase.constraints import FixedPlane, FixedLine, FixAtoms
from potentials.LJ_smooth import LJ_potential_smooth
import numpy as np

def create_stucture(ratio, width, edge, key = 'top', a = 1.42):
    
    vacuum  =   10
    length  =   get_length(ratio, width, edge)
    length_b=   int(length * .8)*2
    width_b =   width * 4
    
    
    print length, width
    
    if edge == 'zz':    
        orig    =   [np.sqrt(3)*a, 2*a] # -> AB-stacking
    if edge == 'ac':    
        orig    =   [2*a, np.sqrt(3)*a] # -> AB-stacking
    
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
                                    C_C=1.42, 
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
        atoms.cell  =   [top.cell[0,0]*1.5, top.cell[0,0] + 2*top.cell[1,1], 2 * vacuum]
        atoms.center()
        atoms.positions[:,2] = 3.4
        atoms.pbc   =   [False, False, False]
        return atoms, L, W, length


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
        length  =   int(ratio * np.sqrt(3) * (width - 1) / 3. + 2./3 )
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
    

def get_constraints(atoms, edge, bond):
    
    # FIXES
    constraints     =   []
    if edge == 'zz':    fixL    =   np.sqrt(3) * bond * 2.05
    if edge == 'ac':    fixL    =   5 * bond

    rend_b, rend_t  =   get_posInds(atoms, 'redge')[1:]
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
    
    return constraints, add_LJ, rend_b, rend_t
    
    
    