'''
Created on 17.9.2015

@author: tohekorh
'''
import numpy as np
from structure_help import nrst_neigh, local_normal
from numpy import zeros, exp, dot, sqrt
import scipy
import time
import multiprocessing
from ase.structure import graphene_nanoribbon

# Taidolla taytyy lisata ->
# os.system("taskset -p 0xff %d" % os.getpid())
# ks. http://bugs.python.org/issue17038#msg180663


def get_setOrig(r, edge, bond):
    
    if edge == 'ac':
        nx  =   int(r[0]/(3*bond))
        ny  =   int(r[1]/(np.sqrt(3)*bond))
        
        set_orig    =   [3*bond*nx, np.sqrt(3)*bond*ny, 0]
        
        return set_orig
        
    elif edge == 'zz':
        
        nx  =   int(r[0]/(np.sqrt(3)*bond))
        ny  =   int(r[1]/(3*bond))
        
        set_orig    =   [np.sqrt(3)*bond*nx, 3*bond*ny, 0]
        
        return set_orig
        
        

def gradN(ri, ind, positions, layer_neighbors, take_t = False):
    # Matrix
    #dn1/dx dn2/dx dn3/dx
    #dn1/dy dn2/dy dn3/dy
    #dn1/dz dn2/dz dn3/dz
    
    #return np.ones((3,3))
    
    def nj(xi, *args):
        i                       =   args[0]
        
        posits_use              =   positions.copy()
        posits_use[ind, i]      =   xi
        #posits_ext_use          =   extend_structure(posits_use.copy(), pbc, cell)
               
        nj                      =   local_normal(ind, posits_use, layer_neighbors)
                
        return nj
        
        
    dni = zeros((3,3))
    
    for i in range(3):
        dni[i,:]    =   scipy.misc.derivative(nj, ri[i], dx=1e-3, n=1, args=[i], order = 3)
    
    return dni
    


def GradV_ij(rij, r, pij, pji, ni, nj, dni, dnj, params):
    
    # Here we calculate the gradient of the KC- potential to obtain force:
    # See notes, for derivation.    
    A, z0, lamna, delta, C, C0, C2, C4  =   params
    
    gpij        =   -2./delta**2*f(pij, delta, C0, C2, C4) + exp(-(pij/delta)**2)* \
                    (2.*C2/delta**2 + 4*C4*pij**2/delta**4)
    gpji        =   -2./delta**2*f(pji, delta, C0, C2, C4) + exp(-(pji/delta)**2)* \
                    (2.*C2/delta**2 + 4*C4*pji**2/delta**4)
    

    njdotrij    =   nj[0]*rij[0] + nj[1]*rij[1] + nj[2]*rij[2]  # nj.rij
    nidotrij    =   ni[0]*rij[0] + ni[1]*rij[1] + ni[2]*rij[2]  # ni.rij
    
    bulkPreFac  =  -lamna*exp(-lamna*(r - z0))*(C + f(pij, delta, C0, C2, C4) + f(pji, delta, C0, C2, C4)) \
                +   6*A*z0**6/r**7
    expPreFac   =   exp(-lamna*(r - z0))

    
    rijDotdni   =   dot(dni, rij)
    
    FIJ         =  -(bulkPreFac*(-rij/r) \
                   + expPreFac*gpji*(njdotrij*nj - rij) \
                   + expPreFac*gpij*(nidotrij*(ni  - rijDotdni) - rij))
    
    
    ###################
    # THIS TERMS IS CORRECTION BECAUSE IN KC FIJ != FJI, see deriv from notes...
    rijDotdnj   =   dot(dnj, rij)

    fji_add     =   exp(-lamna*(r - z0))*(gpij*(nidotrij)*rijDotdni + gpji*(njdotrij)*rijDotdnj)
    #Fij + Fji  =    fij_add
    
    ###################
    
    return FIJ, fji_add - FIJ

def f(p, delta, C0, C2, C4):
    #A, z0, lamna, delta, C, C0, C2, C4 =   params
    
    return np.exp(-(p/delta)**2)*(C0 \
                 +  C2*(p/delta)**2  \
                 +  C4*(p/delta)**4)


def get_potential_ij(ri, rj, ni, nj, positions, i, params, \
                     cutoff):
    
    A, z0, lamna, delta, C, C0, C2, C4 =   params
    
    e_KC            =   0.

    
    rij =   rj - ri
    r   =   sqrt(rij[0]**2 + rij[1]**2 + rij[2]**2)
    
    if r < cutoff:
        
        niDotrij    =   ni[0]*rij[0] + ni[1]*rij[1] + ni[2]*rij[2] 
        njDotrij    =   nj[0]*rij[0] + nj[1]*rij[1] + nj[2]*rij[2] 
        
        if -1e-12 < r**2 - niDotrij**2 < 0.:
            pij     =   0.
        else:
            pij     =   sqrt(r**2 - niDotrij**2)
            
        if -1e-12 < r**2 - njDotrij**2 < 0.:
            pji     =   0.
        else:
            pji     =   sqrt(r**2 - njDotrij**2)

        
        e_KC       +=   exp(-lamna*(r - z0))*(C + f(pij, delta, C0, C2, C4) \
                                                + f(pji, delta, C0, C2, C4)) \
                                                - A*(r/z0)**(-6.)
        
                

    return e_KC

def grad_pot_ij(posits, i, rj, params, \
                cutoff, layer_neighbors):                                
    
    def e(rik, *args):
        k               =   args[0]
        posits_use      =   posits.copy()   
        posits_use[i,k] =   rik
        ri              =   posits_use[i]
        #posits_ext      =   extend_structure(posits_use, pbc, cell)

        ni          =   local_normal(i, posits, layer_neighbors)*(-1)
        nj          =   np.array([0.,0.,1.])#local_normal(j, posits_ext, layer_neighbors)*njfac
        
        
        return get_potential_ij(ri, rj, ni, nj, posits_use, i, params, cutoff)
             
        
    de          =   zeros(3)
    
    for k in range(3):
        de[k]   =   scipy.misc.derivative(e, posits[i,k], dx=1e-6, n=1, \
                                          args=[k], order = 3)    
        
    return de

def get_forces_ij(ri, rj, ni, nj, dni, dnj, positions, i, params, \
                      cutoff, layer_neighbors = None):
        
    Fij             =   zeros(3)
    Fji             =   zeros(3)

    rij         =   rj - ri
    r           =   sqrt(rij[0]**2 + rij[1]**2 + rij[2]**2) 
    
    if r < cutoff:
        niDotrij    =   ni[0]*rij[0] + ni[1]*rij[1] + ni[2]*rij[2] 
        njDotrij    =   nj[0]*rij[0] + nj[1]*rij[1] + nj[2]*rij[2] 
        
        if -1e-12 < r**2 - niDotrij**2 < 0.:
            pij     =   0.
        else:
            pij     =   sqrt(r**2 - niDotrij**2)
            
        if -1e-12 < r**2 - njDotrij**2 < 0.:
            pji     =   0.
        else:
            pji     =   sqrt(r**2 - njDotrij**2)
        
        #F          += - GradV_ij(rij, r, pij, pji, ni, nj, dni, params)
        
        F1, F2      =   GradV_ij(rij, r, pij, pji, ni, nj, dni, dnj, params)
        Fij        +=   F1
        Fji        +=   F2 
                
    
    
    # test that the force is -grad Vij. Very expensive!
    if layer_neighbors != None:
        
        fij = -grad_pot_ij(positions.copy(), i, rj, params, \
                              cutoff, layer_neighbors)
        #fji= -grad_pot_ij(positions.copy(), rj, i, params, \
        #                      cutoff, layer_neighbors, not jtop)

        print 'force to %i' %(i)
        for ll in range(3):
            print Fij[ll] - fij[ll], Fij[ll]
            if 1e-6 < np.abs(fij[ll]) or 1e-6 < np.abs(Fij[ll]): 
                if np.abs(Fij[ll] - fij[ll])/np.linalg.norm(fij) > 1e-1: #
                #if np.abs(F[ll] - f[ll]) > 1e-6:
                
                    print np.abs(Fij[ll] - fij[ll])/np.linalg.norm(fij)
                    print ri - rj, np.linalg.norm(ri - rj)
                    print (Fij - fij)/np.linalg.norm(fij)
                    print Fij
                    print fij
                    raise
        
        #print 'forces from all images of %i to %i' %(i,rj)
        #for ll in range(3):
        #    print Fji[ll] - fji[ll], Fji[ll]
        #    if 1e-6 < np.abs(fji[ll]) or 1e-6 < np.abs(Fji[ll]): 
        #        if np.abs(Fji[ll] - fji[ll]) > 1e-6:
        #        #if np.abs(F[ll] - f[ll]) > 1e-6:
        #        
        #            print i,rj, jtop
        #            print Fij
        #            print f
        #            raise
        #print 
        
    
    return Fij


def get_neigh_layer_indices(layer, layer_indices):
    
    neigh_ind  =   []
    
    if layer > 0:
        neigh_ind.append(['bottom', layer_indices[layer - 1]])
        
    if layer < len(layer_indices) - 1:
        neigh_ind.append(['top', layer_indices[layer + 1]])
        
        
    return neigh_ind   


def make_neighSet(cutoff, edge, bond):
   
    atoms   =   graphene_nanoribbon(50, 50, type='zigzag', C_C = bond)
    atoms.rotate([1,0,0], np.pi/2, rotate_cell = True)
    
    if edge == 'zz':
        atoms.rotate([0,0,1], np.pi/2, rotate_cell = True)
        if cutoff != 15:
            print 'not implemented'
            raise
    if edge == 'ac':
        if cutoff != 15:
            print 'not implemented... cannot be sure of the alingmen for other cutoffs.'
            raise
    
        
    atoms.center()
        
    middle  =   [np.average(atoms.positions[:,0]), \
                 np.average(atoms.positions[:,1]), \
                 np.average(atoms.positions[:,2])]
    
    rem     =   []
    midMin  =   100
    for i, atom in enumerate(atoms):
        if cutoff < np.linalg.norm(atom.position - middle): rem.append(i)
             
    del atoms[rem]
    
    for i, pos in enumerate(atoms.positions):
        if np.linalg.norm(pos - middle) < midMin:
            midMin  =   np.linalg.norm(pos - middle)
            iCentre =   i

    return atoms.positions - atoms.positions[iCentre]
   
class KC_potential_p:
    
    # This is a constraint class to include the registry dependent 
    # interlayer potential to ase. 
    
    def __init__(self, params):
        
        posits          =   params['positions']
        self.chem_symbs =   params['chemical_symbols']
        self.cutoff     =   params['ia_dist']
        self.edge       =   params['edge']
        self.bond       =   params['bond']
        self.neigh_set  =   make_neighSet(15, self.edge, self.bond)
        
        # parameters as reported in original KC paper
        self.delta  =   0.578       # Angstrom   
        self.C0     =   15.71*1e-3  # [meV]
        self.C2     =   12.29*1e-3  # [meV]
        self.C4     =   4.933*1e-3  # [meV]
        self.C      =   3.030*1e-3  # [meV]
        self.lamna  =   3.629       # 1/Angstrom
        self.z0     =   3.34        # Angstrom
        self.A      =   10.238*1e-3 # [meV]     
        
        self.params =   [self.A, self.z0, self.lamna, self.delta, self.C, self.C0, self.C2, self.C4]
        
        self.cores              =   max((2, params['ncores']))
        self.layer_neighbors    =   nrst_neigh(posits, 'layer')    
        
         
    def adjust_positions(self, oldpositions, newpositions):
        pass
    
    
    def adjust_forces(self, posits, forces):
        
        params          =   self.params
        cutoff          =   self.cutoff
        
        layer_neighbors =   self.layer_neighbors
        chem_symbs      =   self.chem_symbs 
        
        dnGreat         =   np.empty(len(forces), dtype = 'object')
        
        
        # PARALLELL BEGIN
        
        # Divide the loop over all atoms into split_len parts.
        
        if self.cores > 1:
            split_len   =   float(len(posits))/self.cores
            split_arr   =   [int((i + 1)*split_len) for i in range(self.cores - 1)]
        else:
            split_len   =   len(posits)
            split_arr   =   [0]
        split_posits    =   np.split(posits, split_arr)
        split_arr       =   np.insert(split_arr, 0, 0)
        
        
        # Arrays for paralell computing
        jobs        =   []
        ar11        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar12        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar13        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar21        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar22        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar23        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar31        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar32        =   multiprocessing.Array('d', np.zeros(len(forces)))
        ar33        =   multiprocessing.Array('d', np.zeros(len(forces)))
            
        
        for i, pos_arr in enumerate(split_posits):
            
            args        =   (i, split_arr[i], pos_arr, \
                            ar11, ar12, ar13,  \
                            ar21, ar22, ar23,  \
                            ar31, ar32, ar33,  \
                            posits, layer_neighbors, chem_symbs)
            
            process     =   multiprocessing.Process(target=self.pre_gradN, 
                                                    args=args)
            jobs.append(process)

        # Start the processes (i.e. calculate the random number lists)        
        for j in jobs:
            j.start()
        
        # Ensure all of the processes have finished
        for j in jobs:
            j.join()
        

        for i in range(len(forces)):
            dnGreat[i]      =   np.zeros((3,3))
            #   SHEAR ADD 
            dnGreat[i][0,0] =   ar11[i]
            dnGreat[i][0,1] =   ar12[i]
            dnGreat[i][0,2] =   ar13[i]
            dnGreat[i][1,0] =   ar21[i]
            dnGreat[i][1,1] =   ar22[i]
            dnGreat[i][1,2] =   ar23[i]
            dnGreat[i][2,0] =   ar31[i]
            dnGreat[i][2,1] =   ar32[i]
            dnGreat[i][2,2] =   ar33[i]
    

        
        jobs        =   []
        af1         =   multiprocessing.Array('d', np.zeros(len(forces)))
        af2         =   multiprocessing.Array('d', np.zeros(len(forces)))
        af3         =   multiprocessing.Array('d', np.zeros(len(forces)))
        
        
        for k, pos_arr in enumerate(split_posits):
            

            args        =   (split_arr[k], pos_arr,  af1, af2, af3, dnGreat,\
                             posits, layer_neighbors, \
                             chem_symbs, params, cutoff, len(forces))
            
            process     =   multiprocessing.Process(target=self.forces_on_posArr, 
                                                    args=args)
            jobs.append(process)

        # Start the processes         
        for j in jobs:
            j.start()
        
        # Ensure all of the processes have finished
        for j in jobs:
            j.join()
        
        
        for i in range(len(forces)):
            forces[i,0]     +=  af1[i]
            forces[i,1]     +=  af2[i]
            forces[i,2]     +=  af3[i]
            
    def pre_gradN(self, i, split_beg, pos_arr,  ar11, ar12, ar13, \
                                                ar21, ar22, ar23, \
                                                ar31, ar32, ar33, \
                                                posits, \
                                                layer_neighbors, chem_symbs):
        
        for j, r in enumerate(pos_arr):
            ind_atom    =   j  + split_beg
            if chem_symbs[ind_atom] == 'C': # SHEAR ADD  (and ind_atom in..)
                dni             =   gradN(r, ind_atom, posits, layer_neighbors)
                ar11[ind_atom]  =   dni[0,0]
                ar12[ind_atom]  =   dni[0,1]
                ar13[ind_atom]  =   dni[0,2]
                ar21[ind_atom]  =   dni[1,0]
                ar22[ind_atom]  =   dni[1,1]
                ar23[ind_atom]  =   dni[1,2]
                ar31[ind_atom]  =   dni[2,0]
                ar32[ind_atom]  =   dni[2,1]
                ar33[ind_atom]  =   dni[2,2]
    
      
        
    def forces_on_posArr(self, split_beg, pos_arr,  af1, af2, af3, dnGreat, \
                             posits, layer_neighbors, \
                             chem_symbs, params, cutoff, natoms):
        
        
        nj          =   np.array([0., 0., 1.])    
        dnj         =   np.array([[0., 0., 0],
                                  [0., 0., 0],
                                  [0., 0., 0]])

        for l, ri in enumerate(pos_arr):
            i   =   l  + split_beg 
            if chem_symbs[i] == 'C': # SHREAR ADD
                
                ni_f        =   local_normal(i, posits, layer_neighbors)
                dni_f       =   dnGreat[i]
                    
                ni          =   -1*ni_f
                dni         =   -dni_f
                

                neigh_orig  =   get_setOrig(ri, self.edge, self.bond)
                neigh_pos   =   self.neigh_set + neigh_orig
                
                
                
                for rj in neigh_pos: 
                    
                    # Force due to atom j on atom i
                    fij     =   get_forces_ij(ri, rj, ni, nj, dni, dnj, \
                                              posits, i, params, \
                                              cutoff) #, layer_neighbors = layer_neighbors)
            
                    
                    af1[i]  +=  fij[0]
                    af2[i]  +=  fij[1]
                    af3[i]  +=  fij[2]
                        
     
    def adjust_potential_energy(self, posits, energy):
        
        params          =   self.params
        cutoff          =   self.cutoff
        layer_neighbors =   self.layer_neighbors
        chem_symbs      =   self.chem_symbs
        
        split_len       =   float(len(posits))/self.cores
        split_arr       =   [int((i + 1)*split_len) for i in range(self.cores - 1)]
        split_posits    =   np.split(posits, split_arr)
        split_arr       =   np.insert(split_arr, 0, 0)
        
        
        jobs            =   []
        e_KCm           =   multiprocessing.Array('d', np.zeros(self.cores))
        
        
        for k, pos_arr in enumerate(split_posits):
            
            args        =   (k, split_arr[k], pos_arr,  e_KCm, \
                             posits, layer_neighbors, \
                             chem_symbs, params, cutoff, len(posits))
            
            process     =   multiprocessing.Process(target=self.energy_on_posArr, 
                                                    args=args)
            jobs.append(process)

        # Start the processes        
        for j in jobs:
            j.start()
        
        # Ensure all of the processes have finished
        for j in jobs:
            j.join()
        
        
        e_KC    =   0.
        for k in range(self.cores):
            e_KC    +=  e_KCm[k]
        
        return e_KC
        
    
    def energy_on_posArr(self, k, split_len, pos_arr, e_KCm, \
                             posits, layer_neighbors, 
                             chem_symbs, params, cutoff, natoms):
        
        nj          =   np.array([0., 0., 1.])    
                
        for l, ri in enumerate(pos_arr):
            i   =   l  + split_len 
            if chem_symbs[i] == 'C': # SHEAR ADD
                ni          =   local_normal(i, posits, layer_neighbors)
                ni          =   -1*ni
                    
                neigh_orig  =   get_setOrig(ri, self.edge, self.bond)
                neigh_pos   =   self.neigh_set + neigh_orig
            
                for rj in neigh_pos: 
                       
                    e       =   get_potential_ij(ri, rj, ni, nj, posits, i, params, \
                                                 cutoff) 
            
                    e_KCm[k]       +=  e

    
    # TETSTEST
    def energy_i(self, posits):
        
        params          =   self.params
        cutoff          =   self.cutoff
        #posits_ext      =   extend_structure(posits.copy(), self.pbc, self.cell)
        #layer_indices   =   self.layer_indices
        layer_neighbors =   self.layer_neighbors
        chem_symbs      =   self.chem_symbs
        
        e_KC            =   np.zeros(len(posits))
        
        for i, ri in enumerate(posits):
            if chem_symbs[i] == 'C':
                ni              =   local_normal(i, posits, layer_neighbors)
                
                ni          =   -1*ni
                nj          =   np.array([0., 0., 1.])    
                

                                         
                neigh_orig  =   get_setOrig(ri, self.edge, self.bond)
                neigh_pos   =   self.neigh_set + neigh_orig
            
                for rj in neigh_pos: 
                        
                    e, new_maps     =   get_potential_ij(ri, rj, ni, nj, posits, i, params, \
                                                 cutoff) 
            
                    if new_maps != None:    self.map_seqs   =   new_maps
                    
                    e_KC[i]       +=  e

        return e_KC
        
        
    # TETSTEST



    