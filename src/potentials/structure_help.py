'''
Created on 28.10.2014

@author: tohekorh
'''
import numpy as np
from numpy import sqrt

def extend_structure(posits, pbc, cell):
    
    '''
    posits_ext = posits.copy()
    
    for i in range(3):
        if pbc[i]:
            
            #n = int(1 + (ia_length + 5)/cell[i])
            n = 1
            add_posits = np.zeros((2*n*len(posits_ext), 3))
            m = 0
            kset = np.concatenate((range(-n,0), range(1,n + 1)))
            for k in kset:
                for pos in posits_ext:
                    add_posits[m] = pos
                    add_posits[m][i] += k*cell[i]
                    m += 1
            posits_ext = np.concatenate((posits_ext, add_posits))
    
    return posits_ext 
    '''
    
    # Extend tthe structure to the periodic directions by one unit cell.
    posits_ext  =   posits.copy()
    
    for i in range(3):
        if pbc[i]:
            pos_plus        =   posits_ext.copy()
            pos_minus       =   posits_ext.copy()
            pos_plus[:,i]  +=   cell[i]
            pos_minus[:,i] -=   cell[i]
            
            posits_ext  = np.concatenate((posits_ext, pos_plus, pos_minus))
    
    return posits_ext
            
    

def which_layer(i, layer_inds):
    
    for j, inds in enumerate(layer_inds):
        if i in inds: return j, layer_inds[j]
    else: raise
    
def nrst_neigh(posits, key, *args):
    
    if key == 'layer':
        d               =   1.5
        n               =   len(posits)
        layer_inds_ext  =   find_layers(posits)[1]
        neighbours      =   np.empty(n, dtype = object)
        for i in range(n):
            
            inds            =   which_layer(i, layer_inds_ext)[1]
            r               =   posits[i]
            neighbours[i]   =   []
            for ind in inds:
                if ind != i:
                    if np.linalg.norm(r - posits[ind]) < d:
                        neighbours[i].append(ind)
        
        return neighbours

   
    
def map_seq(m, n):

    
    perms       =   []    
    if n == 1:
        for i in range(-m, m +1):   #, n
            if np.abs(i) == m:
                perms.append([i])    
            
    
    elif n == 2:
        for i in range(-m, m +1):   #, n
            for j in range(-m, m +1):
                if np.abs(i) + np.abs(j) == m:
                    perms.append([i,j])    
        
    elif n == 3:
        for i in range(-m, m +1):   #, n
            for j in range(-m, m +1):
                for k in range(-m, m +1):
                    if np.abs(i) + np.abs(j) + np.abs(k) == m:
                        perms.append([i,j, k])    
    return perms
                   
def map_rj(rj, map_seq, pbc, cell):      
    
    rjs     =   np.zeros((len(map_seq), 3))
    
    for k, map_s in enumerate(map_seq):
        r   =   rj.copy()
        l   =   0
        for i in range(3):
            if pbc[i]:
                r[i]    +=  map_s[l]*cell[i]  
                l       +=  1
        rjs[k]   =   r
    
    return rjs

def map_rjs(rj, pbc, cell, n, ia_length):      
    
    
    perms       =   []    
    
    if n == 1:
        
        for i in range(3):
            if pbc[i]:
                cell_L  =   cell[i]
        
        m   =   int(ia_length/cell_L) + 1
        for i in range(-m, m +1):
            perms.append([i])    
    
    elif n == 2:
        k       = 0
        cell_L  =   np.zeros(2)
        for i in range(3):
            if pbc[i]:
                cell_L[k]  =   cell[i]
                k         += 1

        m1      =   int(ia_length/cell_L[0]) + 1
        m2      =   int(ia_length/cell_L[0]) + 1

        for i in range(-m1, m1 +1):
            for j in range(-m2, m2 +1):
                perms.append([i,j])    
    
    
    rjs     =   np.zeros((len(perms), 3))
    
    for k, map_s in enumerate(perms):
        r   =   rj.copy()
        l   =   0
        for i in range(3):
            if pbc[i]:
                r[i]    +=  map_s[l]*cell[i]  
                l       +=  1
        rjs[k]   =   r
    
    return rjs


def local_normal(i, posits, layer_neighbors, const = False):
    
    
    #tang_vec    =   np.zeros((len(layer_neighbors[i]), 3))
    #for k, j in enumerate(layer_neighbors[i]):
    #    tang_vec[k]  =  posits_ext[j] - ri 
    
    
    if not const:
        ri          =   posits[i]
        tang_vec    =   posits[layer_neighbors[i]] - ri
        
        if len(tang_vec) == 3:
            '''
            normal  =   np.cross(tang_vec[0], tang_vec[1])/np.linalg.norm(np.cross(tang_vec[0], tang_vec[1])) \
                    +   np.cross(tang_vec[2], tang_vec[0])/np.linalg.norm(np.cross(tang_vec[2], tang_vec[0])) \
                    +   np.cross(tang_vec[1], tang_vec[2])/np.linalg.norm(np.cross(tang_vec[1], tang_vec[2])) 
    
            normal  =   normal/np.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
            '''
            
            n1      =   np.zeros(3)  # = tang_vec[0] x tang_vec[1]
            n1[0]   =   tang_vec[0,1]*tang_vec[1,2] - tang_vec[0,2]*tang_vec[1,1]
            n1[1]   = -(tang_vec[0,0]*tang_vec[1,2] - tang_vec[0,2]*tang_vec[1,0])
            n1[2]   =   tang_vec[0,0]*tang_vec[1,1] - tang_vec[0,1]*tang_vec[1,0]
    
            n1      =   n1/sqrt(n1[0]**2 + n1[1]**2 + n1[2]**2)
            
            if n1[2] < 0:
                n1  =   -n1
            
            n2      =   np.zeros(3) # = tang_vec[2] x tang_vec[0]
            n2[0]   =   tang_vec[2,1]*tang_vec[0,2] - tang_vec[2,2]*tang_vec[0,1]
            n2[1]   = -(tang_vec[2,0]*tang_vec[0,2] - tang_vec[2,2]*tang_vec[0,0])
            n2[2]   =   tang_vec[2,0]*tang_vec[0,1] - tang_vec[2,1]*tang_vec[0,0]
    
            n2      =   n2/sqrt(n2[0]**2 + n2[1]**2 + n2[2]**2)
            
            if n2[2] < 0:
                n2  =   -n2
    
    
            n3      =   np.zeros(3) # = tang_vec[1] x tang_vec[2]
            n3[0]   =   tang_vec[1,1]*tang_vec[2,2] - tang_vec[1,2]*tang_vec[2,1]
            n3[1]   = -(tang_vec[1,0]*tang_vec[2,2] - tang_vec[1,2]*tang_vec[2,0])
            n3[2]   =   tang_vec[1,0]*tang_vec[2,1] - tang_vec[1,1]*tang_vec[2,0]
    
            n3      =   n3/sqrt(n3[0]**2 + n3[1]**2 + n3[2]**2)
    
            if n3[2] < 0:
                n3  =   -n3
    
            
            
            normal  =   n1 + n2 + n3
            normal  =   normal/sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
                    
            return normal
            
        
        elif  len(tang_vec) == 2:
            #normal  =   np.cross(tang_vec[0], tang_vec[1])/np.linalg.norm(np.cross(tang_vec[0], tang_vec[1]))
    
            n1      =   np.zeros(3)  # = tang_vec[0] x tang_vec[1]
            n1[0]   =   tang_vec[0,1]*tang_vec[1,2] - tang_vec[0,2]*tang_vec[1,1]
            n1[1]   = -(tang_vec[0,0]*tang_vec[1,2] - tang_vec[0,2]*tang_vec[1,0])
            n1[2]   =   tang_vec[0,0]*tang_vec[1,1] - tang_vec[0,1]*tang_vec[1,0]
    
            normal  =   n1/sqrt(n1[0]**2 + n1[1]**2 + n1[2]**2)
    
            if normal[2] < 0:
                normal  =   -normal
            
            return normal
            
        else:
            print i, layer_neighbors[i]
            print ri, posits[layer_neighbors[i]]
            raise
    
    else:
        
        normal      =   np.zeros(3)
        normal[2]   =   1.
        return normal 

def find_layers(positions):
    
    zs = positions[:,2]
    layers = [[1, zs[0]]]
    for z in zs:
        add = True
        for layer in layers:
            if abs(z - layer[1]) < .5:
                add = False
                l_idx   = layers.index(layer)
                n       = layer[0]
                layers[l_idx][1] = (n*layer[1] + z)/(n + 1)
                layers[l_idx][0] = n + 1
                
        if add: layers.append([1, z])
    
    layers_new = []
    for layer in layers:
        layers_new.append(layer[1])
    
    layers = layers_new    
    layers.sort()
    inds_L  = []
    for layer in layers:
        inds = []
        for ir, r in enumerate(positions):    
            z   =   r[2]
            if abs(z - layer) < 0.1:
                inds.append(ir)
        inds_L.append(inds)
                
    return layers, inds_L
        