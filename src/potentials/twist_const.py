'''
Created on 27.8.2015

@author: tohekorh
'''
import numpy as np

class twist_const_F:
    
    
    def __init__(self, idxb, idxt, F):
        
        self.idxb   =   idxb      
        self.idxt   =   idxt      
        self.F      =   F
        self.add_set=   'add'
        
    def adjust_positions(self, oldpositions, newpositions):
        pass
    
    
    def adjust_forces(self, posits, forces):
        
        if self.add_set ==   'add':
            forces[self.idxb] +=   self.F  
            forces[self.idxt] +=  -self.F  
        elif self.add_set == 'set':
            forces[self.idxb]  =   self.F  
            forces[self.idxt]  =  -self.F  
     
    def adjust_potential_energy(self, posits, energy):
        
        return 0
    
    def set_F(self, F):
        
        self.F  =   F
        
    
class twistConst_Rod_old:
    
    
    def __init__(self, idxb, idxt, dist):
        
        self.idxb   =   idxb      
        self.idxt   =   idxt
        self.dist   =   dist
        self.F      =   100000.
        self.set_angle(0)
        
    def adjust_positions(self, oldpositions, newpositions):
        
        #step        =   newpositions[self.idxt]  - oldpositions[self.idxt]
        #newpositions[self.idxt] =   oldpositions[self.idxt] + np.dot(self.dir_perp, step) * self.dir_perp 
        
        
        #dist        =   np.linalg.norm(oldpositions[self.idxb] - oldpositions[self.idxt])
        #newpositions[self.idxt] =   newpositions[self.idxb] + self.direc * dist
        pass
        
        
    def adjust_forces(self, posits, forces):
        
        #vec     =   posits[self.idxt] - posits[self.idxb]
        vec_a   =   ( posits[self.idxb] + posits[self.idxt] ) / 2
        vec_z   =   np.array([0., 0., 1.])
        
        M       =   [[self.direc[0], self.dir_perp[0], vec_z[0]],
                     [self.direc[1], self.dir_perp[1], vec_z[1]],
                     [self.direc[2], self.dir_perp[2], vec_z[2]]]
        
        ab      =   np.linalg.solve(M, posits[self.idxb] - vec_a)
        at      =   np.linalg.solve(M, posits[self.idxt] - vec_a)
        
        perp_b  =   ab[1]*self.dir_perp + ab[2]*vec_z
        perp_t  =   at[1]*self.dir_perp + at[2]*vec_z

        r_b     =   np.linalg.norm(perp_b)
        r_t     =   np.linalg.norm(perp_t)
        
        F_b     =   self.F * r_b
        F_t     =   self.F * r_t
        
        #print -F_b * perp_b
        #print -F_t * perp_t
        #print 
        
        forces[self.idxb]  +=   -perp_b * F_b
        forces[self.idxt]  +=   -perp_t * F_t
        
        #pass
        
    def adjust_potential_energy(self, posits, energy):
        
        return 0
    
    def set_angle(self, angle):
        
        self.angle  =   angle
        self.direc  =   np.array([-np.sin(self.angle), \
                                   np.cos(self.angle), \
                                   0.])
        self.dir_perp = np.array([np.cos(self.angle), \
                                  np.sin(self.angle), \
                                  0.])


    def set_dist(self, dist):
        
        self.dist  =   dist

class twistConst_Rod:
    
    
    def __init__(self, atoms, n, edge, bond):
        
        self.pairs  =   get_rightRodInds(atoms, n, edge, bond)
        self.F      =   10.
        self.set_angle(0)
        
    def adjust_positions(self, oldpositions, newpositions):
        
        #step        =   newpositions[self.idxt]  - oldpositions[self.idxt]
        #newpositions[self.idxt] =   oldpositions[self.idxt] + np.dot(self.dir_perp, step) * self.dir_perp 
        
        
        #dist        =   np.linalg.norm(oldpositions[self.idxb] - oldpositions[self.idxt])
        #newpositions[self.idxt] =   newpositions[self.idxb] + self.direc * dist
        pass
        
        
    def adjust_forces(self, posits, forces):
        
        #vec     =   posits[self.idxt] - posits[self.idxb]
        vec_z   =   np.array([0., 0., 1.])
        
        M       =   [[self.direc[0], self.dir_perp[0], vec_z[0]],
                     [self.direc[1], self.dir_perp[1], vec_z[1]],
                     [self.direc[2], self.dir_perp[2], vec_z[2]]]
        
        
        for pair in self.pairs:
            
            idxb, idxt  =   pair
            vec_a   =   ( posits[idxb] + posits[idxt] ) / 2

            ab      =   np.linalg.solve(M, posits[idxb] - vec_a)
            at      =   np.linalg.solve(M, posits[idxt] - vec_a)
            
            perp_b  =   ab[1]*self.dir_perp + ab[2]*vec_z
            perp_t  =   at[1]*self.dir_perp + at[2]*vec_z
        
            
            F_b     =   - self.F * perp_b
            F_t     =   - self.F * perp_t
            
            
            forces[idxb]  +=   F_b
            forces[idxt]  +=   F_t
        
        #pass
        
    def adjust_potential_energy(self, posits, energy):
        
        return 0
    
    def set_angle(self, angle):
        
        self.angle  =   angle
        self.direc  =   np.array([-np.sin(self.angle), \
                                   np.cos(self.angle), \
                                   0.])
        self.dir_perp = np.array([np.cos(self.angle), \
                                  np.sin(self.angle), \
                                  0.])


    def set_dist(self, dist):
        
        self.dist  =   dist

def get_rightRodInds(atoms, n, edge, bond):
    
    atoms_list = []
    for i, atom in enumerate(atoms):
        if atom.number == 6:
            atoms_list.append([i, atom.position[0], atom.position[1], atom.position[2]])
            
    atoms_list  =   np.array(atoms_list)
    pairs       =   np.zeros((n,2))
    
    if edge == 'ac':
        period  =   3. / 4. * bond
    elif edge == 'zz':
        period  =   np.sqrt(3.) / 2. * bond
        
    xmax    =   np.max(atoms_list[:,1])
    intvals =   [xmax + period/2 - period*i for i in range(n + 1)]
    for i in range(n):
        xmin, xmax  =   intvals[i + 1], intvals[i]
        idxs        =   np.where(xmin < atoms_list[:,1])[0]
        idxs2       =   np.where(atoms_list[:,1] < xmax)[0]
        idxs        =   set(idxs).intersection(idxs2)
        ymax, ymin  =   -10, 1000
        
        for idx in idxs:
            if ymax < atoms_list[idx,2]: 
                ymax    =   atoms_list[idx,2]    
                ymax_i  =   idx   
            if atoms_list[idx,2] < ymin: 
                ymin    =   atoms_list[idx,2]    
                ymin_i  =   idx   
        
        pairs[i]    =   [atoms_list[ymin_i, 0], atoms_list[ymax_i, 0]]
        
    return pairs
    
        