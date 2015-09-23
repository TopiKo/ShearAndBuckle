'''
Created on 17.8.2015

@author: tohekorh
'''
from read import get_datas
from misc.solvers import int_toAngst
import numpy as np
import matplotlib.pyplot as plt
from potentials.KC_imagBottom import make_neighSet, get_setOrig

plt.rc('text', usetex=True)
fig_w   =   6
taito   =   True
T       =   10

def plot_posits(atoms, edge, bond):
    
    positions   =   atoms.positions
    pos_used    =   []
    n_set       =   make_neighSet(15, edge, bond)
            
    for i, r in enumerate(positions):
        if atoms[i].number == 6:
            norms   =   [np.linalg.norm(x) for x in positions[pos_used] - r]
            if len(norms) != 0: 
                if 4 < np.min([np.linalg.norm(x) for x in positions[pos_used] - r]):  
                    neigh_set   =   get_setOrig(r, edge, bond) + n_set
                    pos_used.append(i)
                    plt.scatter(neigh_set[:,0], neigh_set[:,1], color = 'red', alpha = .5)
            else:
                neigh_set   =   get_setOrig(r, edge, bond) + n_set
                pos_used.append(i)
                plt.scatter(neigh_set[:,0], neigh_set[:,1], color = 'red', alpha = .5)
            
    plt.scatter(positions[:,0], positions[:,1], color = 'black')
    plt.axis('equal')
    plt.show()
            
def plot_kinkOfBend(edge):
    
    datas   =   get_datas('LJ', T, '%s_twistTaito' %edge)
    data    =   np.array(datas)
    
    
 
    plt.figure(figsize = (fig_w, 2. / 3 * fig_w)) 
    theta_b =   data[:,0] / (2 * data[:,3])
    theta_v =   data[:,0] / (2 * data[:,6])

    LMax    =   np.max(data[:,1])/10
    plt.scatter(data[:,0], theta_b, color = 'b', label = r'bucle born, $L_{max}=%.2f$nm' %LMax)
    plt.scatter(data[:,0], theta_v, color = 'r', label = r'bucle vani')

    plt.legend(frameon = False, loc = 2)
    plt.text(10, .062, r'Experim. $\theta \approx 0.013$')
    plt.title('Buckling, edge=%s, T=%iK' %(edge, T))
    plt.xlabel(r'width \AA')
    plt.ylabel(r'$\theta$')
    plt.show()


def get_data(data, edge):
    
    W,L =   data[0], data[1]
    Dy  =   data[2]
    R   =   data[3]
    Wf  =   int_toAngst(W, edge, key = 'width', C_C = 1.42)
    Lf  =   int_toAngst(L, edge, key = 'length', C_C = 1.42)
    
    theta   =   Wf / (2 * R)  # new journal of physics Bending and buckling -> theta ~ 0.013 (7-ac ribbon).
    print Lf/Wf, Wf, R, theta, L 
    return [Wf, Lf, Dy, R, theta]
    

#plot_kinkOfBend('ac')