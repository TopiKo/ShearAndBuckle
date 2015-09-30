'''
Created on 17.8.2015

@author: tohekorh
'''
from read import get_datas, get_log_data
from misc.solvers import int_toAngst
import numpy as np
import matplotlib.pyplot as plt
from potentials.KC_imagBottom import make_neighSet, get_setOrig

plt.rc('text', usetex=True)
fig_w   =   6
taito   =   True
T       =   10

def plot_posits(atoms, edge, bond, vecs=None):
    
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
    
    if vecs != None:
        for vec in vecs:
            plot_vec    =   np.array([positions[0],
                                      positions[0] + vec])  
            plt.plot(plot_vec[:,0], plot_vec[:,1])
    plt.axis('equal')
    plt.show()
            
def plot_kinkOfBend(edge):
    
    datas   =   get_datas('LJ', T, '%s_twistTaito' %edge)
    data    =   np.array(datas)
    
    plt.figure(figsize = (fig_w, 2. / 3 * fig_w)) 
    theta_b =   data[:,0] / (2 * data[:,3])
    theta_v =   data[:,0] / (2 * data[:,4])
    v       =   data[:,2]/np.max(data[:,2])
    
    LMax    =   np.max(data[:,1])/10
    
    
    set_label   =   True
    for i in range(len(v)):
        if v[i] == 1 and set_label:
            set_label   =   False
            plt.scatter(data[i,0], theta_b[i], color = 'b', alpha = v[i], 
                        label = r'bucle born, $L_{max}=%.2f$nm' %LMax)
            plt.scatter(data[i,0], theta_v[i], color = 'r', alpha = v[i], 
                        label = r'bucle vanish')
        else:
            plt.scatter(data[i,0], theta_b[i], color = 'b', alpha = v[i])
            plt.scatter(data[i,0], theta_v[i], color = 'r', alpha = v[i])
    plt.plot(data[:,0], np.ones(len(data[:,0]))*.04, '-.', color = 'black')
    plt.text(np.min(data[:,0]), .041, 'teor')
    plt.plot(data[:,0], np.ones(len(data[:,0]))*.013, '-.', color = 'black')
    plt.text(np.min(data[:,0]), .014, r'Experim. $\theta \approx 0.013$')
    
    
    plt.legend(frameon = False, loc = 1)
    #plt.text(10, .062, r'Experim. $\theta \approx 0.013$')
    plt.title('Buckling, edge=%s, T=%iK' %(edge, T))
    plt.xlabel(r'width \AA')
    plt.ylabel(r'$\theta$')
    plt.show()


def plot_energy(edge):
    
    datas   =   get_log_data('LJ', T, '%s_twistTaito' %edge)
    k       =   21
    
    def shear_e(W, L, R):
        
        return k*L*(W/R)**2*W/12
        #return k*L*theta**2*W/3.
        
        
    for data in datas:
        energy_table    =   data[2]
        W, L            =   data[:2]
        R               =   energy_table[:,2]
        z               =   energy_table[:,4]
        epot            =   energy_table[:,5] - energy_table[0,5]
        
        thetas          =   W/(2*energy_table[:,2])
        epot_t          =   shear_e(W, L, R)
        
        plt.plot(thetas, epot, label = 'Epot')
        plt.plot(thetas, epot_t, label = 'Epot teor')
        plt.legend(loc = 1)
        plt.ylabel('tot E eV')
        plt.twinx()
        plt.plot(thetas, z, color = 'red', label = 'maxH')
        plt.legend(loc = 2)
        plt.ylabel('height max Angst')
        plt.xlabel('Theta w/(2R)')
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
    

#plot_energy('ac')
#plot_kinkOfBend('ac')