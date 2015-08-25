'''
Created on 17.8.2015

@author: tohekorh
'''
from read import get_datas
from misc.solvers import int_toAngst
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
fig_w   =   6
taito   =   True

def plot_kinkOfBend(edge):
    
    #datas1      =   get_datas('LJ', edge, taito)
    datas2      =   get_datas('LJ', edge + '_ff', taito)

    #full_data1  =   np.zeros((len(datas1),5))
    full_data2  =   np.zeros((len(datas2),5))
    
    #for i, data in enumerate(datas1):
    #    full_data1[i]    =  get_data(data, edge) 
          
    for i, data in enumerate(datas2):
        full_data2[i]    =  get_data(data, edge)   
    
    #print full_data1[:,2]
    print full_data2[:,2]

    plt.figure(figsize = (fig_w, 2. / 3 * fig_w)) 
    LMax    =   np.max(full_data2[:,1])/10
    #plt.scatter(full_data1[:,0], full_data1[:,4], color = 'g', label = r'$\hat{v}$=const, $L_{max}=%.2f$nm' %LMax)
    plt.scatter(full_data2[:,0], full_data2[:,4], color = 'r', label = '$\hat{v}$=variab , $L_{max}=%.2f$nm' %LMax)
    plt.legend(frameon = False, loc = 2)
    plt.text(10, .032, r'Experim. $\theta \approx 0.013$')
    plt.title('Buckling, edge=zz')
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
    

plot_kinkOfBend('ac')