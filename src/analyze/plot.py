'''
Created on 17.8.2015

@author: tohekorh
'''
from read import get_datas
from misc.solvers import int_toAngst
import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
fig_w   =   4

def plot_kinkOfBend(edge):
    
    datas       =   get_datas('LJ', edge)
    full_data   =   np.zeros((len(datas),4))
    
    for i, data in enumerate(datas):
        W,L =   data[0], data[1]
        R   =   data[3]
        Wf  =   int_toAngst(W, edge, key = 'width', C_C = 1.42)
        Lf  =   int_toAngst(L, edge, key = 'length', C_C = 1.42)
        
        theta   =   Wf / (2*R)  # new journal of physics Bending and buckling -> theta ~ 0.13.
        print Lf/Wf, theta, L 
        
        full_data[i]    =   [Wf, Lf, R, theta]
    
    print full_data
    plt.figure(figsize = (fig_w, 2. / 3 * fig_w)) 
    plt.scatter(full_data[:,0], full_data[:,3])
    plt.xlabel(r'width \AA')
    plt.ylabel(r'$\theta$')
    plt.show()
    
    

plot_kinkOfBend('zz')