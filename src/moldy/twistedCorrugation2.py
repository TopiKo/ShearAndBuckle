'''
Created on 23.9.2015

@author: tohekorh
'''
'''
Created on 22.9.2015

@author: tohekorh
'''
import numpy as np
from misc.lammps_help import get_simulParams, get_lammps_params
from potentials.KC_imagBottom import KC_potential_p
from structure import trans_atomsKC
from analyze.plot import plot_posits
from potentials.tests import get_optimal_h
from ase.constraints import FixedLine
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS
from structure import create_stucture
from ase.visualize import view
from ase.structure import graphene_nanoribbon
import matplotlib.pyplot as plt 
import os

path        =   '/space/tohekorh/ShearSlide/files/KC_corrugation/twist2/'
    
def corr_KC(width, edge):
    
    params0 =   get_simulParams(edge)[-1]
    bond    =   params0['bond']
    
    '''
    atoms   =   graphene_nanoribbon(1, 1, type= 'armchair', C_C=bond, saturated = False)
    atoms.rotate([1,0,0], np.pi/2, rotate_cell = True)
    atoms.rotate([0,0,1], -np.pi/2, rotate_cell = True)
    atoms.set_cell([20, 20, 10])
    atoms.center()
    del atoms[[2,3]]
    '''
    atoms   =   create_stucture(2, width, edge, key = 'top', a = bond)[0]
    
    atoms.set_cell([70, 70, 20])
    atoms.center()
    atoms.positions[:,2]    =   3.4
     
    h_t =   []
    for i in range(len(atoms)):
        if atoms[i].number == 1:
            h_t.append(i)
    
    del atoms[h_t]
    
    #view(atoms)
    params  =   {}
    params['positions']         =   atoms.positions
    params['chemical_symbols']  =   atoms.get_chemical_symbols()   
    params['ia_dist']           =   10
    params['edge']              =   edge
    params['bond']              =   bond    
    params['ncores']            =   2
    params['no_edge_neigh']     =   True
    add_KC                      =   KC_potential_p(params)
    
    
    constraints =   []
    for i in range(len(atoms)):
        fix_l   =   FixedLine(i, [0., 0., 1.])
        constraints.append(fix_l)
    
    constraints.append(add_KC)
    
    lamp_parameters =   get_lammps_params(H=False)
    calc            =   LAMMPS(parameters = lamp_parameters) #, files=['lammps.data'])
    atoms.set_calculator(calc)
    atoms.set_constraint(constraints)
    
    #dyn     =   BFGS(atoms, trajectory = 'test.traj')
    #dyn.run(fmax=0.05)
    
    #plot_posits(atoms, edge, bond)
    
    trans_vec   =   trans_atomsKC(atoms.positions[0], edge, bond)
    atoms.translate(trans_vec)
    
    #plot_posits(atoms, edge, bond)
    #exit()
    init_pos    =   atoms.positions.copy()
    r_around    =   init_pos[1]
    
    thetas      =   np.linspace(0, np.pi/3, 61)
    n           =   15
    
    lat_vec1    =   np.array([3./2*bond,  np.sqrt(3)/2*bond, 0.])   
    lat_vec2    =   np.array([3./2*bond, -np.sqrt(3)/2*bond, 0.])
    
    for i, theta in enumerate(thetas):
        fname   =   path + 'corr_%s_theta=%.2f.data' %(edge, theta/(2*np.pi)*360)
        #if not os.path.isfile(fname): 
        print 'Calculating theta = %.2f' %(theta/(2*np.pi)*360)
        atoms.positions =   init_pos
        atoms.rotate([0,0,1], theta, center = r_around)
        
        R               =   np.array([[np.cos(theta), -np.sin(theta), 0.],
                                      [np.sin(theta),  np.cos(theta), 0.],
                                      [0., 0., 1.]])
        
        lat_vec_theta1  =   np.dot(R, lat_vec1.copy())
        lat_vec_theta2  =   np.dot(R, lat_vec2.copy())
                
        #trans_vec1      =   lat_vec_theta1.copy()/n
        trans_vec2      =   lat_vec_theta2.copy()/n
        
        data            =   np.zeros((n,n))
        #plot_posits(atoms, edge, bond, vecs =  [lat_vec_theta1, lat_vec_theta2])
        
        for k in range(n):
            atoms.positions =   init_pos
            atoms.translate(lat_vec_theta1*float(k)/n)
            #plot_posits(atoms, edge, bond, vecs =  [lat_vec_theta1, lat_vec_theta2])
            #print trans_vec1*float(k)/n, k, n, float(k)/n 
            
            for l in range(n):
                atoms.translate(trans_vec2)
                emin, hmin   =   get_optimal_h(atoms, len(atoms), dyn = False)
                #data[k,l,:]  =   [emin, hmin]   
                data[k,l]  =  emin #atoms.get_potential_energy()/len(atoms)
                
        header  =   '%s runs along x-dir, angle measured from x-axis, natoms = %i. x (Angs), e (eV/atom), hmin \n\
the lattice vectors are l1 = [%.5f, %.5f, %.5f] and l2 = [%.5f, %.5f, %.5f], they are divided in %i parts. data[i,j,:] \n\
-> atoms pos += l1/n*i + l2/n*j, initial position is such that atom1 is in middle if hexagon.' \
    %(edge, len(atoms), lat_vec_theta1[0], lat_vec_theta1[1], lat_vec_theta1[2], \
      lat_vec_theta2[0], lat_vec_theta2[1], lat_vec_theta2[2], n)
        np.savetxt(fname, data, header = header)
        
    
def plot_corr_surf(edge):
    
    params0 =   get_simulParams(edge)[-1]
    bond    =   params0['bond']
    
    lat_vec1    =   np.array([3./2*bond,  np.sqrt(3)/2*bond, 0.])   
    lat_vec2    =   np.array([3./2*bond, -np.sqrt(3)/2*bond, 0.])
    
    print lat_vec1
    
    for filen in os.listdir(path):
        if filen[:7] == 'corr_%s' %(edge):
            print filen.split('theta=')[-1].split('.data')[0]
            theta   =   float(filen.split('theta=')[-1].split('.data')[0])/360*2*np.pi
            data    =   np.loadtxt(path + filen)
            data   -=   np.min(data)
            emax    =   np.max(data)
            n,m     =   data.shape
            
            R       =   np.array([[np.cos(theta), -np.sin(theta), 0.],
                                  [np.sin(theta),  np.cos(theta), 0.],
                                  [0., 0., 1.]])
        
            lat_vec_theta1  =   np.dot(R, lat_vec1.copy())
            lat_vec_theta2  =   np.dot(R, lat_vec2.copy())
            X   =   np.zeros((n,m))
            Y   =   np.zeros((n,m))
            Z   =   np.zeros((n,m))
            
            for i in range(n):
                for j in range(m):
                    x,y,_   =   i*lat_vec_theta1/n + j*lat_vec_theta2/n
                    X[i,j]  =   x
                    Y[i,j]  =   y
                    Z[i,j]  =   data[i,j]*1000 #t get meV
    
            origin = 'lower'
            
            xmin, xmax, Dx  =   np.min(X), np.max(X), np.max(X) - np.min(X)
            ymin, ymax, Dy  =   np.min(Y), np.max(Y), np.max(Y) - np.min(Y)
            plt.figure(figsize=(8,6), dpi=80, facecolor='w', edgecolor='k')
            
            plt.arrow(xmin, ymin, 1, .0)
            plt.arrow(xmin, ymin, 1*np.cos(theta), 1*np.sin(theta))
            
            plt.text(xmin - .05*Dx*np.sin(theta), ymin + Dy*(.05 + .05*np.cos(theta)), 'ac-edge of top', rotation=theta/(2*np.pi)*360)
            plt.text(xmin, ymin - .05*Dy, 'ac-edge of bottom')
            CS = plt.contourf(X, Y, Z, 100, # [-1, -0.1, 0, 0.1],
                #alpha=0.5,
                cmap=plt.cm.cool,
                origin=origin)
            
            print theta
            plt.axis('off')
            plt.axis('equal')
            plt.title('Theta = %.1f deg' %(theta/(2*np.pi)*360))
            cbar = plt.colorbar(CS)
            cbar.ax.set_ylabel('corrugation meV/atom')
            plt.show()
            
    '''
    origin = 'lower'

    #origin = 'upper'
    
    delta = .5
    
    x = y = np.arange(-3.0, 3.01, delta)
    X, Y = np.meshgrid(x, y)
    Z1 = plt.mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = plt.mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    Z = 10 * (Z1 - Z2)
    
    print X
    print Y
    print Z
    
    nr, nc = Z.shape
    
    
    
    
    # We are using automatic selection of contour levels;
    # this is usually not such a good idea, because they don't
    # occur on nice boundaries, but we do it here for purposes
    # of illustration.
    CS = plt.contourf(X, Y, Z, emax, # [-1, -0.1, 0, 0.1],
                        #alpha=0.5,
                        cmap=plt.cm.bone,
                        origin=origin)
    
    plt.show()
         
    '''
def study_files(width, edge):
    
    theta_dat   =   []
    for filen in os.listdir(path):
        if filen[:12] == 'corr_w=%02d_%s' %(width, edge):
            print filen.split('theta=')[-1].split('.data')[0]
            theta   =   float(filen.split('theta=')[-1].split('.data')[0])
            data    =   np.loadtxt(path + filen)
            data[:,1] -=    np.min(data[:,1]) 
            
            #plt.plot(data[:,0],data[:,1])
            #plt.show()
            emax    =   np.max(data[:,1])
            theta_dat.append([theta, emax])
    
    # TOPI SORT
    theta_dat   =   np.array(theta_dat) 
    tdat_new    =   np.zeros(theta_dat.shape)
    for i, theta in enumerate(np.sort(theta_dat[:,0])):
        idx = np.where(theta_dat[:,0] == theta)[0][0]
        tdat_new[i,:]   =   [theta, theta_dat[idx, 1]]
        
    return tdat_new

def plot_corr(wedges):
    
    for width, edge in wedges:
        theta_dat   =   study_files(width, edge)
        
        plt.plot(theta_dat[:,0], 1000*theta_dat[:,1], label = 'w=%i, %s' %(width, edge))
        
    plt.xlabel('theta, deg')
    plt.ylabel('corrugation wall, meV/atom')
    plt.legend(frameon = False)
    plt.show()
#corr_KC(5, 'ac')
#plot_corr([[4,'zz'], [5, 'ac']])
plot_corr_surf('ac')
#corr_KC(5, 'ac')

#widths  =   [5,7,9,11]
#for width in widths:
#    corr_KC(width, 'ac') 

#widths  =   [4,6,8,10]
#for width in widths:
#    corr_KC(width, 'zz') 
