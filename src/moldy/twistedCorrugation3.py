'''
Created on 24.9.2015

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
from write_structure import saveAndPrint
from ase.io.trajectory import PickleTrajectory
from datetime import datetime
import matplotlib.pyplot as plt 
import os

path    =   '/space/tohekorh/ShearSlide/files/KC_corrugation/twist4/'

edge    =   'ac'

params0 =   get_simulParams(edge)[-1]
bond    =   params0['bond']
    
if edge == 'ac':
    lat_vec1    =   np.array([3./2*bond,  np.sqrt(3)/2*bond, 0.])   
    lat_vec2    =   np.array([3./2*bond, -np.sqrt(3)/2*bond, 0.])
elif edge == 'zz':
    lat_vec1    =   np.array([np.sqrt(3.)/2*bond, 3./2*bond, 0.])   
    lat_vec2    =   np.array([np.sqrt(3.)*bond, 0., 0.])

    
def Rmat(theta):
    
    return  np.array([[np.cos(theta), -np.sin(theta), 0.],
                      [np.sin(theta),  np.cos(theta), 0.],
                      [0., 0., 1.]])
        

def corr_KC():
    
    atoms   =   graphene_nanoribbon(1, 1, type= 'armchair', C_C=bond, saturated = False)
    atoms.rotate([1,0,0], np.pi/2, rotate_cell = True)
    
    if edge == 'ac':    
        atoms.rotate([0,0,1], -np.pi/2, rotate_cell = True)
        del atoms[[0,3]]
        trans_idx   =   1
    elif edge == 'zz':
        del atoms[[1,0]]
        trans_idx   =   0
    
    atoms.set_cell([20, 20, 10])
    atoms.center()
    
    params  =   {}
    params['positions']         =   atoms.positions
    params['chemical_symbols']  =   atoms.get_chemical_symbols()   
    params['ia_dist']           =   10
    params['edge']              =   edge
    params['bond']              =   bond    
    params['ncores']            =   2
    add_KC                      =   KC_potential_p(params, True)
    
    
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
    trans_vec   =   trans_atomsKC(atoms.positions[trans_idx], edge, bond)
    atoms.translate(trans_vec)
    
    #plot_posits(atoms, edge, bond)
    init_pos    =   atoms.positions.copy()
    r_around    =   init_pos[trans_idx]
    
    #thetas      =   np.linspace(0, np.pi/3, 7) #, endpoint = False)
    #thetas_deg  =   np.array([1,3,5,7,9,11,12,13,15,17,43,45,47,48,49,51,57,55,57,59])
    thetas_deg  =   np.array([.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5])
    
    traj        =   PickleTrajectory(path + '%s_corr_twist_thetas_(%.1f-%.1f).traj' \
                                     %(edge, np.min(thetas_deg), 
                                       np.max(thetas_deg)), 'w', atoms)

    n           =   100
    
    for i, theta_deg in enumerate(thetas_deg):
        fname   =   path + 'corr_%s_theta=%.2f.data' %(edge, theta_deg)
        print 'Calculating theta = %.2f' %(theta_deg)
        theta   =   theta_deg/360*np.pi*2
        print 'time ' + str(datetime.now().time())
        atoms.positions =   init_pos
        atoms.rotate([0,0,1], theta, center = r_around)
        rot_init_pos    =   atoms.positions.copy()
        
        lat_vec_theta1  =   lat_vec1.copy()
        lat_vec_theta2  =   lat_vec2.copy()

                
        trans_vec2      =   lat_vec_theta2.copy()/n
        
        data            =   np.zeros((n,n))
        
        for k in range(n):
            atoms.positions =   rot_init_pos
            atoms.translate(lat_vec_theta1*float(k)/n)
            #plot_posits(atoms, edge, bond, vecs =  [lat_vec_theta1, lat_vec_theta2])
            print '%.1f percent done' %(100*float(k)/n)
            for l in range(n):
                atoms.translate(trans_vec2)
                emin        =   get_optimal_h(atoms, len(atoms), dyn = False)[0]
                data[k,l]   =   emin #atoms.get_potential_energy()/len(atoms) #emin #
                saveAndPrint(atoms, traj, False)
                
                
        header  =   '%s runs along x-dir, angle measured from x-axis, natoms = %i. x (Angs), e (eV/atom), hmin \n\
the lattice vectors are l1 = [%.5f, %.5f, %.5f] and l2 = [%.5f, %.5f, %.5f], they are divided in %i parts. data[i,j,:] \n\
-> atoms pos += l1/n*i + l2/n*j, initial position is such that atom1 is in middle if hexagon.' \
    %(edge, len(atoms), lat_vec_theta1[0], lat_vec_theta1[1], lat_vec_theta1[2], \
      lat_vec_theta2[0], lat_vec_theta2[1], lat_vec_theta2[2], n)
        np.savetxt(fname, data, header = header)
        
    
def plot_corr_surf():
    
    print lat_vec1
    
    for filen in os.listdir(path):
        if filen[:7] == 'corr_%s' %(edge):
            print filen.split('theta=')[-1].split('.data')[0]
            theta   =   float(filen.split('theta=')[-1].split('.data')[0])/360*2*np.pi
            data    =   np.loadtxt(path + filen)
            data   -=   np.min(data)
            n,m     =   data.shape
            
            
            lat_vec_theta1  =   lat_vec1.copy()
            lat_vec_theta2  =   lat_vec2.copy()
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
            
            #plt.text(xmin - .05*Dx*np.sin(theta), ymin + Dy*(.05 + .05*np.cos(theta)), 'ac-edge of top', rotation=theta/(2*np.pi)*360)
            #plt.text(xmin, ymin - .05*Dy, 'ac-edge of bottom')
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

def get_width_f(width_i):
    
    if edge == 'ac':
        return np.sqrt(3) * bond / 2 * ( width_i - 1 )
    if edge == 'zz':
        return 2 * bond * ( width_i )


def ucell_imagsInUcell(width, theta, length, orig, lat_vec_theta1, lat_vec_theta2, test = False):
    
    
    maps    =   []
    width_f =   get_width_f(width)
    for n1 in range(-length, length):
        for n2 in range(-length, length):
            r   =   lat_vec1 * n1 + lat_vec2 * n2
            if np.abs(r[1]) < width_f / 2 + .05:
                maps.append([n1,n2])
                
    
    
    i       =   0
    mset    =   []
    aset    =   np.zeros((len(maps), 3))
    for n1,n2 in maps:
        r   =   lat_vec_theta1*n1 + lat_vec_theta2*n2 + orig
        
        mat =   np.array([[lat_vec1[0], lat_vec2[0]],
                          [lat_vec1[1], lat_vec2[1]]])
        
        mfloat  =   np.linalg.solve(mat, r[:2])
        
        if mfloat[0] < 0:
            m1  =   -int(np.abs(mfloat[0])) - 1 
        else:
            m1  =   int(mfloat[0])
        if mfloat[1] < 0:
            m2  =   -int(np.abs(mfloat[1])) - 1
        else:
            m2  =   int(mfloat[1])
        
        
        mset.append([m1,m2])   
        aset[i] =   r - (lat_vec1*m1 + lat_vec2*m2)
        
        i      +=   1
    
    #print aset
    if test:
        for m1,m2 in mset:
            pos = lat_vec1 * m1 + lat_vec2 * m2
            if edge == 'ac':
                friend  =   [pos[0] - bond, pos[1]]
            elif edge == 'zz':
                friend  =   [pos[0] + np.sqrt(3)*bond/2, pos[1] + bond/2]

            plt.scatter([pos[0], friend[0]], [pos[1], friend[1]], color = 'red', alpha = .5)
        for n1,n2 in maps:
            pos = lat_vec_theta1*n1 + lat_vec_theta2*n2 + orig
            if edge == 'ac':
                friend  =   [pos[0] + bond*np.cos(theta), pos[1] + bond*np.sin(theta)]
            elif edge == 'zz':
                theta_friend    =   theta + np.pi/6.*7.
                friend  =   [pos[0] + bond*np.cos(theta_friend), 
                             pos[1] + bond*np.sin(theta_friend)]

            plt.scatter([pos[0], friend[0]], [pos[1], friend[1]], alpha = .5)
    
        plt.axis('equal')
        plt.show()
    
    return aset

def get_energy(aset, data):
    
    n, m    =   data.shape
    
    mat =   np.array([[lat_vec1[0], lat_vec2[0]],
                      [lat_vec1[1], lat_vec2[1]]]) 
    
    e   =   np.zeros((len(aset), 3))
    for i, a in enumerate(aset):
        [a1, a2]=   np.linalg.solve(mat, a[:2])
        k, l    =   int(a1*(n-1)), int(a2*(m - 1))
        e[i]    =   [data[k,l]*2, a1, a2] #tot enegy of this cell
        #print data[k,l], k, l
    return e

def get_minimal_orig(data, width, theta, length, lat_vec_theta1, lat_vec_theta2):
    
    def eOfOrig(orig):
        positions_inUcell   =   ucell_imagsInUcell(1, theta, 1, 
                                                   orig, lat_vec_theta1, 
                                                   lat_vec_theta2)
        ncells      =   len(positions_inUcell)
        energy      =   get_energy(positions_inUcell, data)
        return np.sum(energy[:,0])/(ncells*2)
    
    n,m     =   data.shape
    emin    =   1000
    
    for i in range(n):
        for j in range(m):
            orig    =   float(i) / n * lat_vec_theta1 + float(j) / m * lat_vec_theta2 
            e       =   eOfOrig(orig)
            if e < emin:
                emin        =   e
                imin, jmin  =   i, j
    
    return float(imin) / n * lat_vec_theta1 + float(jmin) / m * lat_vec_theta2, e
    
def analyze_ribbonCorr(theta, width = 7, length = 25, test = False, plot_pic = False):
    
    
    for filen in os.listdir(path):
        if filen[:7] == 'corr_%s' %(edge):
            #print filen.split('theta=')[-1].split('.data')[0]
            thetaf   =   float(filen.split('theta=')[-1].split('.data')[0])
            if thetaf == theta:
                data    =   np.loadtxt(path + filen)
                data   -=   np.min(data)
                
                theta           =   float(theta)/360*np.pi*2
                R               =   Rmat(theta)
                lat_vec_theta1  =   np.dot(R, lat_vec1.copy())
                lat_vec_theta2  =   np.dot(R, lat_vec2.copy())

                
                mdir            =   np.array([np.cos(theta + np.pi/2), np.sin(theta + np.pi/2), 0.])
                if edge == 'ac':
                    L           =   np.sqrt(3)*bond*4
                elif edge == 'zz':
                    L           =   3.*bond*2
                orign           =   int(np.max(data.shape)/1.5)
                origins         =   np.zeros((orign,3))
                orig_init, _    =   get_minimal_orig(data, width, theta, length, 
                                                     lat_vec_theta1, lat_vec_theta2) #np.array([0,0,0])
                
                if plot_pic:
                    positions_inUcell \
                                =   ucell_imagsInUcell(width, theta, length, 
                                                       orig_init, lat_vec_theta1, 
                                                       lat_vec_theta2, test = True)
                
                print orig_init
                for i in range(orign):
                    origins[i,:]=   float(i)/orign*L*mdir + orig_init
                
                corrdata        =   np.zeros((orign, 2))
                energy          =   np.ndarray(orign, dtype = 'object')
                for i, orig in enumerate(origins):
                    positions_inUcell \
                                =   ucell_imagsInUcell(width, theta, length, 
                                                       orig, lat_vec_theta1, 
                                                       lat_vec_theta2, test)
                    ncells      =   len(positions_inUcell)
                    energy[i]   =   get_energy(positions_inUcell, data) #E/cell
                    
                    corrdata[i] =   [orig[1], np.sum(energy[i][:,0])/(ncells*2)] #e/atom
                
                
                if plot_pic:
                    esave           =   np.zeros((len(energy[0]), orign))
                    for i, eset in enumerate(energy):
                        esave[:,i]    =   eset[:,0]/2    
                        
                    for eset in esave:
                        plt.plot(corrdata[:,0], eset)
                    plt.show()
                    
                    
                    plt.plot(corrdata[:,0], corrdata[:,1])
                    plt.show()
    return np.max(corrdata[:,1]) - np.min(corrdata[:,1])
                
                
def study_files(width):
    
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

def plot_corr_ribbon():
    

    thetas  =       [0,.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9, 9.5,
                     10, 10.5, 11,11.5,12,13,15,17,
                     20,24,26,28,
                     30,32,34,36,
                     40,43,45,47,48,49,
                     50,51,52,54,55,56,57,58,59,60] #[0, 10, 20 ,30, 40, 50, 60]
    ws  =   [3,5,7,9,11]#[5,7,9,11] #[5, 7, 9, 11, 13, 15, 17, 19, 21]
    ls  =   [3,5,7,9,11,13,15,17,19,21,23,25] #[3,5,7,9,11,13,15] 
    inset_e =   np.zeros((len(ws), len(ls), 3))
    theta_pic   =   10
    n           =   0
    #thetas  =   [10, 15]
    cor_e   =   np.zeros((len(thetas),2))
    coll_data   =   np.zeros((len(ws)*len(ls), 2 + len(thetas)))
    
    for j, w in enumerate(ws):
        for k, l in enumerate(ls):

            for i, theta in enumerate(thetas):
                cor_e[i]    =   [theta, analyze_ribbonCorr(theta, width = w, length = l, test = False)]
                
                if theta == theta_pic:
                    #print j,k,w, l, cor_e[i,1]
                    inset_e[j,k]  =   [w, l, cor_e[i,1]]
            
            
            #print cor_e[:,1]
            coll_data[n, :2] =  w, l
            coll_data[n, 2:] =  cor_e[:,1]
            n           +=  1
        
        plt.plot(cor_e[:,0], cor_e[:,1], '-o')
    
    
    #np.savetxt('/space/tohekorh/ShearSlide/files/corr_data_%s.txt' %edge, 
    #           coll_data, header= 'w, l, max cor on given theta, thetas = ' + str(thetas))
    
    #print np.loadtxt('/space/tohekorh/ShearSlide/files/corr_data_%s.txt' %edge)
    
    plt.xlabel('Rotation Deg')
    plt.ylabel('Corrugation barrier eV/atom')
    plt.title('Corrugation barrier for rotated ribbon \n edge = %s' %edge)
    plt.savefig('/space/tohekorh/ShearSlide/pictures/corr_barr_%s.png' %edge)
    plt.show()
    
    print inset_e
    print inset_e[:,0,2]
    
    #for i in range(len(inset_e[0,:,0])):
    #    plt.plot(inset_e[:,i,0], inset_e[:,i,2], label = inset_e[0,i,1])
    #plt.xlabel('width')
    
    for i in range(len(inset_e[:,0,0])):
        plt.plot(inset_e[i,:,1], inset_e[i,:,2], label = inset_e[i,0,0])
    plt.xlabel('length')
    
    plt.legend(frameon = False)
    plt.title('Corr barrier with different lengths w.r.t width')
    plt.ylabel('corr barrier')
    
    plt.show()

def plot_fromFile():
    thetas  =       [0,1,2,3,4,5,6,7,8,9,
                     10, 11,12,13,15,17,
                     20,24,26,28,
                     30,32,34,36,
                     40,43,45,47,48,49,
                     50,51,52,54,55,56,57,58,59,60]
    data    =   np.loadtxt('/space/tohekorh/ShearSlide/files/corr_data_%s.txt' %edge)
    new_dat =   np.zeros((len(data), 2))
    
    for i, line in enumerate(data):
        w,l     =   line[:2]
        corrs   =   np.array(line[2:])
        
        n = 0
        new_corrs   =   np.zeros(len(corrs) - 9)
        for j in range(len(corrs)):
            if j not in [1,3,5,7,9,11,13,15,17]:
                new_corrs[n]    =   corrs[j]
                n += 1
        
        corrs = new_corrs
        plt.plot(thetas, corrs)
        if w not in [3,5]:
            new_dat[i]  =   [w*l,np.average(corrs)]        
        
        
    plt.show()
        
    plt.scatter(new_dat[:,0], new_dat[:,1])
    plt.xlabel('Area Angst2')
    plt.ylabel('Average corrugation eV')
    plt.show()
        
#

plot_fromFile()
#plot_corr_ribbon()
#exit()
#corr_KC()
#plot_corr([[4,'zz'], [5, 'ac']])
#plot_corr_surf()


#widths  =   [5,7,9,11]
#for width in widths:
#    corr_KC(width, 'ac') 

#widths  =   [4,6,8,10]
#for width in widths:
#    corr_KC(width, 'zz') 
