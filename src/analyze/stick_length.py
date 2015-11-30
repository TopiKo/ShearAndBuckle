'''
Created on 30.10.2015

@author: tohekorh
'''

from moldy.structure import create_bend_stucture, trans_atomsKC
import os
from misc.ui import query_yes_no
from ase.io.trajectory import PickleTrajectory
from ase.visualize import view
from plot import plot_posits
from potentials.KC_imagBottom import make_neighSet, get_setOrig
import numpy as np
import matplotlib.pyplot as plt
pathf   =   '/space/tohekorh/ShearSlide/files/KC_10/ac_stickTaito/'

width   =   7
ratio   =   20
Ld_i    =   1
edge    =   'ac'
bond    =   1.39778034758


for filen in os.listdir(pathf + 'w=%i/r=%i/' %(width, ratio)):
    if filen[-5:] == '.traj':
        print filen 
        if query_yes_no('Take this file?', default="yes"):
            

            bend, straight, [matchL_idx, matchR_idx, vec], [L_bend, L_straight], [left_idxs, right_idxs]\
                =   create_bend_stucture(width, ratio, Ld_i, edge, bond)

            shift_v =   -straight.positions[matchR_idx] + (bend.positions[matchL_idx] + vec) 
            straight.translate(shift_v)
            print matchR_idx
            trans_vec   =    trans_atomsKC(straight.positions[matchR_idx], edge, bond)
            print straight.positions[matchR_idx] + trans_vec
            
            
            rlim        =   .3
            asum_dev    =   .0 * np.array([np.sqrt(3)/2., 1./2])
            m = 0

            y1,y2,y3    =   0, np.sqrt(3)*bond/2, np.sqrt(3)*bond
            nset    =   np.array([[0, y3], [bond, y3], [2*bond, y3], [3*bond, y3],
                         [bond/2, y2], [bond/2 + bond, y2], [bond/2 + 2*bond, y2], 
                         [0, y1], [bond, y1], [2*bond, y1], [3*bond, y1]])
                  
            
            traj    =   PickleTrajectory(pathf + 'w=%i/r=%i/' %(width, ratio) + filen)
            xmin, xmax  =   np.min(traj[0].positions[:,0]), np.max(traj[0].positions[:,0])
            ymin, ymax  =   np.min(traj[0].positions[:,1]), np.max(traj[0].positions[:,1])
            
            for j in range(len(traj)):
                if j%10 == 0:
                    atoms   =   traj[j]
                    coll_ids=   []
                    r_ids   =   []
                    for i, r in enumerate(atoms.positions):
                        if atoms[i].number == 6:
                            ro      =   get_setOrig(r, edge, bond)[:2] + asum_dev
                            alset   =   np.array(ro) + nset 
                            
                            
                            for dev in alset - r[:2]:
                                if np.linalg.norm(dev) < rlim:
                                    coll_ids.append(i)
                                    r_ids.append(np.linalg.norm(dev))
                                    
                                
                            '''
                            if np.abs(ro[1] - r[1]) < rlim or \
                               np.abs(ro[1] + np.sqrt(3)/2*bond - r[1]) < rlim:
                                
                                if  np.abs(ro[0] - r[0]) < rlim or \
                                    np.abs(ro[0] + bond - r[0]) < rlim or \
                                    np.abs(ro[0] + 2*bond - r[0]) < rlim or \
                                    np.abs(ro[0] + bond/2 - r[0]) < rlim or \
                                    np.abs(ro[0] + (1./2 + np.sqrt(3))*bond - r[0]) < rlim:
                                    
                                    coll_ids.append(i)
                                    print i
                                    print get_setOrig(r, edge, bond)
                                    print r  
                                    print get_setOrig(r, edge, bond) - r
                                    
                                if np.abs(get_setOrig(r, edge, bond)[0] + 2*bond - r[0]) < .2 or \
                                   np.abs(get_setOrig(r, edge, bond)[0] + 1*bond - r[0]) < .2:
                                
                                    print i
                                    print get_setOrig(r, edge, bond)
                                    print r  
                                    print get_setOrig(r, edge, bond) - r
                                '''
                    
                    
                    conf_ids    =   []
                    conf_r      =   []
                    for i in range(len(coll_ids)):
                        ids     =   coll_ids[i]
                        for ids2 in coll_ids:
                            if ids != ids2 and ids not in conf_ids:
                                if np.linalg.norm(atoms.positions[ids] - atoms.positions[ids2]) < 1.3*bond:
                                    conf_ids.append(ids)
                                    conf_r.append(r_ids[i])
                                    break
                    
                    rgba_colors = np.zeros((len(conf_ids), 4))
                    # for red the first column needs to be one
                    rgba_colors[:,1] = 1.0
                    # the fourth column needs to be your alphas
                    rgba_colors[:, 3] = .5*np.array(conf_r)/np.max(conf_r)
                    
                    plt.scatter(atoms.positions[:,0], atoms.positions[:,1], color = 'blue')
                    plt.scatter(atoms.positions[conf_ids,0], atoms.positions[conf_ids,1], 
                                color = rgba_colors)
                    
                    #plt.scatter(atoms.positions[coll_ids,0], atoms.positions[coll_ids,1])
                    #plt.show()
                    
                    plt.gca().set_xlim((xmin - 10, xmax + 30))
                    #plt.axis('equal')
                    #plt.axis('off')
                    plt.savefig('/space/tohekorh/Dropbox/Work/ShearAndSlide/corrugationPlots/fig%04d.png' %m)
                    m += 1
                    plt.clf()
                    #plot_posits(atoms, edge, bond)
                    
                    #view(atoms)