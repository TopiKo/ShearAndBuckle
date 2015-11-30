'''
Created on 17.8.2015

@author: tohekorh
'''

from misc.ui import query_yes_no, read_simul_params_file, read_stick_simul_params_file
#import numpy as np
import os
import numpy as np
from ase.io.trajectory import PickleTrajectory

def get_datas(ia, T, key):
    
    
    path    =   '/space/tohekorh/ShearSlide/files/%s_%i/%s/' %(ia, T, key) 
    data    =   []
    print path
    for x in os.walk(path):
        
        folder = x[0]
        if folder != path:
            print folder
            if query_yes_no('Take this folder?', default="yes"):
                W   =   int(folder.split('=')[1])
                for filen in os.listdir(folder):
                    if filen[-6:] == '.simul':
                        print filen 
                        if query_yes_no('Take this file?', default="yes"):
                            fPath   =   folder + '/' + filen
                            wf, L, width_i, _, v, _, _, _, _, \
                            _, _, _, _,_, _   =   \
                                        read_simul_params_file(fPath)
                            
                            if width_i != W: 
                                print width_i, W 
                                raise
                            
                            logpath =   fPath[:-6] + '.log'
                            with open(logpath, 'r') as f:
                                lines = f.readlines()
                                for i in range(0, len(lines)):
                                    line = lines[i]
                                    if line[:7] == '# Kink!':
                                        ne = lines[i - 1] # you may want to check that i < len(lines)
                                        shift_b, Rad_b, theta_b  =   ne.split(' ')[1:4]
                                        Rad_b   =   float(L)/float(theta_b)
                                        
                                    elif line[:15] == '# Kink vanished':
                                        ne = lines[i - 1] # you may want to check that i < len(lines)
                                        shift_d, Rad_d, theta_d  =   ne.split(' ')[1:4]
                                        Rad_d   =   float(L)/float(theta_d)
                                        
                                        break
                                data.append([wf, L, v, Rad_b, Rad_d])
    
                                    
    return data

def get_log_data(ia, T, key, Wis = range(1,20)):
    
    path    =   '/space/tohekorh/ShearSlide/files/%s_%i/%s/' %(ia, T, key) 
    data    =   []
    print path
    for x in os.walk(path):
        
        folder = x[0]
        if folder != path:
            print folder
            if query_yes_no('Take this folder?', default="yes"):
                W   =   int(folder.split('=')[1])
                if W in Wis:
                    for filen in os.listdir(folder):
                        if filen[-6:] == '.simul':
                            print filen 
                            if query_yes_no('Take this file?', default="yes"):
                                fPath   =   folder + '/' + filen
                                Wval, L, wi, li, v =   read_simul_params_file(fPath)[:5]
                                
                                
                                if wi != W: raise
                                
                                logpath =   fPath[:-6] + '.log'
                                trajpath=   fPath[:-6] + '.traj'
                                atoms   =   PickleTrajectory(trajpath, 'r')[-1]
                                natoms  =   len(atoms)
                                
                                data.append([wi, li, Wval, L, natoms, v, np.loadtxt(logpath)])
                            
    return data


def get_stick_data(ia, T, key, wi):
    
    path    =   '/space/tohekorh/ShearSlide/files/%s_%i/%s/w=%i' %(ia, T, key, wi) 
    #path    =   '/home/topi/workspace/ShearAndBuckle/data/stick/w=7/'
    
    data    =   []
    print path
    for x in os.walk(path):
        
        folder = x[0]
        if folder != path:
            print folder
            if query_yes_no('Take this folder?', default="yes"):
                r   =   int(folder.split('=')[-1]) 
                if r%10 == 0 or r in [12,14,16,18]:
                    for filen in os.listdir(folder):
                        if filen[-6:] == '.simul':
                            print filen 
                            fPath   =   folder + '/' + filen
                            width_i, length, length_d, temp, _, _, _, _, _, stick  \
                                =   read_stick_simul_params_file(fPath)
                            if stick:
                                print 'stick at L_b/L_t = %.2f' %(length_d/length)
                                data.append([length, length_d])
                            
                            if width_i != wi: raise
                            
                            #trajpath=   fPath[:-6] + '.traj'
                            #atoms   =   PickleTrajectory(trajpath, 'r')[-1]
                            
                        
                        
    return np.array(data)