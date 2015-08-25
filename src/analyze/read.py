'''
Created on 17.8.2015

@author: tohekorh
'''

from misc.ui import query_yes_no, read_simul_params_file
import numpy as np
import os

def get_datas(ia, edge, btaito = True):
    
    if btaito:  taito = 'taito/'
    else:       taito = ''
    
    path    =   '/space/tohekorh/ShearSlide/files/%s/%s%s/' %(ia, taito, edge) 
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
                            _, L, width_i, _, _, _, _, _, _, \
                            thres_Z, _, _, _, _   =   \
                                        read_simul_params_file(fPath)
                            
                            print thres_Z
                            if width_i != W: raise
                            
                            logpath =   fPath[:-6] + '.log'
                            with open(logpath, 'r') as f:
                                lines = f.readlines()
                                for i in range(0, len(lines)):
                                    line = lines[i]
                                    if line[:7] == '# Kink!':
                                        ne = lines[i - 1] # you may want to check that i < len(lines)
                                        shift, Rad  =   ne.split(' ')[1:3]
                                        data.append([W, L, float(shift), float(Rad)])
                                        break
    
    return data