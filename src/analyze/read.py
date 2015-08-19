'''
Created on 17.8.2015

@author: tohekorh
'''

from misc.ui import query_yes_no
import numpy as np
import os

def get_datas(ia, edge):
    
    path    =   '/space/tohekorh/ShearSlide/files/%s/%s/' %(ia, edge) 
    data    =   []
    
    for x in os.walk(path):
        
        folder = x[0]
        print folder
        if query_yes_no('Take this folder?', default="yes"):
            W   =   int(folder.split('=')[1])
            for filen in os.listdir(folder):
                if filen[-4:] == '.log':
                    print filen 
                    if query_yes_no('Take this file?', default="yes"):
                        fPath   =   folder + '/' + filen
                        L       =   int(filen.split('=')[1].split('.')[0])
                        #np.loadtxt(folder + '/' + filen)
                        with open(fPath, 'r') as f:
                            lines = f.readlines()
                            for i in range(0, len(lines)):
                                line = lines[i]
                                if line[:7] == '# Kink!':
                                    ne = lines[i - 1] # you may want to check that i < len(lines)
                                    shift, Rad  =   ne.split(' ')[1:3]
                                    data.append([W, L, float(shift), float(Rad)])
                                    break
    
    return data