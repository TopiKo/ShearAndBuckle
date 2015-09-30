'''
Created on 29.9.2015

@author: tohekorh
'''

from study_corrSticking import runAndStudy
from misc.lammps_help import get_simulParams
import sys

#width, edge, ratio  =  int(sys.argv[1]), sys.argv[2], int(sys.argv[3])

width, edge, ratio  =   7, 'zz', 10
taito   =   False


params_dic          =   get_simulParams(edge)[-1]
params_dic['taito'] =   taito
params_dic['ratio'] =   ratio
params_dic['width'] =   width
params_dic['ncores']=   2

pot_key             =   'KC'

for i in range(1, 6):
    params_dic['LdildeL_ratio'] =   .2*i
    folder  =   runAndStudy(params_dic, pot_key, save = True)
        

