'''
Created on 31.8.2015

@author: tohekorh
'''

from twist_topByRod import shearDyn
from misc.lammps_help import get_simulParams
import sys

#width, edge, ratio  =  int(sys.argv[1]), sys.argv[2], int(sys.argv[3])

width, edge, ratio =   8, 'zz', 4
vs      =   [0.0001, 0.0002, 0.0003, 0.0004, 0.0005]   
taito   =   False


params_dic          =   get_simulParams(edge)[-1]
params_dic['taito'] =   taito
params_dic['ratio'] =   ratio
params_dic['width'] =   width
params_dic['vMAX']  =   .0008    
pot_key             =   'LJ'

for _ in range(0, 5):
    for v in vs:
        params_dic['vmax']  =   v
        folder  =   shearDyn(params_dic, pot_key, True)