'''
Created on 31.8.2015

@author: tohekorh
'''

from twist_topByRod import shearDyn
from misc.lammps_help import get_simulParams
from hysteresis import hyster_study

taito               =   False
widths, edge, ratio =   [7], 'ac', 7


params_dic          =   get_simulParams(edge)[-1]
params_dic['taito'] =   taito
params_dic['ratio'] =   ratio
#params_dic['vmax']  =   .01

#folder = '/space/tohekorh/ShearSlide/files/LJ_10/ac_twistRod/w=5/'
for width in widths:
    params_dic['width'] =   width
    folder  =   shearDyn(params_dic, True)
    #hyster_study(edge, folder)
    
    