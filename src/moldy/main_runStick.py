'''
Created on 29.9.2015

@author: tohekorh
'''

from study_corrSticking import runAndStudy
from misc.lammps_help import get_simulParams
from misc.solvers import get_Ld_quess
import sys

#width, edge, ratio, ncores =  int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) 

width, edge, ratio, ncores   =   7, 'ac', 26, 2
taito   =   False


params_dic          =   get_simulParams(edge)[-1]
params_dic['taito'] =   taito
params_dic['ratio'] =   ratio
params_dic['width'] =   width
params_dic['ncores']=   ncores

pot_key             =   'KC'
Ld_apprx            =   get_Ld_quess(width, ratio, edge)


def half_intervalAlg(Ld_apprx):
    
    n           =   10
    acc         =   Ld_apprx / ( 5 * n )
    a_1         =   Ld_apprx / n
    
    int_val     =   a_1
    fdir        =   1
    change_found=   False
    while acc < int_val or .5 < a_1:
        
        print 'Ld/L = %.6f' %(a_1 + fdir*int_val)
        print 'interval = %.6f' %int_val
        print 'Direction = %i' %fdir

        a_1     =   a_1 + fdir*int_val
        params_dic['LdildeL_ratio'] =   a_1   
        stick   =   runAndStudy(params_dic, pot_key, save = True) #test_f(a_1 + fdir*int_val/2) #
        fdir_old=   fdir
        
        if stick:   fdir    =  -1
        else:       fdir    =   1
        
        if fdir != fdir_old:    change_found    =   True    
        
        if change_found:
            int_val     =   int_val/2
        
        print a_1, int_val, fdir
        
half_intervalAlg(Ld_apprx)
            
'''
stick       =   False
i           =   0
maxLd       =   .25
while not stick:
    Ld      =   i * .01
    print 'Ld %.2f' %Ld
    
    params_dic['LdildeL_ratio'] =   Ld
    stick   =   runAndStudy(params_dic, pot_key, save = True)
    if stick or maxLd <= Ld:
        print 'Stick at Ld = %.2f' %Ld
        break
    i   +=  1
    print i
'''