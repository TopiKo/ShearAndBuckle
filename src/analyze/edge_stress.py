'''
Created on 12.10.2015

@author: tohekorh
'''

from ase.structure import graphene_nanoribbon
from ase.visualize import view
from moldy.structure import create_stucture
from ase import Atom, Atoms
import numpy as np
from scipy.optimize import fmin, fmin_cg, minimize
from scipy.misc import derivative
from ase.calculators.lammpsrun import LAMMPS
from ase.optimize import BFGS
from ase.constraints import FixAtoms, FixedPlane
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
                
                                                    
edge        =   'ac'                               
#                             # opm_bond_lammps2 = 1.39778034758 #1.39768096869
bond        =   1.39778034758 # opm_bond_lammps1 = 1.39687536806


if edge == 'ac':    l0  =   3*bond
elif edge == 'zz':  l0  =   np.sqrt(3)*bond
    
fmax        =   .0005#5
perVal      =   .05 #.03    # strech percentage to consider when fitting Y
k           =   21.
dx          =   l0 / 400 #1e-2
maxSteps    =   200
N           =   5 #7
bad_wis     =   [28, 45,61, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97]
parameters  =   {'pair_style':'rebo',
                 'pair_coeff':['* * CH.airebo C H'],
                 'mass'      :['1 12.0', '2 1.0'],
                 'units'     :'metal', 
                 'boundary'  :'p f f'}


def get_atoms(wi):
    if wi % 2 == 0:
        if edge == 'ac':
            atoms   =   graphene_nanoribbon(wi/2,1, type = 'armchair', C_C=bond, saturated = True)
        elif edge == 'zz':
            atoms   =   graphene_nanoribbon(wi,1, type = 'zigzag', C_C=bond, saturated = True)

        atoms.rotate([1,0,0], -np.pi/2, rotate_cell=True)
        atoms.rotate([0,0,1], -np.pi/2, rotate_cell=True)
        
    else:
        if edge == 'ac':
            atoms   =   graphene_nanoribbon(wi/2 + 1, 1, type = 'armchair',  C_C=bond, saturated = True)
        elif edge == 'zz':
            atoms   =   graphene_nanoribbon(wi + 1, 1, type = 'zigzag',  C_C=bond, saturated = True)

        atoms.rotate([1,0,0], -np.pi/2, rotate_cell=True)
        atoms.rotate([0,0,1], -np.pi/2, rotate_cell=True)
        ymax        =   np.max(atoms.positions[:,1])
        delis, hs   =   [], []
        for i in range(len(atoms)):
            if ymax - 2 < atoms.positions[i,1]: delis.append(i) 
            else:
                if atoms[i].number == 1:
                    hs.append(i)
        
        
        atoms_add   =   Atoms()
        for h in hs:
            if edge == 'ac':   add_h_vec    =   [0, (wi-1)*np.sqrt(3)/2*bond + 2*.943, 0]
            elif edge == 'zz': add_h_vec    =   [-np.sqrt(3)/2*bond, (wi*3./2 - 1)*bond + 2*1.09, 0]
            
            atoms_add   +=  Atom('H', position=atoms[h].position + add_h_vec)
        
        del atoms[delis]
    
        atoms   +=  atoms_add
            
    
    
    if edge == 'ac':
        atoms.set_cell([3*bond, np.sqrt(3)/2*(wi + 8)*bond, 8])
    elif edge == 'zz':
        atoms.set_cell([np.sqrt(3)*bond, 3*bond * (wi*2),  8])
    
    ymax    =   np.max(atoms.positions[:,1])
    yav     =   np.average(atoms.positions[:,1])
    
    miny    =   100
    minyidx =   0
    for i, a in enumerate(atoms):
        if np.abs(a.position[1] - yav) < miny:  
            miny    =   np.abs(a.position[1] - yav)
            minyidx =   i
    
    fixId   =   [np.where(atoms.positions == ymax)[0][0], minyidx]
        
    atoms.set_pbc([True, False, False])
    atoms.center() 
    #view(atoms)
    return atoms, fixId



def energy(lx, *args):
    
    atoms, dyn      =   args
    cell_old        =   atoms.cell
    cell_new        =   cell_old
    cell_new[0,0]   =   lx #3*bond*delta
    atoms.set_cell(cell_new, scale_atoms = True)
    
    #print atoms.cell  
    #atoms.center()
    
    dyn.run(fmax=fmax, steps=maxSteps)
    e   =   atoms.get_potential_energy()
    
    return e




def get_YoungAndTau3():
    wis     =   range(4, 120) #[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,25,28,30]
    #wYt     =   np.zeros((len(wis), 4))
    path_f  =   '/space/tohekorh/ShearSlide/pictures/Ytau3/%s/' %edge
    
    wYt =   np.zeros((len(wis) - len(bad_wis), 4))
    j   =   0
    
    for wi in wis:
        if wi not in bad_wis:    
            print 'width %i' %wi
            atoms, fid  =   get_atoms(wi)
            constraints =   []
            constraints.append(FixAtoms(indices = [0]))
            for k in range(len(atoms)):
                constraints.append(FixedPlane(k, [0,0,1]))
            
            atoms.set_constraint(constraints)
            
            # CALCULATOR
            calc    =   LAMMPS(parameters=parameters) 
            atoms.set_calculator(calc)
        
        
            dyn     =   BFGS(atoms, logfile = '/space/tohekorh/log_useless.log', 
                             trajectory='/space/tohekorh/ShearSlide/opm_dx=%.3f_%s.traj' %(dx, edge))
            dyn.run(fmax=fmax, steps=maxSteps)
            
            dedx    =   derivative(energy, l0, dx=dx, n=1, order=21, args = (atoms, dyn))
            
            tau     =   dedx/2.
            cparams.set_tau(tau)
            print 'edge stress %s eV/Angst = %.3f' %(edge, tau)
            
            
            qopm_l  =   (-2*tau/(21*cparams.get_w()) + 1)*l0
            
            energy(qopm_l, atoms, dyn)
            init_pos=   atoms.positions.copy() 
            
            lxs_stre=   np.linspace(1.*qopm_l, (1. + perVal)*qopm_l, N)
            lxs_comp=   np.linspace(1.*qopm_l, (1 - perVal)*qopm_l,  N)
            
            lxs     =   np.zeros(2*N - 1)
            energies=   np.zeros(2*N - 1)
            
            cparams.set_w(wi)
            
            n       =   0
            for lx in lxs_stre:
                print lx/qopm_l
                energies[n] =   energy(lx, atoms, dyn)
                lxs[n]      =   lx
                n          +=   1
            
            atoms.positions =   init_pos
            for lx in lxs_comp[1:]:
                print lx/qopm_l
                energies[n] =   energy(lx, atoms, dyn)
                lxs[n]      =   lx
                n          +=   1
                
            
            eps         =   (lxs - l0)/l0
            
            atoms.positions =   init_pos
            e0      =   energy(l0, atoms, dyn)
            energies-=  e0
            
            
            
            Y_opm   =   curve_fit(energy_teor3, eps, energies, p0=[21])[0][:1]
            
            l_opm   =   -2*tau/(Y_opm*cparams.get_w())*l0 + l0
            eps0    =   l_opm/l0 - 1
            print 'Youngs mod eV/Angst2 = %.3f' %(Y_opm)
            
            
            plt.scatter(eps, energies)
            plt.plot(eps, 2*tau*l0*eps, label = 'tau = %.3f' %tau)
            #plt.plot(eps, dedx*eps*l0, label = 'tau = %.3f' %tau)
            
            plot_eps=   np.linspace(np.min(eps), np.max(eps), 200)
            plt.plot(plot_eps, energy_teor3(plot_eps, Y_opm), '-', color = 'black', 
                    label= 'fit, Y = %.3f eV/angst2' %(Y_opm))
            
    #        plt.plot([l_opm,l_opm], [-.1*np.max(energies), .1*np.max(energies)], '--', 
    #                color = 'black', label = 'eps_opm = %.4f' %(l_opm/l0 - l0))
            plt.plot([0,0], [-.1*np.max(energies), .1*np.max(energies)], '--', 
                    color = 'black')
            plt.plot([eps0, eps0], [-.1*np.max(energies), .1*np.max(energies)], '--', 
                    color = 'black', label = 'eps0 = %.5f' %(eps0))
            plt.legend(frameon = False)
            #plt.show()
            plt.savefig(path_f + 'w=%i_dx=%.3f.png' %(wi, dx))
            plt.clf()
            
            wYt[j]  =   [cparams.get_w(), Y_opm, tau, eps0]
            
            
            plt.scatter(wYt[:j+1,0], wYt[:j+1,1], label = 'Y', color = 'black')
            plt.legend(loc = 2, frameon = False)
            plt.twinx()
            plt.plot(wYt[:j+1,0], wYt[:j+1,2], '-o', label = 'tau')
            plt.xlabel('width')
            plt.title('Young s modulus and edge stress')
            plt.legend(frameon = False)
            plt.savefig(path_f  + 'Ytau.png')
            plt.clf()
    
            np.savetxt(path_f  + 'Ytau.txt', wYt, 
                       header = 'width, Y, tau, eps_opm')
    
            #plt.show()   
            j += 1
        
    

class params():
    
    def __init__(self, eps0, wi):
        self.eps0   =   eps0
        self.set_w(wi)
    
    def set_eps0(self, eps0):
        self.eps0   =   eps0
    def set_tau(self, tau):
        self.tau   =   tau
    def set_w(self, wi):
        if edge == 'ac':
            self.w   =   (wi - 1)*np.sqrt(3)/2*bond
        if edge == 'zz':
            self.w   =   wi*3./2*bond - bond

    def get_eps0(self):
        return self.eps0 
    def get_tau(self):
        return self.tau
    def get_w(self):
        return self.w        
        
def energy_teor(eps, Y, eps0):
    
    #eps0, w =   cparams.get_eps0(), cparams.get_w()
    w   =   cparams.get_w()
    return (.5*eps**2 - eps0*eps)*l0*w*Y

def energy_teor2(eps, Y):
    
    eps0, w =   cparams.get_eps0(), cparams.get_w()
    return (.5*eps**2 - eps0*eps)*l0*w*Y

def energy_teor3(eps, Y):
    
    tau, w =   cparams.get_tau(), cparams.get_w()
    return (.5*eps**2 + 2*tau/(Y*w)*eps)*l0*w*Y

def oneOverR(R, a, b):
    
    return a/R + b

def get_YoungAndTau():
    
    wis =   range(5, 120)
    wYt =   np.zeros((len(wis) - len(bad_wis), 4))
    j   =   0
    
    for wi in wis:
        if wi not in bad_wis:
            print wi, edge
            atoms, fid       =   get_atoms(wi)
            constraints =   []
            constraints.append(FixAtoms(indices = [0]))
            for k in range(len(atoms)):
                constraints.append(FixedPlane(k, [0,0,1]))
            #constraints.append(FixedPlane(fid[0], [0,0,1]))
            #constraints.append(FixedPlane(fid[1], [0,0,1]))
        
            atoms.set_constraint(constraints)
            
            # CALCULATOR
            calc    =   LAMMPS(parameters=parameters) 
            atoms.set_calculator(calc)
        
        
            dyn     =   BFGS(atoms, 
                             trajectory='/space/tohekorh/ShearSlide/opm_dx=%.3f_%s.traj' %(dx, edge),
                             maxstep=.01, logfile = '/space/tohekorh/log_useless.log')
            dyn.run(fmax=fmax, steps=10*maxSteps)
            
            qopm_l  =   (-2*(-.2)/(21*cparams.get_w()) + 1)*l0
            
            energy(qopm_l, atoms, dyn)
            init_pos=   atoms.positions.copy()
    
            
            lxs_stre=   np.linspace(1.*qopm_l, (1. + perVal)*qopm_l, N)
            lxs_comp=   np.linspace(1.*qopm_l, (1 - perVal)*qopm_l,  N)
            
            lxs     =   np.zeros(2*N - 1)
            energies=   np.zeros(2*N - 1)
            
            cparams.set_w(wi)
            
            n       =   0
            for lx in lxs_stre:
                print lx/qopm_l
                energies[n] =   energy(lx, atoms, dyn)
                lxs[n]      =   lx
                n          +=   1
            
            atoms.positions =   init_pos
            for lx in lxs_comp[1:]:
                print lx/qopm_l
                energies[n] =   energy(lx, atoms, dyn)
                lxs[n]      =   lx
                n          +=   1
            
            energies=   energies - np.min(energies)
            eps     =   (lxs - l0)/l0
            
            Y_opm, eps0 =   curve_fit(energy_teor, eps, energies, p0=[21, qopm_l/l0 - 1])[0][:2]
            tau     =   -Y_opm*cparams.get_w()*eps0/2.
            
            plt.scatter(eps, energies)
            plt.plot([eps0, eps0], [-.1*np.max(energies), .1*np.max(energies)], '--', 
                     color = 'black', label = 'eps_opm = %.4f' %eps0)
        
            
            plot_eps=   np.linspace(np.min(eps), np.max(eps), 200)
            plt.plot(plot_eps, energy_teor(plot_eps, Y_opm, eps0), '-', color = 'black', 
                     label= 'fit, Y = %.3f eV/angst2, \n tau = %.3f eV/Angst' %(Y_opm, tau))
            plt.legend(frameon = False)
            plt.xlabel('eps')
            plt.ylabel('energy eV')
            plt.title('Energy and sheet elasticity fit to it')
            
            wYt[j]  =   [cparams.get_w(), Y_opm, tau, eps0]
            plt.savefig('/space/tohekorh/ShearSlide/pictures/Ytau/%s/w=%i.png' %(edge, wi))
            plt.clf()
        
            #plt.show()
    
            plt.scatter(wYt[:j+1,0], wYt[:j+1,1], label = 'Y', color = 'black')
            plt.legend(loc = 2, frameon = False)
            plt.twinx()
            plt.plot(wYt[:j+1,0], wYt[:j+1,2], '-o', label = 'tau')
            plt.xlabel('width')
            plt.title('Young s modulus and edge stress')
            plt.legend(frameon = False)
            plt.savefig('/space/tohekorh/ShearSlide/pictures/Ytau/%s/Ytau.png' %edge)
            plt.clf()
    
            np.savetxt('/space/tohekorh/ShearSlide/pictures/Ytau/%s/Ytau.txt' %edge, wYt, 
                       header = 'width, Y, tau, eps_opm')
            
            j += 1
def get_YoungAndTau2():
    
    wis =   range(5, 120)
    
    path_f  =   '/space/tohekorh/ShearSlide/pictures/Ytau2/%s/' %edge
    wYt =   np.zeros((len(wis) - len(bad_wis), 5))
    j   =   0
    fa,fb   =   0.09, 0
    
    for wi in wis:
        if wi not in bad_wis:
            print wi, edge
            atoms, fid       =   get_atoms(wi)
            constraints =   []
            constraints.append(FixAtoms(indices = [0]))
            for k in range(len(atoms)):
                constraints.append(FixedPlane(k, [0,0,1]))
            #constraints.append(FixedPlane(fid[0], [0,0,1]))
            #constraints.append(FixedPlane(fid[1], [0,0,1]))
        
            atoms.set_constraint(constraints)
            
            # CALCULATOR
            calc    =   LAMMPS(parameters=parameters) 
            atoms.set_calculator(calc)
        
        
            dyn     =   BFGS(atoms, 
                             trajectory='/space/tohekorh/ShearSlide/opm_dx=%.3f_%s.traj' %(dx, edge),
                             maxstep=.01, logfile = '/space/tohekorh/log_useless.log')
            dyn.run(fmax=fmax, steps=10*maxSteps)
            cparams.set_w(wi)
            
            #derivative(energy, l0, dx=dx, n=1, order=5)
            
            #opm_l   =   fmin(energy, l0, xtol = 1e-6, ftol = 1e-6, args = (atoms, dyn))[0]
            cont    =   True
            cont    =   False
            try:
                #opm_l   =   fmin_cg(energy, l0, args=(atoms, dyn), gtol=1e-03, 
                #                    epsilon=1.4901161193847656e-03, 
                #                    full_output=0, disp=1)
                #qopm_l  =   (-2*(-.2)/(21*cparams.get_w()) + 1)*l0
                qopm_l  =   oneOverR(cparams.get_w(), fa, fb)*l0 + l0
                #qopm_l  =   l0
                
                res     =   minimize(energy, [qopm_l], args=(atoms, dyn), method='L-BFGS-B', 
                                     jac=None, bounds=[(l0*.995, l0*1.005)], tol=None, 
                                     callback=None,  
                                     options={'disp': None, 'iprint': -1, 
                                              'gtol': 1e-012*wi, 'eps': 1e-04, 
                                              'maxiter': 15,
                                              'ftol': 2.220446049250313e-09, 
                                              'maxcor': 10, 
                                              'maxfun': 100}) #['x'][0] 'factr': 10000
                
                '''                
                options={'disp': None, 'iprint': -1, 
                                              'gtol': 1e-03*wi, 'eps': 1e-04, 
                                              'maxiter': 15,
                                              'ftol': 2.220446049250313e-09, 
                                              'maxcor': 10, 
                                              'maxfun': 100})
                '''
                opm_l   =   res['x'][0]
                eopm    =   res['fun']
                print 'number of iterations = %i' %res['nit']
                init_pos=   atoms.positions.copy()
                    
                cont    =   True
            except Exception as e:
                print e
            
            if cont:
                eps0    =   (opm_l - l0)/l0 
                print wi, opm_l, l0, opm_l - l0, eps0 
                #qopm_l  =   (-2*(-.2)/(21*cparams.get_w()) + 1)*l0
                
                
                lxs_stre=   np.linspace(1.*opm_l, (1. + perVal)*opm_l, N)
                lxs_comp=   np.linspace(1.*opm_l, (1 - perVal)*opm_l,  N)
                
                lxs     =   np.zeros(2*N - 1)
                
                energies=   np.zeros(2*N - 1)
                
                cparams.set_eps0(eps0)
                
                n       =   0
                for lx in lxs_stre:
                    print lx/opm_l
                    energies[n] =   energy(lx, atoms, dyn)
                    lxs[n]      =   lx
                    n          +=   1
                
                atoms.positions =   init_pos
                for lx in lxs_comp[1:]:
                    print lx/opm_l
                    energies[n] =   energy(lx, atoms, dyn)
                    lxs[n]      =   lx
                    n          +=   1
                
                #min_e   =   energy(opm_l, atoms, dyn)
                energies=   energies - np.min(energies)#- min_e
                eps     =   (lxs - l0)/l0
                
                
                Y_opm   =   curve_fit(energy_teor2, eps, energies, p0=[21])[0][:1]
                tau     =   -Y_opm*cparams.get_w()*eps0/2
                tau2    =   (eopm - energy(l0, atoms, dyn))/(eps0*l0) 
                plt.scatter(eps, energies)
                plt.plot([eps0, eps0], [-.1*np.max(energies), .1*np.max(energies)], '--', 
                         color = 'black', label = 'eps_opm = %.4f' %eps0)
            
                
                plot_eps=   np.linspace(np.min(eps), np.max(eps), 200)
                plt.plot(plot_eps, energy_teor(plot_eps, Y_opm, eps0), '-', color = 'black', 
                         label= 'fit, Y = %.3f eV/angst2, \n tau = %.3f eV/Angst' %(Y_opm, tau))
                plt.legend(frameon = False)
                plt.xlabel('eps')
                plt.ylabel('energy eV')
                plt.title('Energy and sheet elasticity fit to it')
                
                wYt[j]  =   [cparams.get_w(), Y_opm, tau, eps0, tau2]
                plt.savefig(path_f  + 'w=%i.png' %wi)
                plt.clf()
            
                #plt.show()
                if j > 2:
                    fa, fb  =   curve_fit(oneOverR, wYt[:j,0], wYt[:j,3], p0=[1, 0])[0][:2]
                    #if j%5 == 0:
                    #    plt.scatter(wYt[:j,0], wYt[:j,3])
                    #    plt.plot(wYt[:j,0], oneOverR(wYt[:j,0], fa, fb))
                    #    plt.show()
                    #print fa, fb
                
                plt.scatter(wYt[:j+1,0], wYt[:j+1,1], label = 'Y', color = 'black')
                plt.legend(loc = 2, frameon = False)
                plt.twinx()
                plt.plot(wYt[:j+1,0], wYt[:j+1,2], '-o', label = 'tau')
                plt.plot(wYt[:j+1,0], wYt[:j+1,4], '-o', label = 'tau2')
                plt.xlabel('width')
                plt.title('Young s modulus and edge stress')
                plt.legend(frameon = False)
                plt.savefig(path_f  + 'Ytau.png')
                plt.clf()
        
                #np.savetxt(path_f  + 'Ytau.txt', wYt, 
                #           header = 'width, Y, tau, eps_opm')
                
                j   +=  1

class plotParams():
    
    def __init__(self):
        pass
    
    def set_params(self, a,b):
        self.a  =   a
        self.b  =   b
    def get_params(self):
        return self.a, self.b
        

def epsF(w, a, b):
    
    fay, fby    =   pparams.get_params()
    return -2*(a/w + b)/(fay + fby*w)

def plot_YoungAndTau(edge):

    yWt =   np.loadtxt('/space/tohekorh/ShearSlide/pictures/Ytau2/%s/Ytau.txt' %edge)
    
    if np.min(yWt[:,0]) == 0:
        zeros   =   np.where(yWt[:,0] == 0)[0][0]
        yWt     =   yWt[:zeros]
    
    w, Y, t, e, t2  =   [],[],[],[], []
    for val in yWt:
        if 19 < val[1] < 23:
            if len(Y) > 0:
                if val[0] < 40:
                    if np.abs(val[1] - Y[-1]) < .5:
                        w.append(val[0])  
                        Y.append(val[1])
                        t.append(val[2])
                        e.append(val[3])
                        #t2.append(val[4])
                else:
                    if np.abs(val[1] - Y[-1]) < .1:
                        w.append(val[0])  
                        Y.append(val[1])
                        t.append(val[2])
                        e.append(val[3])
                        #t2.append(val[4])
            else:
                w.append(val[0])  
                Y.append(val[1])
                t.append(val[2])
                e.append(val[3])
                #t2.append(val[4])
    ws      =   np.array(w)
    Ys      =   np.array(Y)
    ts      =   np.array(t)
    es      =   np.array(e)
    
    
    fay, fby=   curve_fit(oneOverR, ws, Ys, p0=[1, 0])[0][:2]
    pparams.set_params(fay, fby)
    a,b     =   curve_fit(epsF, ws, es, p0=[.01, 0])[0][:2]

    Ysn     =   oneOverR(ws, fay, fby)
    
    #tau rom 1/(wY) fit
    def oneOverX(x, a):
        return a/x
    
    testa   =   curve_fit(oneOverX, ws*Ysn, es, p0=[.01])[0][0]
    print testa
    plt.plot(1/(ws*Ysn), es)
    plt.plot(1/(ws*Ysn), testa*1/(ws*Ysn))
    print 'tau = %.4feV/Angst' %(-testa/2)
    plt.show()
    
    #print fay, fby
    print a, b
    print fay, fby
    esn     =   epsF(ws, a, b) 
    
    #print tau
    tsn     =   -esn*ws*Ysn/2
    
    
    plt.scatter(ws, Ys, label='Y') 
    plt.plot(ws, Ysn, '-')
    plt.ylabel('Youngs mod eV/Angst2')
    plt.xlabel('width Angst')
    
    plt.legend(frameon = False, loc = 2)
     
    plt.twinx()
    #plt.plot(ws, tsn, '-', label='tau', color = 'black')
    #####plt.plot(ws, t2, '-o', label='tau', color = 'black')
    
    #plt.scatter(1/(ws*Ysn), es, label='eps', color = 'black')
    #plt.plot(1/(ws*Ysn), esn, '-', label='eps', color = 'black')
    plt.scatter(ws, es, label='eps', color = 'black')
    plt.plot(ws, esn, '-', label='eps', color = 'black')
    
    plt.legend(frameon = False)
    plt.ylabel('Tau eV/Angst')
    plt.show() 
    
pparams     =   plotParams()
cparams     =   params(1, 5)
#get_YoungAndTau2()
plot_YoungAndTau('ac')
#get_YoungAndTau()
#get_YoungAndTau2()
#get_YoungAndTau3()
