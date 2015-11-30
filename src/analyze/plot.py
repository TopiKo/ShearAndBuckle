'''
Created on 17.8.2015

@author: tohekorh
'''
from read import get_datas, get_log_data, get_stick_data
from misc.solvers import int_toAngst
import numpy as np
import matplotlib.pyplot as plt
from potentials.KC_imagBottom import make_neighSet, get_setOrig
from scipy.integrate import dblquad
from scipy.optimize import fmin
from ase.utils.eosase2 import curve_fit

plt.rc('text', usetex=True)
fig_w   =   6
taito   =   True
T       =   10
k       =   21
tau     =   -.2
    

def plot_posits(atoms, edge, bond, vecs=None):
    
    positions   =   atoms.positions
    pos_used    =   []
    n_set       =   make_neighSet(15, edge, bond)
            
    for i, r in enumerate(positions):
        if atoms[i].number == 6:
            norms   =   [np.linalg.norm(x) for x in positions[pos_used] - r]
            if len(norms) != 0: 
                if 12 < np.min([np.linalg.norm(x) for x in positions[pos_used] - r]):  
                    neigh_set   =   get_setOrig(r, edge, bond) + n_set
                    pos_used.append(i)
                    plt.scatter(neigh_set[:,0], neigh_set[:,1], color = 'red', alpha = .5)
            else:
                neigh_set   =   get_setOrig(r, edge, bond) + n_set
                pos_used.append(i)
                plt.scatter(neigh_set[:,0], neigh_set[:,1], color = 'red', alpha = .5)
            
    plt.scatter(positions[:,0], positions[:,1], color = 'black')
    
    if vecs != None:
        for vec in vecs:
            plot_vec    =   np.array([positions[0],
                                      positions[0] + vec])  
            plt.plot(plot_vec[:,0], plot_vec[:,1])
    plt.axis('equal')
    plt.show()
            
def plot_kinkOfBend(edge):
    
    datas   =   get_datas('LJ', T, '%s_twistTaito' %edge)
    data    =   np.array(datas)
    
    plt.figure(figsize = (fig_w, 2. / 3 * fig_w)) 
    theta_b =   data[:,0] / (2 * data[:,3])
    theta_v =   data[:,0] / (2 * data[:,4])
    v       =   data[:,2]/np.max(data[:,2])
    
    LMax    =   np.max(data[:,1])/10
    
    
    set_label   =   True
    for i in range(len(v)):
        if v[i] == 1 and set_label:
            set_label   =   False
            plt.scatter(data[i,0], theta_b[i], color = 'b', alpha = v[i], 
                        label = r'bucle born, $L_{max}=%.2f$nm' %LMax)
            plt.scatter(data[i,0], theta_v[i], color = 'r', alpha = v[i], 
                        label = r'bucle vanish')
        else:
            plt.scatter(data[i,0], theta_b[i], color = 'b', alpha = v[i])
            plt.scatter(data[i,0], theta_v[i], color = 'r', alpha = v[i])
    plt.plot(data[:,0], np.ones(len(data[:,0]))*.04, '-.', color = 'black')
    plt.text(np.min(data[:,0]), .041, 'teor')
    plt.plot(data[:,0], np.ones(len(data[:,0]))*.013, '-.', color = 'black')
    plt.text(np.min(data[:,0]), .014, r'Experim. $\theta \approx 0.013$')
    
    
    plt.legend(frameon = False, loc = 1)
    #plt.text(10, .062, r'Experim. $\theta \approx 0.013$')
    plt.title('Buckling, edge=%s, T=%iK' %(edge, T))
    plt.xlabel(r'width \AA')
    plt.ylabel(r'$\theta$')
    plt.show()

def plot_kinkOfBend_deg(edge):
    
    plot_mus    =   True #False #
    magic_par   =   0.0228171 #0.0203724
    datas_ac    =   get_log_data('LJ', T, '%s_twistTaito' %'ac')
    datas_zz    =   get_log_data('LJ', T, '%s_twistTaito' %'zz')
    
    coll_data_ac=   np.zeros((len(datas_ac), 6))    
    coll_data_zz=   np.zeros((len(datas_zz), 6))    

    thresz      =   5
    vmax    =   -1.
    eps     =   0.002843732471143 #[eV]
    kappa   =   .95 #[eV]
    n       =   1/(3*np.sqrt(3)*bond**2/4)      
    pot     =   1
    
    def theta_teor(W):

        Y   =   11.43/W + 19.88
        return 24./Y*np.sqrt(eps*kappa*np.pi)*n
    
    class alpha():
        
        def __init__(self):
            pass
        
        def set_s(self, s):
            self.s  =   s
        def set_t(self, t):
            self.t  =   t
        def set_pot(self, pot):
            self.pot=   pot
        def set_params(self, s,t, pot):
            self.set_s(s)
            self.set_t(t)
            self.set_pot(pot)

        def get_alpha(self, w):
            #Y       =   11.43/W + 19.88
            #A, B    =   48.7, .0055
            #Length  =   12.8 

            #B   =   4/(Y*np.pi**2)*(A/Length**2 + B*Length**2)
            return magic_par/(self.s/w**self.pot + self.t)
            #return self.s/w**self.pot + self.t 
    
    
    alpha   =   alpha()
    def theta_teorFitPar(W, *args):
        
        if len(args) == 1:
            s       =   args[:1]
            alpha.set_s(s)
            
        if len(args) == 2:
            s,t     =   args[:2]
            alpha.set_s(s)
            alpha.set_t(t)
        if len(args) == 3:
            s,t,pot =   args[:3]
            alpha.set_s(s)
            alpha.set_t(t)
            alpha.set_pot(pot)
        
        #Y       =   11.43/W + 19.88
        #A, B    =   48.7, .0055
        #Length  =   12.8 
        #print np.average(4/(Y*np.pi**2)*(A/Length**2 + B*Length**2))
        return magic_par/alpha.get_alpha(W)
        #return 4/(Y*np.pi**2*alpha.get_alpha(W))*(A/Length**2 + B*Length**2)
    
    def overR(W, a1, a2):
        
        return a1/W + a2
    
    def overRandconst(W, a1, a2):
        
        return a1/(W + a2)

    def overR2(W, a1, a2):
        
        return a1/W**2 + a2

    
    for i, data in enumerate(datas_ac):
        energy_table    =   data[6]
        v               =   data[5]
        natoms          =   data[4]
        Wi, Li, W, L    =   data[:4]
        phi             =   energy_table[:,3]
        epot            =   (energy_table[:,5] - energy_table[-1,5])/natoms
        
        
        heights     =   energy_table[:,4]
        inds        =   np.where(thresz < heights)[0]
        indv        =   max(inds)
        indb        =   min(inds)
        
        W           =   int_toAngst(Wi, 'ac', key='width', C_C = bond)
        
        thetas      =   W/(2*(L/phi))
        mus         =   thetas*2./W*360/(2*np.pi)*10 
        
        coll_data_ac[i]    =   [W, v, thetas[indv], thetas[indb], mus[indv], mus[indb]]
        
        if v > vmax: 
            vmax   = v
            
    for i, data in enumerate(datas_zz):
        energy_table    =   data[6]
        v               =   data[5]
        natoms          =   data[4]
        Wi, Li, W, L    =   data[:4]
        phi             =   energy_table[:,3]
        epot            =   (energy_table[:,5] - energy_table[-1,5])/natoms
        
        
        heights     =   energy_table[:,4]
        inds        =   np.where(thresz < heights)[0]
        indv        =   max(inds)
        indb        =   min(inds)

        W           =   int_toAngst(Wi, 'zz', key='width', C_C = bond)
        thetas      =   W/(2*(L/phi))
        mus         =   thetas*2./W*360/(2*np.pi)*10 
        
        coll_data_zz[i]    =   [W, v, thetas[indv], thetas[indb], mus[indv], mus[indb]]
        
        if v > vmax: 
            vmax   = v
        
    
    if plot_mus:
        for i in range(len(coll_data_ac)):
            plt.scatter(coll_data_ac[i,0], coll_data_ac[i,4], 
                        alpha = coll_data_ac[i,1]/vmax, color = 'red')
            plt.scatter(coll_data_ac[i,0], coll_data_ac[i,5], 
                        alpha = coll_data_ac[i,1]/vmax, color = 'blue')

        #for i in range(len(coll_data_zz)):
        #    plt.scatter(coll_data_zz[i,0], coll_data_zz[i,4], 
        #                alpha = coll_data_zz[i,1]/vmax, color = 'red')
        #    plt.scatter(coll_data_zz[i,0], coll_data_zz[i,5], 
        #                alpha = coll_data_zz[i,1]/vmax, color = 'blue')
            
        alpha.set_pot(1)
        alpha.set_t(magic_par)
    
        sopm1_ac       =   curve_fit(theta_teorFitPar, coll_data_ac[:,0], coll_data_ac[:,2], 
                                  p0 = [1.])[0][:1]
        
        plot_widths = np.linspace(np.min(coll_data_ac[:,0]) - .5, np.max(coll_data_ac[:,0]) + 2, 50)
        
        plt.plot(plot_widths, 2/plot_widths*3600/(2*np.pi)*theta_teorFitPar(plot_widths, sopm1_ac, magic_par, 1.), '--', 
                 label = r'$\Theta \approx \frac{%.2f}{w} + %.3f$' %(sopm1_ac, magic_par),
                 color = 'green')
        
        plt.legend(loc = 2, frameon = False)
        plt.ylabel(r'$^{\circ}/$\AA')
        plt.xlabel(r'Width \AA')
        
        #plt.twinx()
        # 
        #alpha.set_params(sopm1_ac, magic_par, 1)
        #plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
        #         label = r'$\alpha = \frac{%.3f}{%.2f/w + %.3f} \rightarrow %.3f$' 
        #         %(magic_par, sopm1_ac, magic_par, 1),
        #         color = 'green')
        #plt.legend(loc  = 1, frameon = False)
        #plt.ylabel(r'$\alpha$')
        plt.show()
            
    else:
            
        for i in range(len(coll_data_ac)):
            plt.scatter(coll_data_ac[i,0], coll_data_ac[i,2], 
                        alpha = coll_data_ac[i,1]/vmax, color = 'red')
            plt.scatter(coll_data_ac[i,0], coll_data_ac[i,3], 
                        alpha = coll_data_ac[i,1]/vmax, color = 'blue')
        
        
        #pot = 1.5
        #alpha.set_pot(pot)
        #sopm15, topm15  =   curve_fit(theta_teorFitPar, coll_data[:,0], coll_data[:,2], p0 = [1., 1])[0][:2]
        
        pot = 2.
        alpha.set_pot(pot)
        sopm2, topm2    =   curve_fit(theta_teorFitPar, coll_data_ac[:,0], coll_data_ac[:,2], p0 = [1., 1])[0][:2]
        
        alpha.set_pot(1)
        alpha.set_t(magic_par)
    
        sopm1       =   curve_fit(theta_teorFitPar, coll_data_ac[:,0], coll_data_ac[:,2], p0 = [1.])[0][:1]
        
        
        sopm3, topm3, potopm    =   curve_fit(theta_teorFitPar, coll_data_ac[:,0], coll_data_ac[:,2], p0 = [1., 1, 1])[0][:3]
        
        
        
        #plt.plot(coll_data[:,0], np.ones(len(coll_data[:,0]))*.04, '-.', color = 'black')
        #plt.text(np.min(coll_data[:,0]), .041, 'teor')
        
        #plt.scatter(coll_data[:,0], theta_teor(coll_data[:,0]))
        #plt.text(np.min(coll_data[:,0]), .041, 'teor')
    
    
        plot_widths = np.linspace(np.min(coll_data_ac[:,0]) - .5, np.max(coll_data_ac[:,0]) + 2, 50)
        
        #plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm15, topm15, 1.5), '--', 
        #         label = 'alpha = %.2f/w**%.1f + %.3f' %(sopm15, 1.5, topm15))
        
        '''
        plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm2, topm2, 2.), '--', 
                 label = r'$\Theta \approx \frac{%.2f}{w^2} + %.3f$' %(sopm2, topm2),
                 color = 'red')
        '''
        plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm1, magic_par, 1.), '--', 
                 label = r'$\Theta \approx \frac{%.2f}{w} + %.3f$' %(sopm1, magic_par),
                 color = 'green')
        
        #plt.plot(plot_widths, overR2(plot_widths, sopm2a, topm2a), 
        #         label = r'$\Theta \approx \frac{1}{%.2fw^2} + %.3f$' %(sopm2a, topm2a), )
        '''
        plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm3, topm3, potopm), '.-', 
                 label = r'$\Theta \approx \frac{%.2f}{w^{%.2f}} + %.3f$' %(sopm3, potopm, topm3),
                 color = 'red')    
        '''
        plt.scatter(8.46, .026, marker = 'D', color = 'black')
        plt.text(np.min(coll_data_ac[:,0]), .028, r'Experim. $\Theta \approx 0.026 = 4deg/nm$')
        plt.legend(frameon = False, loc = 2)
        plt.xlabel('width Angst')
        plt.ylabel(r'$\Theta$')
        plt.twinx()
        
        #alpha.set_params(sopm15, topm15, 1.5)
        #plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
        #         label = 'alpha -> %.3f' %alpha.get_alpha(10000))
        
        '''
        alpha.set_params(sopm2, topm2, 2)
        plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
                 label = r'$\alpha = \frac{%.3f}{%.2f/w^2 + %.3f} \rightarrow %.3f$' 
                 %(magic_par, sopm2, topm2, magic_par/topm2),
                 color = 'red')
        '''
        alpha.set_params(sopm1, magic_par, 1)
        plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
                 label = r'$\alpha = \frac{%.3f}{%.2f/w + %.3f} \rightarrow %.3f$' 
                 %(magic_par, sopm1, magic_par, 1),
                 color = 'green')
        '''
        print sopm1/magic_par, sopm1, magic_par
        plt.plot(plot_widths, plot_widths/(sopm1/magic_par + plot_widths), '-o',               
                 color = 'green')
        
        
        alpha.set_params(sopm3, topm3, potopm)
        plt.plot(plot_widths, alpha.get_alpha(plot_widths), '.-', 
                 label = r'$\alpha = \frac{%.3f}{%.2f/w^{%.2f} + %.3f} \rightarrow  %.3f$' 
                 %(magic_par, sopm3, potopm, topm3, magic_par/topm3),
                 color = 'red')
        '''
        
        plt.ylabel(r'$\alpha$')
        plt.legend(frameon = False)
        plt.show()

def plot_kinkOfBend3(edge):
    
    plot_mus    =   False
    magic_par   =   0.0203724
    datas       =   get_log_data('LJ', T, '%s_twistTaito' %edge)
    coll_data   =   np.zeros((len(datas), 6))    
    thresz      =   5
    vmax    =   -1.
    eps     =   0.002843732471143 #[eV]
    kappa   =   .95 #[eV]
    n       =   1/(3*np.sqrt(3)*bond**2/4)      
    pot     =   1
    
    def theta_teor(W):

        Y   =   11.43/W + 19.88
        return 24./Y*np.sqrt(eps*kappa*np.pi)*n
    
    class alpha():
        
        def __init__(self):
            pass
        
        def set_s(self, s):
            self.s  =   s
        def set_t(self, t):
            self.t  =   t
        def set_pot(self, pot):
            self.pot=   pot
        def set_params(self, s,t, pot):
            self.set_s(s)
            self.set_t(t)
            self.set_pot(pot)

        def get_alpha(self, w):
            #Y       =   11.43/W + 19.88
            #A, B    =   48.7, .0055
            #Length  =   12.8 

            #B   =   4/(Y*np.pi**2)*(A/Length**2 + B*Length**2)
            return self.s/w**self.pot + self.t
            #return self.s/w**self.pot + self.t 
    
    
    alpha   =   alpha()
    def theta_teorFitPar(W, *args):
        
        if len(args) == 1:
            s       =   args[:1]
            alpha.set_s(s)
            
        if len(args) == 2:
            s,t     =   args[:2]
            alpha.set_s(s)
            alpha.set_t(t)
        if len(args) == 3:
            s,t,pot =   args[:3]
            alpha.set_s(s)
            alpha.set_t(t)
            alpha.set_pot(pot)
        
        #Y       =   11.43/W + 19.88
        #A, B    =   48.7, .0055
        #Length  =   12.8 
        #print np.average(4/(Y*np.pi**2)*(A/Length**2 + B*Length**2))
        return magic_par/alpha.get_alpha(W)
        #return 4/(Y*np.pi**2*alpha.get_alpha(W))*(A/Length**2 + B*Length**2)
    
    def overR(W, a1, a2):
        
        return a1/W + a2
    
    def overRandconst(W, a1, a2):
        
        return a1/(W + a2)

    def overR2(W, a1, a2):
        
        return a1/W**2 + a2

    
    for i, data in enumerate(datas):
        energy_table    =   data[6]
        v               =   data[5]
        natoms          =   data[4]
        Wi, Li, Wd, L   =   data[:4]
        phi             =   energy_table[:,3]
        epot            =   (energy_table[:,5] - energy_table[-1,5])/natoms
        
        heights     =   energy_table[:,4]
        inds        =   np.where(thresz < heights)[0]
        indv        =   max(inds)
        indb        =   min(inds)
        
        W           =   Wd #int_toAngst(Wi, edge, key='width', C_C = bond) + 1.2
        print W, Wd, Wd - W  
        
        thetas      =   W/(2*(L/phi))
        
        mus         =   thetas*2./W*360/(2*np.pi)*10 
        
        coll_data[i]=   [W, v, thetas[indv], thetas[indb], mus[indv], mus[indb]]
        
        
        if v > vmax: 
            vmax   = v
        
        
    if plot_mus:
        for i in range(len(coll_data)):
            plt.scatter(coll_data[i,0], coll_data[i,4], 
                        alpha = coll_data[i,1]/vmax, color = 'red')
            plt.scatter(coll_data[i,0], coll_data[i,5], 
                        alpha = coll_data[i,1]/vmax, color = 'blue')
        

        pot = 2.
        alpha.set_pot(pot)
        
        sopm2, topm2    =   curve_fit(theta_teorFitPar, coll_data[:,0], coll_data[:,2], 
                                      p0 = [1., 1])[0][:2]
        
        alpha.set_pot(1)
        alpha.set_t(magic_par)
    
        sopm1       =   curve_fit(theta_teorFitPar, coll_data[:,0], coll_data[:,2], 
                                  p0 = [1.])[0][:1]
        
        print sopm2, topm2
        print sopm1

        plot_widths = np.linspace(np.min(coll_data[:,0]) - .5, np.max(coll_data[:,0]) + 2, 50)
        
        plt.plot(plot_widths, 2/plot_widths*3600/(2*np.pi)*theta_teorFitPar(plot_widths, sopm2, topm2, 2.), '--', 
                 label = r'$\Theta \approx \frac{%.2f}{w^2} + %.3f$' %(sopm2, topm2),
                 color = 'red')
        
        plt.plot(plot_widths, 2/plot_widths*3600/(2*np.pi)*theta_teorFitPar(plot_widths, sopm1, magic_par, 1.), '--', 
                 label = r'$\Theta \approx \frac{%.2f}{w} + %.3f$' %(sopm1, magic_par),
                 color = 'green')
        
        
        plt.scatter(8.46, 4, marker = 'D', color = 'black')
        plt.text(np.min(coll_data[:,0]), 4.1, r'Experim. $\Theta \approx 0.026 = 4deg/nm$')
        plt.legend(frameon = False, loc = 2)
        plt.xlabel('width Angst')
        plt.ylabel(r'$\Theta$')
        plt.twinx()
        
        alpha.set_params(sopm2, topm2, 2)
        plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
                 label = r'$\alpha = \frac{%.3f}{%.2f/w^2 + %.3f} \rightarrow %.3f$' 
                 %(magic_par, sopm2, topm2, magic_par/topm2),
                 color = 'red')
        
        alpha.set_params(sopm1, magic_par, 1)
        plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
                 label = r'$\alpha = \frac{%.3f}{%.2f/w + %.3f} \rightarrow %.3f$' 
                 %(magic_par, sopm1, magic_par, 1),
                 color = 'green')
        
        
        plt.ylabel(r'$\alpha$')
        plt.legend(frameon = False)
        plt.show() 

    
    else:
        for i in range(len(coll_data)):
            plt.scatter(coll_data[i,0], coll_data[i,2], 
                        alpha = coll_data[i,1]/vmax, color = 'red')
            plt.scatter(coll_data[i,0], coll_data[i,3], 
                        alpha = coll_data[i,1]/vmax, color = 'blue')
        
        
        #pot = 1.5
        #alpha.set_pot(pot)
        #sopm15, topm15  =   curve_fit(theta_teorFitPar, coll_data[:,0], coll_data[:,2], p0 = [1., 1])[0][:2]
        
        pot = 2.
        alpha.set_pot(pot)
        sopm2, topm2    =   curve_fit(theta_teorFitPar, coll_data[:,0], coll_data[:,2], 
                                      p0 = [1., 1])[0][:2]
        
        alpha.set_pot(1)
        alpha.set_t(magic_par)
    
        sopm1       =   curve_fit(theta_teorFitPar, coll_data[:,0], coll_data[:,2], 
                                  p0 = [1.])[0][:1]
        
        sopm1 = .16
        #sopm3, topm3, potopm    =   curve_fit(theta_teorFitPar, coll_data[:,0], coll_data[:,2], p0 = [1., 1, 1])[0][:3]
        
        print sopm2, topm2
        print sopm1
        
        
        #plt.plot(coll_data[:,0], np.ones(len(coll_data[:,0]))*.04, '-.', color = 'black')
        #plt.text(np.min(coll_data[:,0]), .041, 'teor')
        
        #plt.scatter(coll_data[:,0], theta_teor(coll_data[:,0]))
        #plt.text(np.min(coll_data[:,0]), .041, 'teor')
    
        plot_widths = np.linspace(np.min(coll_data[:,0]) - .5, np.max(coll_data[:,0]) + 2, 50)
        
        #plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm15, topm15, 1.5), '--', 
        #         label = 'alpha = %.2f/w**%.1f + %.3f' %(sopm15, 1.5, topm15))
        
        plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm2, topm2, 2.), '--', 
                 label = r'$\Theta \approx \frac{%.2f}{w^2} + %.3f$' %(sopm2, topm2),
                 color = 'red')
        
        plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm1, magic_par, 1.), '--', 
                 label = r'$\Theta \approx \frac{%.2f}{w} + %.3f$' %(sopm1, magic_par),
                 color = 'green')
        
        #plt.plot(plot_widths, overR2(plot_widths, sopm2a, topm2a), 
        #         label = r'$\Theta \approx \frac{1}{%.2fw^2} + %.3f$' %(sopm2a, topm2a), )
        
        #plt.plot(plot_widths, theta_teorFitPar(plot_widths, sopm3, topm3, potopm), '.-', 
        #         label = r'$\Theta \approx \frac{%.2f}{w^{%.2f}} + %.3f$' %(sopm3, potopm, topm3),
        #         color = 'red')    
        
        plt.scatter(8.46, .026, marker = 'D', color = 'black')
        plt.text(np.min(coll_data[:,0]), .028, r'Experim. $\Theta \approx 0.026 = 4deg/nm$')
        plt.legend(frameon = False, loc = 2)
        plt.xlabel('width Angst')
        plt.ylabel(r'$\Theta$')
        plt.twinx()
        
        #alpha.set_params(sopm15, topm15, 1.5)
        #plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
        #         label = 'alpha -> %.3f' %alpha.get_alpha(10000))
        alpha.set_params(sopm2, topm2, 2)
        plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
                 label = r'$\alpha = \frac{%.3f}{%.2f/w^2 + %.3f} \rightarrow %.3f$' 
                 %(magic_par, sopm2, topm2, magic_par/topm2),
                 color = 'red')
        
        alpha.set_params(sopm1, magic_par, 1)
        plt.plot(plot_widths, alpha.get_alpha(plot_widths), 
                 label = r'$\alpha = \frac{%.3f}{%.2f/w + %.3f} \rightarrow %.3f$' 
                 %(magic_par, sopm1, magic_par, 1),
                 color = 'green')
        
        
        #alpha.set_params(sopm3, topm3, potopm)
        #plt.plot(plot_widths, alpha.get_alpha(plot_widths), '.-', 
        #         label = r'$\alpha = \frac{%.3f}{%.2f/w^{%.2f} + %.3f} \rightarrow  %.3f$' 
        #         %(magic_par, sopm3, potopm, topm3, magic_par/topm3),
        #         color = 'red')
    
        
        plt.ylabel(r'$\alpha$')
        plt.legend(frameon = False)
        plt.show() 
        
def plot_energy(edge):
    
    Wis     =   [5,7,9,11,13]
    datas   =   get_log_data('LJ', T, '%s_twistTaito' %edge, Wis)
    wili_used   =   []
    n = 0
    colors  =   ['red', 'blue', 'green', 'black', 'cyan', 'yellow']
    aopms   =   []
    teor_aopms  =   []
    W0s     =   []
    
    
    def shear_e(W0, L0, theta):
        #return k*L*(W/R)**2*W/24
        tau =   1.54 #-.13 # -.064/W0 - .13 #- .13
        Y   =   18.75 #19.89 # 11.43/W0 + 19.89 # 19.89 #   
        return 1./6*Y*L0*W0*theta**2*(1. - 2*tau/(k*W0))**2 #- 2*tau**2/(k*W0)*L0
    
    def shear_eNoedge(W0, L0, theta):
        
        return 1./6*k*L0*W0*theta**2

    def x2(thetas, a):
        
        return a*thetas**2
    
    ang_vs  =   np.zeros(len(datas))
    for i, dat in enumerate(datas):
        ang_vs[i] = dat[5]
    
    ang_vs  =   np.max(ang_vs)
    
    
    for data in datas:
        Wi, Li, W, L=   data[:4]
        energy_table=   data[6]
        natoms      =   data[4]
        v           =   data[5]
        phi         =   energy_table[:,3]
        z           =   energy_table[:,4]
        epot        =   (energy_table[:,5] - energy_table[-1,5])
        
        heights     =   energy_table[:,4]
        inds        =   np.where(5 < heights)[0]
        indv        =   max(inds)
        indb        =   min(inds)
        
        thetas      =   W/(2*(L/phi))
        
        if edge == 'ac':
            W0      =   (Wi - 1)*np.sqrt(3)/2*bond
            L0      =   Li*3.*bond - 1.*bond
            Area    =   W0 * L0
            natoms2 =   W0 * L0 / (3*np.sqrt(3)/4*bond**2)
        else: raise
        
        
        
        if [Wi, Li] not in wili_used:
            wili_used.append([Wi, Li])
            weff    =   1 * W0 #- (bond*np.sqrt(3))*.3 # .02/(.16/W0 + .02)*W0 #
            
            epot_t  =   shear_e(weff, L0, thetas)
            plt.plot(thetas, epot_t/Area, '.-', color = colors[n], alpha = .4) #, label = 'Epot teor')
            W0s.append(W0)
            teor_aopms.append(shear_e(weff, L0, 1)/Area)
            n +=    1
        #epot_t2         =   shear_eNoedge(W0, L0, thetas)
        
        
        ind_last = np.where(.02 < thetas)[0][-1]
        
        aopm    =   curve_fit(x2, thetas[ind_last:], epot[ind_last:]/Area, 1)[0][0]
        #aopm    =   curve_fit(x2, thetas[:indb], epot[:indb]/Area, 1)[0][0]
        
        if len(aopms) < n:
            aopms.append([])
        aopms[n - 1].append(aopm)
        
        
        erange  =   np.max(epot/Area)
        
        #plt.plot(thetas, x2(thetas, aopm), color = colors[n -1], alpha = .4) #, label = 'Epot')
        plt.plot(thetas, epot/Area, color = colors[n -1], alpha = .1) #, label = 'Epot')
        
        '''
        plt.plot([thetas[indv], thetas[indv]], [epot[indv]/Area - erange/7, epot[indv]/Area + erange/7], 
                    '--', color = colors[n -1], alpha = .3)
        plt.plot([thetas[indb], thetas[indb]], [epot[indb]/Area - erange/7, epot[indb]/Area + erange/7],
                    '--', color = colors[n -1], alpha = .3)
        '''
        #plt.plot(thetas, z, color =  colors[n -1], alpha = v/ang_vs) #, label = 'maxH')
        #plt.scatter(thetas[indv], heights[indv], color = colors[n -1], alpha = .5)
        #plt.scatter(thetas[indb], heights[indb], color = colors[n -1], alpha = .5)
        
        #plt.plot(thetas, epot_t2, label = 'Epot teor2')
        


    
    plt.xlabel('Theta w/(2R)')
    plt.legend(loc = 1)
    #plt.ylabel('tot Edens eV/angst2')
    #plt.title('width = %.2fAngst, edge = %s' %(W, edge))
    plt.legend(loc = 2)
    plt.ylabel('height max Angst')
    
    plt.show()
    
    
    aopm_av =   np.zeros(len(aopms))
    
    for i in range(n):
        aopm_av[i]  =   np.average(aopms[i])    
        
    plt.scatter(W0s, aopm_av, color = 'red')
    plt.scatter(W0s, teor_aopms, color = 'green', label = 'teor')
    plt.ylabel('aopm')
    plt.xlabel('width')
    plt.legend(loc = 2, frameon = False)
    
    plt.twinx()
    plt.scatter(W0s, np.array(aopm_av)/np.array(teor_aopms), 
                color = 'black', label = 'teor/expr')
    
    plt.ylabel('aopm teor over aopm expr')
    plt.legend(loc = 4, frameon = False)
    plt.show()
    
    
    

def get_data(data, edge):
    
    W,L =   data[0], data[1]
    Dy  =   data[2]
    R   =   data[3]
    Wf  =   int_toAngst(W, edge, key = 'width', C_C = 1.42)
    Lf  =   int_toAngst(L, edge, key = 'length', C_C = 1.42)
    
    theta   =   Wf / (2 * R)  # new journal of physics Bending and buckling -> theta ~ 0.013 (7-ac ribbon).
    print Lf/Wf, Wf, R, theta, L 
    return [Wf, Lf, Dy, R, theta]


def plot_stick():

    #width_f = 7.5
    #print 4./360*2*np.pi*width_f/10/2.
    #exit()
    plot_mus    = False
    edge    =   'ac'
    T, wi   =   10, 11
    LbLt    =   get_stick_data('KC', T, '%s_stickTaito' %edge, wi)
    LbLt.view('i8,i8').sort(order=['f0'], axis=0)
    width_f =   int_toAngst(wi, edge, key='width', C_C = bond)
    
    Rs      =   LbLt[:,0]/(np.pi/3.)
    thetas  =   width_f/(2*Rs) 
    mus     =   thetas*2./width_f*360/(2*np.pi)*10 
    
    #thetas =   width_f/(2*360/mus*1./(2*np.pi)*10)
    def teor(theta, W0):
        #return 2*np.sqrt(3)/9*1000*bond**3/W0*k*(1- 2*tau/(k*W0))**2*theta**3
        return 2*2*np.sqrt(3)/9*1000*bond**3/W0*k*(1- 2*tau/(k*W0))**2*theta**3

    
    if plot_mus:
                
        for i in range(len(mus)):
            if LbLt[i,1] < 3:
                print 'hoos'
                plt.scatter(mus[i], LbLt[i,1]/LbLt[i,0], marker = 'D',color = 'red')
            else:
                plt.scatter(mus[i], LbLt[i,1]/LbLt[i,0], color = 'blue') 
        
        plt.scatter(mus, teor(thetas, width_f), marker = 'v', label = 'teor')
        
        plt.plot([2,2], [0, .4], '--', color = 'black')
        plt.plot([4,4], [0, .4], '--', color = 'black')
        plt.text(1.03, .4, 'experimetl stick below 2deg/angst')
        plt.text(1.03, -.05, 'red squares only one unit cell \n enought to hold bend')
        plt.text(3.7, .4, 'buckling')
        plt.title('Reguired tail length for bend to stick ac, w=7')
        plt.legend(loc = 2, frameon = False)
        plt.xlabel('Deg/angst')
        plt.ylabel('Ltail/Lbend')
        plt.show()

    else:

        for i in range(len(mus)):
            if LbLt[i,1] < 3:
                print 'hoos'
                plt.scatter(thetas[i], LbLt[i,1]/LbLt[i,0], marker = 'D',color = 'red')
                #plt.scatter(thetas[i], LbLt[i,1], marker = 'D',color = 'red')
                
            else:
                plt.scatter(thetas[i], LbLt[i,1]/LbLt[i,0], color = 'blue') 
                #plt.scatter(thetas[i], LbLt[i,1], color = 'blue') 
            print LbLt[i,1]/(np.pi/3)*2*thetas[i]/width_f, LbLt[i,1], 1/(np.pi/3)*2*thetas[i]/width_f, LbLt[i,0]
        #plt.scatter(thetas, teor(thetas, width_f), marker = 'v', label = 'teor')
        #plt.scatter(thetas, teor(thetas, width_f)*LbLt[:,0], marker = 'v', label = 'teor')
        
        plt.twinx()
        plt.scatter(thetas, LbLt[:,0])
        plt.scatter(thetas, LbLt[:,1])
        
        plt.legend(loc = 2, frameon = False)
        plt.xlabel('Theta')
        plt.ylabel('Ltail Angst')
        plt.show()

def plot_stick2(edge):

    
    plot_mus    = True # False #
    #edge    =   'ac'
    T, wis  =   10, [5, 7, 9, 11, 13]
    
    data    =   np.empty(len(wis), dtype = 'object')
    width_f =   np.zeros(len(wis))
    F       =   .006 # eV/angst
    fmax    =   F*1./(3*np.sqrt(3)*bond**2/4)*.4
    a       =   .3
    colors  =   ['blue', 'red', 'black', 'green', 'yellow']
    #print fmax
    
    for i, wi in enumerate(wis):
        data[i]     =   get_stick_data('KC', T, '%s_stickTaito' %edge, wi)
        width_f[i]  =   int_toAngst(wi, edge, key='width', C_C = bond)
    
    print width_f
    
    def axplusb(x, a, b):
        
        return a*x + b
    
    def Mmax(Lt, width_f, f):
        
        def gfun(x):
            return 0
        def hfun(x):
            return width_f/2
        
        def func(y,x):
            return np.sqrt(x**2 + y**2)
            
        return 2*f*dblquad(func, 0, Lt, gfun, hfun)[0]
    
    
    
    #thetas =   width_f/(2*360/mus*1./(2*np.pi)*10)
    def teor(theta, W0):
        #return 2*np.sqrt(3)/9*1000*bond**3/W0*k*(1- 2*tau/(k*W0))**2*theta**3
        return 2*2*np.sqrt(3)/9*1000*bond**3/W0*k*(1- 2*tau/(k*W0))**2*theta**3

    plot2 = []
    plot3 = []
    
    fit_vals    =   [[.0007, .46],
                     [.0008, .4],
                     [.00075, .64],
                     [.00084, .705],
                     [.00085, .868],]

    fit_vals    =   [[.00077, .6],
                     [.00077, .6],
                     [.00077, .6],
                     [.00077, .6],
                     [.00077, .6],]

    
    for j, wi in enumerate(wis):
        LbLt    =   data[j]
        LbLt.view('i8,i8').sort(order=['f0'], axis=0)
        Rs      =   LbLt[:,0]/(np.pi/3.)
        thetas  =   width_f[j]/(2*Rs) 
        mus     =   thetas*2./width_f[j]*360/(2*np.pi)*10 
        print thetas
        
        '''
        if plot_mus:
                    
            for i in range(len(mus)):
                if LbLt[i,1] < 3:
                    print 'hoos'
                    plt.scatter(mus[i], LbLt[i,1]/LbLt[i,0], marker = 'D',color = 'red')
                else:
                    plt.scatter(mus[i], LbLt[i,1]/LbLt[i,0], color = 'blue') 
            
            plt.scatter(mus, teor(thetas, width_f), marker = 'v', label = 'teor')
            
            plt.plot([2,2], [0, .4], '--', color = 'black')
            plt.plot([4,4], [0, .4], '--', color = 'black')
            plt.text(1.03, .4, 'experimetl stick below 2deg/angst')
            plt.text(1.03, -.05, 'red squares only one unit cell \n enought to hold bend')
            plt.text(3.7, .4, 'buckling')
            plt.title('Reguired tail length for bend to stick ac, w=7')
            plt.legend(loc = 2, frameon = False)
            plt.xlabel('Deg/angst')
            plt.ylabel('Ltail/Lbend')
            plt.show()
    
        else:
        '''
        if plot_mus:
            #plt.plot(mus, LbLt[:,1], '-o', color = colors[j], 
            #         alpha = .8, label = r'width = %.2f\AA' %(width_f[j] + 2))
            
            #a,b = curve_fit(axplusb, mus, LbLt[:,1]/LbLt[:,0], p0 = [1.,11])[0][:2]
            
            #plt.plot(mus, LbLt[:,1]/LbLt[:,0], '-o', color = colors[j], 
            #         alpha = .8, label = r'width = %.2f\AA' %(width_f[j] + 2)) 
            #plt.plot(mus, a*mus + b, '--', color = colors[j], 
            #         alpha = .8) 
            #plot2.append([width_f[j], -b/a])
            #plt.plot(mus, LbLt[:,1]/width_f[j], '-o', color = colors[j], 
            #         alpha = .8, label = r'width = %.2f\AA' %(width_f[j] + 2)) 
            R_exp   =   width_f[j]/(2*thetas)
            plt.plot(mus, LbLt[:,1]/R_exp, '-o', color = colors[j], 
                     alpha = .8, label = r'width = %.2f\AA' %(width_f[j] + 2)) 
            
        else:
            plt.plot(thetas, LbLt[:,1], '-o', color = colors[j], alpha = .8)
            
            '''
            a,b = curve_fit(axplusb, thetas, LbLt[:,1]/LbLt[:,0], p0 = [1.,11])[0][:2]
            
            plt.plot(thetas, LbLt[:,1]/LbLt[:,0], '-o', color = colors[j], 
                     alpha = .8, label = r'width = %.2f\AA' %(width_f[j] + 2)) 
            plt.plot(thetas, a*thetas + b, '--', color = colors[j], 
                     alpha = .8) 
            
            
            plot2.append([width_f[j], -b/a])
            '''
        #plt.plot(thetas, LbLt[:,1]/LbLt[:,0], '-o', color = colors[j], alpha = .8) 
        
        print  np.sqrt(width_f[j]/.01*a + a**2)
            
        def gfunc(Lt, theta, f, a):
            L       =   np.sqrt(width_f[j]/theta*a + a**2) + Lt 
            return (1./6.*theta*k*width_f[j]**2 - Mmax(L, width_f[j], f))**2
        
        def Ltail(thetas, f, a):
            
            Lts =   np.zeros(len(thetas))
            print f, a 
            for i, theta in enumerate(thetas):
                Lts[i]  =   fmin(gfunc, 1, args=([theta, f, a]), disp = 0)[0]
            return Lts
        
        def Ltail2(theta, f, a):
            
            return fmin(gfunc, 1, args=([theta, f, a]), disp = 0)[0]**2
        
        Lmax        =   35
        Lts         =   np.linspace(0, Lmax, 15)
        thts        =   np.linspace(0.005, .04, len(Lts))
        teor_Lt     =   np.zeros(len(Lts))
         
        if False:
            fmax, a =   curve_fit(Ltail, thetas, LbLt[:,1], p0 = [fmax, .2])[0]
        else:
            fmax, a =   fit_vals[j]
            
        print fmax, a
        
        for i, Lt in enumerate(Lts):
            #teor_theta[k]   =   6*Mmax(Lt)/(k*width_f**2) 
            teor_Lt[i]      =   Ltail([thts[i]], fmax ,a) #fmin(gfunc, 1, args=([thts[i]]), disp = 0)[0]
            #print thts[i], teor_Lt[i], 6*Mmax(teor_Lt[i], width_f[j])/(k*width_f[j]**2) - thts[i]
        
        
        
        if plot_mus:
            degs        =   thts*2/width_f[j]*3600/(2*np.pi)
            Lbend_teor  =   np.pi/3*width_f[j]/(2*thts)
            #plt.plot(degs, teor_Lt/Lbend_teor, '--', color = colors[j], label = 'width = %i' %wi)
            #plt.plot(degs, teor_Lt/width_f[j], '--', color = colors[j], label = 'width = %i' %wi)
            
            R_teor  =   width_f[j]/(2*thts)
            
            plt.plot(degs, teor_Lt/R_teor, '--', color = colors[j], label = 'width = %i' %wi)
            
            theta0  =   fmin(Ltail2, .001, args=([fmax, a]), disp = 0)[0]
            deg0    =   theta0*2/width_f[j]*3600/(2*np.pi)
            plot3.append([width_f[j], deg0, theta0])
        else:
            plt.plot(thts, teor_Lt, '--', color = colors[j], label = 'width = %i' %wi)
            
        #plt.legend(loc = 2, frameon = False)
    if not plot_mus:
        plt.xlabel('Theta')
        plt.ylabel('Ltail Angst')
    else:
        plt.xlabel('Curve deg/Angst')
        plt.ylabel('Ltail Angst')
     
    plt.legend(loc = 2, frameon = False)
    plt.title('Required tail length for pinning')
    plt.show()
    
    plot3 = np.array(plot3)
    plt.plot(plot3[:,0], plot3[:,1], '-o')
    plt.show()

    plot3 = np.array(plot3)
    plt.plot(plot3[:,0], plot3[:,2], '-o')
    plt.show()

    plot2 = np.array(plot2)
    plt.plot(plot2[:,0], plot2[:,1], '-o')
    plt.show()

bond        =   1.39767578125

#plot_energy('ac')
#plot_kinkOfBend3('ac')
#plot_kinkOfBend_deg('ac')
plot_stick2('ac')

#plot_stick2()
#plot_kinkOfBend3('ac')