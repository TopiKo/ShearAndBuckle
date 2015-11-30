'''
Created on 16.11.2015

@author: tohekorh
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from read import get_datas, get_log_data, get_stick_data
from scipy.optimize import curve_fit
from misc.solvers import int_toAngst
from scipy.optimize import fmin
from scipy.integrate import dblquad
from ase.io.trajectory import PickleTrajectory
edge    =   'ac'
T       =   10
tau     =   1.54    #-.13 # -.064/W0 - .13 #- .13
k       =   18.75   #19.89 # 11.43/W0 + 19.89 # 19.89 #   
bond    =   1.39767578125

rc('text', usetex=True)
rcParams.update({'font.size': 12})

fig_w   =   4

def get_axes(n,m, hPerw):
    
    #fig 2
    #h,w        =   .15, .17       
    # fig 3
    h,w        =   .25, .15       
    
    fig         =   plt.figure(figsize = (fig_w, hPerw*fig_w))
    
    b           =   (1. - h/2*(m - 1) - h + h/3)/m
    #b           =   (1. - h/2*(m - 1) - h)/m
    
    yfrac       =   []
    axs         =   []
    
    for j in range(m):
        k       =   m - j - 1
        yfrac.append(h/2 + k*(b + h/2))
        
        a       =   (1. - 1.3*w)/n
        xfrac   =   []
        xw      =   []
        for i in range(n):
            xfrac.append(w + a*i)
            xw.append(a)
    
        for i in range(len(xfrac)):
            axs.append(fig.add_axes([xfrac[i], yfrac[j], xw[i], b]))   
        
    return fig, axs

def plot_fig2():
    
    
    Wis     =   [5,7,9,11,13] #,9,11,13]
    datas   =   get_log_data('LJ', T, '%s_twistTaito' %edge, Wis)
    wili_used   =   []
    colors  =   ['red', 'blue', 'green', 'black', 'cyan', 'yellow']
    W0s     =   []
    offset  =   .4
    thresz  =   5
    
    coll_data_ac=   np.zeros((len(datas), 7))    
    vmax    =   -1
    fig, axs    =   get_axes(1,3, 1.5)
    
    def shear_eDens(W0, theta):
        return 1./6*k*theta**2*(1. - 2*tau/(k*W0))**2 #- 2*tau**2/(k*W0)*L0
    
    for i, data in enumerate(datas):
        Wi, Li, W, L=   data[:4]
        energy_table=   data[6]
        natoms      =   data[4]
        v           =   data[5]
        phi         =   energy_table[:,3]
        z           =   energy_table[:,4]
        epot        =   (energy_table[:,5] - energy_table[-1,5])*1000
        cutn        =   int(len(epot)/2)
        
        if Wi == 5: m = 0
        if Wi == 7: m = 1
        if Wi == 9: m = 2
        if Wi == 11: m = 3
        if Wi == 13: m = 4

        epot        =   epot[cutn:]
        
        if edge == 'ac':
            W0      =   (Wi - 1)*np.sqrt(3)/2*bond
            L0      =   Li*3.*bond - 1.*bond
            Area    =   W0 * L0
        else: raise
        
        
        #print W0
        
        thetash     =   W/(2*(L/phi))
        thetas      =   thetash[cutn:]
        
        
        heights     =   energy_table[:,4]
        inds        =   np.where(thresz < heights)[0]
        indv        =   max(inds)
        indb        =   min(inds)
        
        W           =   int_toAngst(Wi, 'ac', key='width', C_C = bond)
        
        mus         =   thetash*2./W*360/(2*np.pi)*10 
        coll_data_ac[i]    =   [W, v, thetash[indv], thetash[indb], mus[indv], mus[indb], Wi]
        
        if v > vmax: 
            vmax    =   v
        
        
        if [Wi, Li] not in wili_used:
            wili_used.append([Wi, Li])
            
            epot_tDens      =   shear_eDens(W0, thetas)*1000
            axs[0].plot(thetas, epot_tDens + m*offset, '-', color = colors[m], alpha = 1) #, label = 'Epot teor')
            W0s.append(W0)
            axs[0].text(.01 - (m+1)*.0017, (m+1)*offset*.85 + offset/10, r'N=%i' %Wi, 
                        color = colors[m])
            axs[1].text(.08 - float(m)/4*.04, float(m)/4*6. + 9, r'N=%i' %Wi, 
                        color = colors[m])
            
        axs[0].plot(thetas, epot/Area + m*offset, color = colors[m], alpha = .05) #, label = 'Epot')
        axs[1].plot(thetash, heights, color = colors[m], alpha = .15) #, label = 'Epot')
        


    axs[0].set_xlim([0, .03])
    axs[0].set_ylim([-.2, 3.5])
    axs[0].set_xlabel(r'$\Theta$')
    axs[0].set_ylabel(r'Energy density (m$eV/$\AA$^2$)')
    
    axs[1].set_ylim([2, 17])
    axs[1].set_xlabel(r'$\Theta$')
    axs[1].set_ylabel(r'Max Height (\AA)')
    
    axs[0].legend(loc = 2)
    
    # KINK
    magic_par   =   0.0228171
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
            return magic_par/(self.s/w**self.pot + self.t)
            
    
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
        
        return magic_par/alpha.get_alpha(W)
    
    ac_curvature    =   np.zeros(5)
    amount          =   np.zeros(5)
    for i in range(len(coll_data_ac)):
        Wi  =   coll_data_ac[i, -1]
        
        if Wi == 5: 
            m = 0
            ac_curvature[0]   +=  coll_data_ac[i,4]
            amount[0]  +=  1
            
        if Wi == 7: 
            m   = 1
            ac_curvature[1]   +=  coll_data_ac[i,4]
            amount[1]  +=  1
            
        if Wi == 9: 
            m = 2
            ac_curvature[2]   +=  coll_data_ac[i,4]
            amount[2]  +=  1
            
        if Wi == 11: 
            m = 3
            ac_curvature[3]   +=  coll_data_ac[i,4]
            amount[3]  +=  1
            
        if Wi == 13: 
            m = 4
            ac_curvature[4]   +=  coll_data_ac[i,4]
            amount[4]  +=  1
            
        
        axs[2].scatter(coll_data_ac[i,0], coll_data_ac[i,4], 
                    alpha = coll_data_ac[i,1]/vmax, color = colors[m])
        axs[2].scatter(coll_data_ac[i,0], coll_data_ac[i,5], 
                    alpha = coll_data_ac[i,1]/vmax, marker = 'D', 
                    color = colors[m])
    
    for i in range(5):    
        print 'average curvature %i-ac = %.4f deg/nm' %(Wis[i], ac_curvature[i]/amount[i])
    

    alpha.set_pot(1)
    alpha.set_t(magic_par)

    sopm1_ac       =   curve_fit(theta_teorFitPar, coll_data_ac[:,0], coll_data_ac[:,2], 
                              p0 = [1.])[0][:1]
    
    plot_widths = np.linspace(np.min(coll_data_ac[:,0]) - .5, np.max(coll_data_ac[:,0]) + 2, 50)
    
    axs[2].plot(plot_widths, 2/plot_widths*3600/(2*np.pi)*theta_teorFitPar(plot_widths, sopm1_ac, magic_par, 1.), 
                '--', label = r'$\Theta \approx \frac{%.2f}{w} + %.3f$' %(sopm1_ac, magic_par),
                color = 'black')
    
    #axs[2].legend(loc = 1, frameon = False)
    axs[2].set_ylabel(r'Curvature $\kappa$ (deg/nm)')
    axs[2].set_xlabel(r'Width (\AA)')
    
    #xticks  =   (np.array(Wis) - 1)*np.sqrt(3)/2*1.4 
    #xticklabels =   Wis
    #axs[2].set_xticks(xticks)
    #axs[2].set_xticklabels(xticklabels)    
    plt.show()

def plot_fig3():
    
    plot_mus    = True # False #
    edge    =   'ac'
    T, wis  =   10, [5, 7, 9, 11, 13]
    _, axs    =   get_axes(1,1, .8)
    
    
    data    =   np.empty(len(wis), dtype = 'object')
    width_f =   np.zeros(len(wis))
    F       =   .006 # eV/angst
    fmax    =   F*1./(3*np.sqrt(3)*bond**2/4)*.4
    a       =   .3
    colors  =   ['red', 'blue', 'green', 'black', 'cyan', 'yellow']
    
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
        
        if wi == 5: m = 0
        if wi == 7: m = 1
        if wi == 9: m = 2
        if wi == 11: m = 3
        if wi == 13: m = 4
        
        
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
            axs[0].scatter(mus, LbLt[:,1]/R_exp, color = colors[m])
            axs[0].text(8 - m*1.5, .55 - m*.04, r'w=%i' %wi, color = colors[m]) 
            
        else:
            axs[0].plot(thetas, LbLt[:,1], '-o', color = colors[m], alpha = .8)
            
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
            
            axs[0].plot(degs, teor_Lt/R_teor, '--', color = colors[m], label = 'width = %i' %wi)
            
            theta0  =   fmin(Ltail2, .001, args=([fmax, a]), disp = 0)[0]
            deg0    =   theta0*2/width_f[j]*3600/(2*np.pi)
            plot3.append([width_f[j], deg0, theta0])
            
        else:
            axs[0].plot(thts, teor_Lt, '--', color = colors[m], label = 'width = %i' %wi)
            
        #plt.legend(loc = 2, frameon = False)
    if not plot_mus:
        axs[0].set_xlabel(r'$\Theta$')
        axs[0].set_ylabel(r'$L_t/R$')
    else:
        axs[0].set_xlabel(r'Curvature $\kappa$ (deg/nm)')
        axs[0].set_ylabel(r'$L_t/R$')
     
    #axs[0].legend(loc = 2, frameon = False)
    #axs[0].set_title('Required tail length for pinning')
    axs[0].set_xlim([0,12.5])
    
    
    a = plt.axes([.67, .27, .27, .25])
    plot3 = np.array(plot3)
    plt.plot(plot3[:,0], plot3[:,1], '--', color = 'black', alpha = .5)
    for i in range(5):
        plt.scatter(plot3[i,0], plot3[i,1], color = colors[i])
    plt.yticks(np.linspace(.6, 2.2, 5))
    plt.xticks(np.linspace(4, 16, 4))
    plt.xlabel(r'Width (\AA)')
    plt.ylabel(r'$\kappa$ (deg/nm)')
    
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

def plot_fig3a():
    
    plot_mus    = True # False #
    edge    =   'ac'
    T, wis  =   10, [5, 7, 9, 11, 13]
    _, axs    =   get_axes(1,1, .8)
    
    
    data    =   np.empty(len(wis), dtype = 'object')
    width_f =   np.zeros(len(wis))
    F       =   .006 # eV/angst
    fmax    =   F*1./(3*np.sqrt(3)*bond**2/4)*.4
    a       =   .3
    colors  =   ['red', 'blue', 'green', 'black', 'cyan', 'yellow']
    
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
        
        if wi == 5: m = 0
        if wi == 7: m = 1
        if wi == 9: m = 2
        if wi == 11: m = 3
        if wi == 13: m = 4
        
        
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
            axs[0].scatter(LbLt[:,1]/R_exp , mus, color = colors[m])
            axs[0].text(.55 - m*.04, 9 - m*1.5, r'w=%i' %wi, color = colors[m]) 
            
        else:
            axs[0].plot(thetas, LbLt[:,1], '-o', color = colors[m], alpha = .8)
            
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
        
        Lmax        =   55
        Lts         =   np.linspace(0, Lmax, 15)
        thts        =   np.linspace(0.005, .045, len(Lts))
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
            
            axs[0].plot(teor_Lt/R_teor, degs, '--', color = colors[m], label = 'width = %i' %wi)
            
            theta0  =   fmin(Ltail2, .001, args=([fmax, a]), disp = 0)[0]
            deg0    =   theta0*2/width_f[j]*3600/(2*np.pi)
            plot3.append([width_f[j], deg0, theta0])
            
        else:
            axs[0].plot(teor_Lt, thts, '--', color = colors[m], label = 'width = %i' %wi)
            
        #plt.legend(loc = 2, frameon = False)
    if not plot_mus:
        axs[0].set_ylabel(r'$\Theta$')
        axs[0].set_xlabel(r'$L_t/R$')
    else:
        axs[0].set_ylabel(r'Curvature $\kappa$ (deg/nm)')
        axs[0].set_xlabel(r'$L_t/R$')
     
    #axs[0].legend(loc = 2, frameon = False)
    #axs[0].set_title('Required tail length for pinning')
    axs[0].set_ylim([0,13])
    axs[0].set_xticks(np.linspace(0, .7, 8))
    
    
    a = plt.axes([.3, .675, .27, .23])
    plot3 = np.array(plot3)
    plt.plot(plot3[:,0], plot3[:,1], '--', color = 'black', alpha = .5)
    for i in range(5):
        plt.scatter(plot3[i,0], plot3[i,1], color = colors[i])
    plt.yticks(np.linspace(.6, 2.2, 5))
    plt.xticks(np.linspace(4, 16, 4))
    plt.xlabel(r'Width (\AA)')
    plt.ylabel(r'$\kappa$ (deg/nm)')
    
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

def plot_fig3b():
    
    plot_mus    = True # False #
    edge    =   'ac'
    T, wis  =   10, [5, 7, 9, 11, 13]
    _, axs    =   get_axes(1,1, .8)
    
    
    data    =   np.empty(len(wis), dtype = 'object')
    width_f =   np.zeros(len(wis))
    F       =   .006 # eV/angst
    fmax    =   F*1./(3*np.sqrt(3)*bond**2/4)*.4
    a       =   .3
    colors  =   ['red', 'blue', 'green', 'black', 'cyan', 'yellow']
    
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
        
        if wi == 5: m = 0
        if wi == 7: m = 1
        if wi == 9: m = 2
        if wi == 11: m = 3
        if wi == 13: m = 4
        
        
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
            axs[0].scatter(LbLt[:,1]/10, mus, color = colors[m])
            axs[0].text(3 + 1.5/4*m, 11 - 7./4*m, r'N=%i' %wi, color = colors[m]) 
            
        else:
            axs[0].plot(thetas, LbLt[:,1], '-o', color = colors[m], alpha = .8)
            
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
        
        Lmax        =   55
        Lts         =   np.linspace(0, Lmax, 15)
        thts        =   np.linspace(0.005, .045, len(Lts))
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
            
            axs[0].plot(teor_Lt/10, degs, '--', color = colors[m], label = 'width = %i' %wi)
            
            theta0  =   fmin(Ltail2, .001, args=([fmax, a]), disp = 0)[0]
            deg0    =   theta0*2/width_f[j]*3600/(2*np.pi)
            plot3.append([width_f[j], deg0, theta0])
            
        else:
            axs[0].plot(teor_Lt/10, thts, '--', color = colors[m], label = 'width = %i' %wi)
            
        #plt.legend(loc = 2, frameon = False)
    if not plot_mus:
        axs[0].set_ylabel(r'$\Theta$')
        axs[0].set_xlabel(r'$L_t (nm)$')
    else:
        axs[0].set_ylabel(r'Curvature $\kappa$ (deg/nm)')
        axs[0].set_xlabel(r'$L_t$ (nm)')
     
    #axs[0].legend(loc = 2, frameon = False)
    #axs[0].set_title('Required tail length for pinning')
    axs[0].set_ylim([0,13])
    axs[0].set_xlim([-.1,5])
    axs[0].set_xticks(np.linspace(0, 5, 6))
    
    
    a = plt.axes([.3, .675, .27, .23])
    plot3 = np.array(plot3)
    plt.plot(plot3[:,0], plot3[:,1], '--', color = 'black', alpha = .5)
    for i in range(5):
        plt.scatter(plot3[i,0], plot3[i,1], color = colors[i], marker = 'D')
    plt.yticks(np.linspace(.6, 2.2, 5))
    plt.xticks(np.linspace(4, 16, 4))
    plt.xlabel(r'Width (\AA)')
    plt.ylabel(r'$\kappa$ (deg/nm)')
    
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
    
    
def plot_fig3c():
    
    plot_mus    = True # False #
    edge    =   'ac'
    T, wis  =   10, [5, 7, 9, 11, 13]
    _, axs    =   get_axes(1,1, .8)
    
    deva    =   .6
    data    =   np.empty(len(wis), dtype = 'object')
    width_f =   np.zeros(len(wis))
    F       =   .006 # eV/angst
    fmax    =   F*1./(3*np.sqrt(3)*bond**2/4)*.4
    a       =   .3
    colors  =   ['red', 'blue', 'green', 'black', 'cyan', 'yellow']
    
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
        
        if wi == 5: m = 0
        if wi == 7: m = 1
        if wi == 9: m = 2
        if wi == 11: m = 3
        if wi == 13: m = 4
        
        
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
            
            add_L   =   np.sqrt(2*R_exp*deva + deva**2)/10.
            
            axs[0].scatter(LbLt[:,1]/10 + add_L, mus, color = colors[m])
            axs[0].text(3 + 1.5/4*m, 11 - 7./4*m, r'w=%i' %wi, color = colors[m]) 
            
        else:
            axs[0].plot(thetas, LbLt[:,1], '-o', color = colors[m], alpha = .8)
            
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
        
        Lmax        =   55
        Lts         =   np.linspace(0, Lmax, 15)
        thts        =   np.linspace(0.005, .045, len(Lts))
        teor_Lt     =   np.zeros(len(Lts))
         
        if False:
            fmax, a =   curve_fit(Ltail, thetas, LbLt[:,1], p0 = [fmax, .2])[0]
        else:
            fmax, a =   fit_vals[j]
            
        print fmax, a
        
        for i, Lt in enumerate(Lts):
            #teor_theta[k]   =   6*Mmax(Lt)/(k*width_f**2) 
            teor_Lt[i]      =   Ltail([thts[i]], fmax ,a) + np.sqrt(width_f[j]/thts[i]*deva + deva**2) #fmin(gfunc, 1, args=([thts[i]]), disp = 0)[0]
            #print thts[i], teor_Lt[i], 6*Mmax(teor_Lt[i], width_f[j])/(k*width_f[j]**2) - thts[i]
        
        
        
        if plot_mus:
            degs        =   thts*2/width_f[j]*3600/(2*np.pi)
            Lbend_teor  =   np.pi/3*width_f[j]/(2*thts)
            #plt.plot(degs, teor_Lt/Lbend_teor, '--', color = colors[j], label = 'width = %i' %wi)
            #plt.plot(degs, teor_Lt/width_f[j], '--', color = colors[j], label = 'width = %i' %wi)
            
            R_teor  =   width_f[j]/(2*thts)
            
            axs[0].plot(teor_Lt/10, degs, '--', color = colors[m], label = 'width = %i' %wi)
            
            theta0  =   fmin(Ltail2, .001, args=([fmax, a]), disp = 0)[0]
            deg0    =   theta0*2/width_f[j]*3600/(2*np.pi)
            plot3.append([width_f[j], deg0, theta0])
            
        else:
            axs[0].plot(teor_Lt/10, thts, '--', color = colors[m], label = 'width = %i' %wi)
            
        #plt.legend(loc = 2, frameon = False)
    if not plot_mus:
        axs[0].set_ylabel(r'$\Theta$')
        axs[0].set_xlabel(r'$L (nm)$')
    else:
        axs[0].set_ylabel(r'Curvature $\kappa$ (deg/nm)')
        axs[0].set_xlabel(r'$L$ (nm)')
     
    #axs[0].legend(loc = 2, frameon = False)
    #axs[0].set_title('Required tail length for pinning')
    axs[0].set_ylim([0,13])
    axs[0].set_xlim([1,6])
    axs[0].set_xticks(np.linspace(0, 7, 8))
    
    
    a = plt.axes([.3, .675, .27, .23])
    plot3 = np.array(plot3)
    plt.plot(plot3[:,0], plot3[:,1], '--', color = 'black', alpha = .5)
    for i in range(5):
        plt.scatter(plot3[i,0], plot3[i,1], color = colors[i])
    plt.yticks(np.linspace(.6, 2.2, 5))
    plt.xticks(np.linspace(4, 16, 4))
    plt.xlabel(r'Width (\AA)')
    plt.ylabel(r'$\kappa$ (deg/nm)')
    
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



def plot_fig1():
    
    mdfile  =   '/space/tohekorh/ShearSlide/files/LJ_10/ac_twistTaito/w=7/md_L=24_v=0.20_00.traj'
    #mdfile  =   '/space/tohekorh/ShearSlide/files/KC_10/ac_stickTaito/w=7/r=20/md_L=145_stL=27_00.traj'
    
    traj    =   PickleTrajectory(mdfile)
    
    z_init  =   np.average(traj[0].positions[:,2])
    xrange  =   [np.min(traj[0].positions[:,0]), np.max(traj[0].positions[:,0])]
    yrange  =   [np.min(traj[0].positions[:,1]), np.max(traj[0].positions[:,1])]
    
    from ase.structure import graphene_nanoribbon
    from ase.visualize import view
    from ase import Atoms
    base    =   graphene_nanoribbon(70, 50, type='armchair', saturated=False,
                                    C_C=bond, main_element='N')
    
    base.rotate([1,0,0], np.pi/2)
    base.positions[:,2]     =   z_init - 3.8 
    atoms_v =   Atoms()
    n       =   int(len(traj)/1.3)
    #n       =   
    
    nsnap   =   3
    
    atoms_use   =   [traj[i] for i in np.array(range(nsnap))*n/nsnap]
    
    for i, atoms in enumerate(atoms_use):
        
        atoms.positions[:,1] = -atoms.positions[:,1] - i*20 
        #atoms.positions[:,0] = -atoms.positions[:,0] 
        
        atoms_v +=  atoms
    
    cent_atoms  =   np.array([np.average(atoms_v.positions[:,0]), np.average(atoms_v.positions[:,1]), 0])
    cent_base   =   np.array([np.average(base.positions[:,0]), np.average(base.positions[:,1]), 0])
    base.translate(cent_atoms - cent_base)
        
    #atoms_v +=  base
    view(atoms_v, viewer = 'vmd')
    
    
    
#plot_fig1()
plot_fig2()    
#plot_fig3b()
#plot_fig3a()