'''
Created on 14.9.2015

@author: tohekorh
'''
import numpy as np
import matplotlib.pyplot as plt
evdw    =   .017
kappa   =   .95
k       =   21
sigma   =   3.4
eps     =   0.002843732471143

ws      =   [1] #,4,5,6,7,8,9,10,11,12,13,14,15]
n       =   4./(3*np.sqrt(3)*1.42**2)

print (16*np.sqrt(eps)*np.sqrt(kappa)*n*np.sqrt(3*np.pi))/k
print 24./k*np.sqrt(eps*np.pi*kappa)*n


#print (9*eps*L**4*n**2 + kappa*np.pi**3)/(2*k*L**2*np.pi)
def theta(Lin):
    
    return (40*(3*eps*Lin**4*n**2 + kappa*np.pi**3)*W)/(k*Lin**2*np.pi*(3*W + 2*0))
    
    # Peka
    #return 40./(3*k*np.pi**2)*(15*evdw*Lin**2/sigma**2 + kappa*np.pi**4/Lin**2)

def theta2(Lin, W, Win):
    
    #return (40*(3*eps*Lin**4*n**2 + kappa*np.pi**3)*W)/(k*Lin**2*np.pi*(3*W + 2*Win))
    
    #return (9*eps*Lin**4*n**2 + kappa*np.pi**3)/(2*k*Lin**2*np.pi)
    
    
    #return (40*(3*eps*Lin**4*n**2 + kappa*np.pi**3)*W*(W - Win)**4)/(k*Lin**2*np.pi*(5*W - 2*Win)*Win**4)
    #return 1./k * (kappa*np.pi**2/Lin**2 + 9*eps*Lin**2*n**2/np.pi)    
    
    #return 2*W/k*1./(W- Win)*(kappa*np.pi**2/Lin**2+ 9*eps*Lin**2*n**2/np.pi)
    
    return 2/(k*Lin**2*np.pi*(1-Win/W))*(9*eps*Lin**4*n**2 + kappa*np.pi**3)
    #return 20*W/(k*Lin**2*np.pi*(5*W-2*Win))*(3*eps*Lin**4*n**2 + kappa*np.pi**3)
    
    #return (40*(3*eps*Lin**4*n**2 + kappa*np.pi**3)*(W - 
    #Win)**4)/(k*Lin**2 *np.pi* W**2*(3*W**2 - 10*W*Win + 10*Win**2))
    
    
def wder(Lin, theta):
    
    A   =   15*evdw/sigma**2/6
    B   =   -3*kappa*np.pi**4 
    C   =   .1*np.pi**2*k*theta

    return A*Lin**4 + C*Lin**2 + B 

W       =   1.

Wins    =   np.linspace(0,W,10, endpoint = False)

#Wins    =   [.99]
for Wi in Wins:
    for w in ws:
        Lins    =   np.linspace(w/5, 24*w, 100)
        thetas  =   np.zeros(len(Lins))
        thetas2 =   np.zeros(len(Lins))
        wders   =   np.zeros(len(Lins))

        for i, Lin in enumerate(Lins):
            thetas2[i]  =   theta2(Lin, W, Wi)
            thetas[i]   =   theta(Lin)
            
            #wders[i]    =   wder(Lin, thetas[i])
        plt.plot(Lins, thetas2, '--')
    plt.plot(Lins, thetas)
        
    #plt.plot(Lins, wders)


plt.xlabel('Lins')

plt.show()
'''

Win =   np.linspace(0, W, 100)
A   =   1
theta_val   =   .015

Lins    =   np.linspace(0, 20, 10)
for Lin in Lins:
    Es = lambda Win: -k*A**2*np.pi**2/(80*Lin*W)*theta_val*(5*W - 2*Win)*Win 
    Eb = lambda Win: kappa*A**2*np.pi**4*Win/(4*Lin**3) 
    Evdw = lambda Win: A**2*3./4*Lin*np.pi*eps*n**2*Win
    
    Es2 = lambda Win: -2*A**2*theta_val*k*np.pi**2/4/Lin*Win*(1-Win/W) 
    Eb2 = lambda Win: kappa*A**2*np.pi**4*Win/Lin**3 
    Evdw2 = lambda Win: 9*A**2*eps*Lin*n**2*np.pi*Win
    
    
    #plt.plot(Win, Es2(Win), label = 'Es')
    #plt.plot(Win, Eb2(Win), label = 'Eb')
    #plt.plot(Win, Evdw2(Win), label = 'Evdw')
    plt.plot(Win, Eb2(Win) + Evdw2(Win) + Es2(Win), label = 'Etot')
    
plt.show()
'''