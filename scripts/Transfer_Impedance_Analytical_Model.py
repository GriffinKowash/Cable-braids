# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:42:52 2020

@author: kellan.kremer
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ellipk
from scipy.special import ellipe

NORMALIZE = False

######################################################
#Braid Parameters
######################################################

### Yazaki parameters ###

u=1.257E-6
c=16 # number of carriers
n= 5 #number of strands per carrier
d=0.12e-3 # strand diameter (36 gauge = 0.127mm)
a=0.84e-3 #braid radius
Dm = 2*a + 2*d  #mean braid diameter
x = 19.15 # weave angle
sigma = 2.12e7 # strand conductivity 

"""
### Vance Figure 8 K98 parameters ###
u=1.257E-6
c=58 # number of carriers
n=10 # number of strands per carrier
d=0.16e-3 # strand diameter
a=10.0e-3 #braid radius
Dm = 2*a + 2*d  #2*d # mean braid diameter
x = 30 # weave angle
sigma = 6e7 # strand conductivity 
"""
"""
### Sandbox parameters ###
u=1.257E-6
c=42 # number of carriers
n=9 # number of strands per carrier
d=0.01578e-3 # strand diameter
a=1e-3 #braid radius
Dm = 2*a #+ 2*d  #2*d # mean braid diameter
x = 30 # weave angle
sigma = 6e6 # strand conductivity 
"""

freq_low, freq_high, step = int(1e4), int(1.0001e10), int(1e4)


def calc_braid_impedance(model='Vance', **kwargs):  
    f20=(n*d*c) / (2*np.pi*Dm*np.cos(x*(np.pi/180)))
    #f20=(n*d*c) / (4*np.pi*a*np.cos(x*(np.pi/180)))
    k20 =2*f20-f20*f20
    opt20 = 1-k20
    if x < 45:
        e20 = np.sqrt(1-(np.tan(x*(np.pi/180))**2))
        e_sq20 = (1-(np.tan(x*(np.pi/180))**2))
    else:
        e20 = np.sqrt(1 - np.arctan(x * np.pi / 180)**2)
        e_sq20 = (1 - np.arctan(x * np.pi / 180)**2)
    #z = complex(v,j)
    t20 = 9.6 * f20 *(f20**2*(2-f20)**2*d/Dm)**1.5   #check fill
    te20 = 12 * f20 *(f20**2*(2-f20)**2*d/Dm)**1.5   #check fill
    k1 =(np.pi/4)*np.power(((2/3)*f20*np.cos(x*(np.pi/180))+(np.pi/10)),-1) 
    k220 =(np.pi/4)*np.power(((2/3)*f20*np.cos(x*(np.pi/180))+(3/8)),-1) 
    w1 = ((2*np.pi*Dm)/c)*np.cos(x*(np.pi/180))-n*d #hole width
    h1 = (2*d**2)/(w1 + d) #radial spindle separation
    bigk20 = ellipk(e_sq20)
    bige20 = ellipe(e_sq20)
    
    
    ######################################################
    #Vance Model
    ######################################################
    
    
    if model.lower() == "vance": #Hole Inductance
        #m20 = 0.875*((np.pi*u)/(6*c))*np.power(opt20,1.5)*(e_sq20/(bige20-(1-e_sq20)*bigk20))*np.exp(-t20)
        m20 = ((u*2*c)/(np.pi*np.cos(x*(np.pi/180))))*(w1/(np.pi*Dm))**2*np.exp(((-np.pi*d)/w1)-2)
        #m20 = ((np.pi*u)/(6*c))*np.power(opt20,1.5)*(e_sq20/(bige20-(1-e_sq20)*bigk20))
    
        mt20 = m20 # This is the transfer inductance
        #print('The transfer Inductance is:',mt20)
        
        dataX20 = []
        dataf = []
        rt1 = []
        zt20 = []
        rt20v = []  # first value in this list is the transfer resistance
        Lt20v = []
        for i in range(freq_low, freq_high, step):
            X20 = 1J*2*np.pi*i*mt20
            dataX20.append(X20)
            dataf.append(i)
            dlt = np.sqrt(2/(2*np.pi*i*u*sigma))
            r20 = 4/(np.pi*d**2*n*c*sigma*np.cos(x*(np.pi/180)))*(((1+1J)*d/dlt))/np.sinh(((1+1J)*d/dlt))
            zt1 = abs(r20 + X20)
            rt20 = abs(r20)
            Lt20 = abs(X20)
            rt1.append(r20)
            zt20.append(zt1)
            rt20v.append(rt20)
            Lt20v.append(Lt20)
            
        #zt_mag20 =np.sqrt(np.power(rt1,2) + np.power(dataX20,2) + np.power(ptv,2)) 
        #print('The transfer resistance is:', rt20v[0])
        title = '%s Transfer Impedance'% model 
        plt.figure(0)
        if NORMALIZE:
            plt.loglog(dataf,zt20/zt20[0], label = 'Vance Transfer Impedance (normalized)')
        else:
            plt.loglog(dataf,zt20, label = 'Vance Transfer Impedance')
        #plt.loglog(dataf,zt_mag20, label = '38.6 Degrees Simulated')
        #plt.loglog(dataf,zt_sim, label = '38.6 Degrees Simulated')
        plt.title(title,fontsize=16)
        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Transfer Impedance (Ohm/m)")
        plt.xlim(1e4,2e10)  
        plt.ylim(1e-6,1e2)  
        plt.grid(True)
        #plt.savefig(title, dpi=200)
    
    
        title = 'Components of Transfer Impedance'
        plt.figure(1)
        plt.loglog(dataf,rt20v, label = 'Transfer Resistance')
        plt.loglog(dataf,Lt20v, label = 'Transfer Inductance')
        plt.loglog(dataf,zt20, label = 'Simulated Transfer Impedance')
        #plt.loglog(dataf,zt_magsim, label = '38.6 Degrees Simulated')
        #plt.loglog(dataf,zt_sim, label = '38.6 Degrees Simulated')
        plt.title(title,fontsize=16)
        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Transfer Impedance (Ohm/m)")
        plt.xlim(1e4,2e10)  
        plt.ylim(1e-6,1e2)  
        plt.grid(True)
        #plt.savefig(title, dpi=200)
        
        if __name__ != '__main__':
            return dataf, zt20
    
    
        
    ######################################################
    #Tyni Model
    ######################################################
    
    
    elif model.lower() == "tyni": #Hole Inductance + Tyni Braid Inductance
        #m20 = 0.875*((np.pi*u)/(6*c))*np.power(opt20,1.5)*(e_sq20/(bige20-(1-e_sq20)*bigk20))*np.exp(-t20)
    
        m20 = ((u*2*c)/(np.pi*np.cos(x*(np.pi/180))))*(w1/(np.pi*Dm))**2*np.exp(((-np.pi*d)/w1)-2)
        mb20 = -u * (d/(4*np.pi *Dm))*(h1/d)*(1-np.tan(x*(np.pi/180))**2) # Tyni braid inductance
    
        
        #mb20 = u * (d/(4*np.pi *Dm))*(0.22/(f20*np.cos(x*(np.pi/180))))*np.cos(2*k1*x*(np.pi/180)) #Kley braid Inductance
    
        
        mt20 = m20 - mb20 # This is the transfer inductance
    
        
        dataX20 = []
        dataf = []
        rt1 = []
        zt20 = []
        rt20v = [] # first value in this list is the transfer resistance
        Lt20v = []
        for i in range(freq_low, freq_high, step):
            X20 = 1J*2*np.pi*i*mt20
            dataX20.append(X20)
            dataf.append(i)
            dlt = np.sqrt(2/(2*np.pi*i*u*sigma))
            r20 = 4/(np.pi*d**2*n*c*sigma*np.cos(x*(np.pi/180)))*(((1+1J)/dlt)*d)/np.sinh(((1+1J)/dlt)*d)
            zt1 = abs(r20 + X20)
            rt20 = abs(r20)
            Lt20 = abs(X20)
            rt1.append(r20)
            zt20.append(zt1)
            rt20v.append(rt20)
            Lt20v.append(Lt20)
    
        #zt_mag20 =np.sqrt(np.power(rt1,2) + np.power(dataX20,2) + np.power(ptv,2)) 
        
        title = '%s Transfer Impedance'% model 
        plt.figure(0)
        plt.loglog(dataf,zt20, label = 'Tyni Transfer Impedance')
        #plt.loglog(dataf,zt_mag20, label = '38.6 Degrees Simulated')
        #plt.loglog(dataf,zt_sim, label = '38.6 Degrees Simulated')
        plt.title(title,fontsize=16)
        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Transfer Impedance (Ohm/m)")
        plt.xlim(1e4,1e8)  
        plt.ylim(1e-4,1e0)  
        plt.grid(True)
        #plt.savefig(title, dpi=200)
    
    
        title = 'Components of Transfer Impedance'
        plt.figure(1)
        plt.loglog(dataf,rt20v, label = 'Transfer Resistance')
        plt.loglog(dataf,Lt20v, label = 'Transfer Inductance')
        plt.loglog(dataf,zt20, label = 'Simulated Transfer Impedance')
        #plt.loglog(dataf,zt_magsim, label = '38.6 Degrees Simulated')
        #plt.loglog(dataf,zt_sim, label = '38.6 Degrees Simulated')
        plt.title(title,fontsize=16)
        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Transfer Impedance (Ohm/m)")
        plt.xlim(1e4,1e8)  
        plt.ylim(1e-4,1e0)  
        plt.grid(True)
        #plt.savefig(title, dpi=200)
        
        if __name__ != '__main__':
            return dataf, zt20
    
    
    
    ######################################################
    #Demoulin Model
    ######################################################
    
    elif model.lower() == "demoulin": #Hole Inductance +  Tyni Braid Inductance + Porpoising 
        #m20 = 0.875*((np.pi*u)/(6*c))*np.power(opt20,1.5)*(e_sq20/(bige20-(1-e_sq20)*bigk20))*np.exp(-t20)
        #m20 = ((np.pi*u)/(6*c))*np.power(opt20,1.5)*(e_sq20/(bige20-(1-e_sq20)*bigk20))
        m20 = ((u*2*c)/(np.pi*np.cos(x*(np.pi/180))))*(w1/(np.pi*Dm))**2*np.exp(((-np.pi*d)/w1)-2)
        #m20 = ((u*2*c)/(np.pi*np.cos(x*(np.pi/180))))*(w1/(np.pi*Dm))**2*np.exp(((-np.pi*d)/w1)-2)
        #m20 = ((3.1415*u)/(6*c))*np.power(opt20,1.5)*(e_sq20/(bige20-(1-e_sq20)*bigk20))
        #m20 = 0.875*((np.pi*u)/(6*c))*np.power(opt20,1.5)*(e_sq20/(bige20-(1-e_sq20)*bigk20))*np.exp(-t20)
    
        mb20 = -u * (d/(4*np.pi *Dm))*(h1/d)*(1-np.tan(x*(np.pi/180))**2)
        #mb20 = ((u*h1)/(4*np.pi*Dm))*(1-np.tan(x*(np.pi/180))**2) # Tyni braid inductance
    
        k = -1*(1.16/(n*c*d))*np.arctan(n/3)*np.sin((np.pi/2)-(2*(x*(np.pi/180))))*(u/sigma)**(1/2) # when x<45 k is negative @ x=45 k=0
        #k = -2.461e-6
        mt20 = m20 - mb20
        #mt20 = mb20
        #mt20 = m20
        
        dataX20 = []
        dataf = []
        rt1 = []
        ztv = []
        rt20v = []  # first value in this list is the transfer resistance
        Lt20v = []
        ptv = []
        Lpv = []
        kt20v = []
        #zt_mag20 = []
        for i in range(freq_low, freq_high, step):
            X20 = 1J*2*np.pi*i*mt20
            dataX20.append(X20)
            dataf.append(i)
            dlt = np.sqrt(2/(2*np.pi*i*u*sigma))
            r20 = 4/(np.pi*d**2*n*c*sigma*np.cos(x*(np.pi/180)))*(((1+1J)*d/dlt))/np.sinh(((1+1J)*d/dlt))
            p = k*(2*np.pi*i)**(1/2)*np.exp(1J*(np.pi/4))
            rt20 = abs(r20)
            Lt20 = abs(X20)
            pt = abs(p)
            zr = r20
            rt1.append(rt20)
            rt20v.append(rt20)
            Lt20v.append(Lt20)
            ptv.append(pt)
            zt = abs(r20 + X20 + p)
            #zt_test = np.sqrt(np.power(r20,2)+np.power(X20,2)+np.power(p,2))
            Lp = abs(X20 + p)
            Lpv.append(Lp)
            ztv.append(zt)
            kt20 = rt20 + Lp
            kt20v.append(kt20)
            #zt_mag20.append(zt_test)
            
        #zt_mag20 =np.sqrt(np.power(rt1,2) + np.power(dataX20,2) + np.power(ptv,2)) 
        slope = (ztv[-1]-ztv[0])/(dataf[-1]-dataf[0])*(1/(2*np.pi))
        print(slope)
        
        fit = [rt20v[i] + 2*np.pi*slope*dataf[i] for i in range(0,len(dataf))]
        title = '%s Transfer Impedance'% model 
        plt.figure(0)
        plt.loglog(dataf,ztv, label = 'Demoulin Transfer Impedance')
        plt.loglog(dataf,fit, label = 'fit')
        #plt.loglog(dataf,zt_mag20, label = '38.6 Degrees Simulated')
        #plt.loglog(dataf,zt_sim, label = '38.6 Degrees Simulated')
        plt.title(title,fontsize=16)
        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Transfer Impedance (Ohm/m)")
        plt.xlim(1e4,1e8)  
        plt.ylim(1e-4,1e0)  
        plt.grid(True)
        #plt.savefig(title, dpi=200)
    
    
        title = 'Components of Transfer Impedance'
        plt.figure(1)
        plt.loglog(dataf,rt20v, label = 'Transfer Resistance')
        plt.loglog(dataf,Lpv, label = 'Transfer Inductance')
        plt.loglog(dataf,ztv, label = 'Simulated Transfer Impedance')
        #plt.loglog(dataf,zt_magsim, label = '38.6 Degrees Simulated')
        #plt.loglog(dataf,zt_sim, label = '38.6 Degrees Simulated')
        plt.title(title,fontsize=16)
        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Transfer Impedance (Ohm/m)")
        plt.xlim(1e4,1e8)  
        plt.ylim(1e-4,1e0)  
        plt.grid(True)
        #plt.savefig(title, dpi=200)
        
        if __name__ != '__main__':
            return dataf, ztv
        
    print('The transfer resistance is:', rt20v[0])
    print('The transfer Inductance is:',mt20)
    
    



if __name__ == '__main__':
    model = input('Which analytical model would you like to use to compute the transfer impedance for a braided shield? (Vance/Tyni/Demoulin):')

    if model.lower() in ['vance', 'tyni', 'demoulin']:
        calc_braid_impedance(model)
    else:
        print('Enter a valid analytical model.')