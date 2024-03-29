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

braid = 'Vance k80'

if braid == 'Steve':
    ### Steve parameters ###
    u=1.257E-6
    c=48 # number of carriers
    n= 18 #number of strands per carrier
    d=0.127e-3 # strand diameter
    a=9.92e-3 #braid radius
    Dm = 2*a + 2*d  #mean braid diameter
    x = 30 # weave angle
    sigma = 1.45e7 # strand conductivity vance

elif braid == 'Yazaki':
    ### Yazaki parameters ###
    u=1.257E-6
    c=16 # number of carriers
    n= 5 #number of strands per carrier
    d=0.12e-3 # strand diameter
    a=0.82e-3 #braid radius
    Dm = 2*a + 2*d  #mean braid diameter
    x = 19.15 # weave angle
    sigma = 6e7 # strand conductivity 

elif braid == 'Vance k80':
    ### Vance K80 parameters ###
    u=1.257E-6
    c=42 # number of carriers
    n=9 # number of strands per carrier
    d=0.16e-3 # strand diameter
    a=9.92e-3 #braid radius
    Dm = 2*a + 2*d # mean braid diameter
    x = 30 # weave angle
    sigma = 6e7 # strand conductivity 

elif braid == 'Vance k69':
    ### Griffin k69 parameters ###
    u=1.257E-6
    c=24 # number of carriers
    n=12 # number of strands per carrier
    d=0.16e-3 # strand diameter
    a=8.00e-3 #braid radius
    Dm = 2*a + 2*d # mean braid diameter
    x = 20 # weave angle
    sigma = 6e7 # strand conductivity 

elif braid == 'Yazaki':
    ### Schippers 12 parameters ###
    u=1.257E-6
    c=32 # number of carriers
    n=5 # number of strands per carrier
    d=0.202e-3 # strand diameter
    a=4.00e-3 #braid radius
    Dm = 2*a + 2*d # mean braid diameter
    x = 35.2 # weave angle
    sigma = 6e7 # strand conductivity 

elif braid == 'Schippers 8':
    ### Schippers 8 parameters ###
    u=1.257E-6
    c=24 # number of carriers
    n=7 # number of strands per carrier
    d=0.160e-3 # strand diameter
    a=3.00e-3 #braid radius
    Dm = 2*a + 2*d # mean braid diameter
    x = 38.6 # weave angle
    sigma = 6e7 # strand conductivity 

elif braid == 'sandbox':
    ### Sandbox parameters ###
    u=1.257E-6
    c=32 # number of carriers
    n=5 # number of strands per carrier
    d=0.202e-3 # strand diameter
    a=4e-3 #braid radius
    Dm = 2*a + 2*d  #2*d # mean braid diameter
    x = 35.2 # weave angle
    sigma = 5e7 # strand conductivity 

elif braid == 'Otin 1':
    # Sample 1 from Otin & Isanta 2013
    u=1.257E-6
    c=16 # number of carriers
    n=3 # number of strands per carrier
    d=0.100e-3 # strand diameter
    a=0.750e-3 #braid radius
    Dm = 2*a + 2*d  # mean braid diameter
    x = 21.44 # weave angle
    sigma = 4.12e7 # strand conductivity 

elif braid == 'Otin 2':
    # Sample 2 from Otin & Isanta 2013
    u=1.257E-6
    c=16 # number of carriers
    n=6 # number of strands per carrier
    d=0.101e-3 # strand diameter
    a=1.475e-3 #braid radius
    Dm = 2*a + 2*d  # mean braid diameter
    x = 30.08 # weave angle
    sigma = 4.12e7 # strand conductivity 

elif braid == 'Glenair 100-002A1000':
    # Glenair braid 100-002A1000
    u=1.257E-6
    c=64 # number of carriers
    n=9 # number of strands per carrier
    d=0.127e-3 # strand diameter
    a=12.7e-3 #braid radius
    Dm = 2*a + 2*d  # mean braid diameter
    x = 47.9 # weave angle
    sigma = 6e7 # strand conductivity 
    
elif braid == 'Glenair 100-002A125':
    # Glenair braid 100-002A1000
    u=1.257E-6
    c=24 # number of carriers
    n=5 # number of strands per carrier
    d=0.127e-3 # strand diameter
    a=1.6e-3 / 2 #braid radius
    Dm = 2*a + 2*d  # mean braid diameter
    x = 1 # weave angle
    sigma = 4.12e7 # strand conductivity 
    
elif braid == 'Glenair 100-002A250':
    # Glenair braid 100-002A1000
    u=1.257E-6
    c=24 # number of carriers
    n=16 # number of strands per carrier
    d=0.127e-3 # strand diameter
    a=3.2e-3 #braid radius
    Dm = 2*a + 2*d  # mean braid diameter
    x = 27.6 # weave angle
    sigma = 4.12e7 # strand conductivity 


freq_low, freq_high, step = int(1e4), int(1.0001e10), int(1e4)


def calc_braid_impedance(model='Vance', corrections=False, **kwargs): 
    Dm = 2*a + 2*d # mean braid diameter
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
        
        if False:
            w = 2*np.pi*Dm / c - n*d / np.cos(x*np.pi/180)
            l = w / np.tan(x*np.pi/180)
            if x < 45:
                e = np.sqrt(1 - (w / l)**2)
            else:
                e = np.sqrt(1 - (l / w)**2)
            F = n * c * d / (2 * np.pi * Dm * np.cos(x*np.pi/180))
            K = 2*F - F**2
            m20 = np.pi*u/(6*c) * (1-K)**(3/2) * (e**2 / (ellipe(e) - (1 - e**2) * ellipk(e)))
        
    
        mt20 = m20 # This is the transfer inductance
        print('The transfer Inductance is:',mt20)
        
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
        
        return rt20v, Lt20v, zt20
        
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
    
        
        if corrections:
            mt20 = m20 + mb20 # This is the transfer inductance
        else:
            mt20 = m20 - mb20
        
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
        Dm = 2*a + 2*d
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
        if corrections:
            mt20 = m20 + mb20
        else:
            mt20 = m20 + mb20
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
        

    elif model.lower() == 'kley':

        Dm = 2*a + 2.5*d
        G0 = c*n*d / (2*np.pi*Dm)  # minimal filling factor
        G = G0 / np.cos(x*np.pi/180)  # filling factor
        B = G*(2-G)  # optical coverage
        Rgs = 4 / (sigma*c*n*(d**2)*np.pi*np.cos(x*np.pi/180))  # DC resistance
        
        k1 = (np.pi/4) / ((2/3)*G0 + np.pi/10)
        k2 = (np.pi/4) / ((2/3)*G0 + 3/8)
        te = 12*G*np.cbrt((B**2)*d/Dm)
        th = 9.6*G*np.cbrt((B**2)*d/Dm)
        dR = 0.67*d/np.sqrt(np.cos(x*np.pi/180))
        
        dataf = []
        data_wLs = []
        data_Lt = []
        data_Zr = []
        data_Zt = []
        
        dataf = np.arange(freq_low, freq_high, step)
            
        delta = np.sqrt(2/(2*np.pi*dataf*u*sigma))
        wLs = (1/(np.pi*sigma*delta*Dm)) * (10*np.pi*(G0**2)*np.cos(x*np.pi/180)*(1-G)*np.exp(-te) - (3.3/(2*np.pi*G0))*np.cos(2*k2*x*np.pi/180))
        Lt = (u/c) * (0.875*np.pi/6 * (2 - np.cos(x*np.pi/180)*((1-G)**3)*np.exp(-th)) - 0.11/n * np.cos(2*k1*x*np.pi/180))
        Zr = (2*dR*(1+1j)/delta) / (sigma*G0*np.cos(x*np.pi/180)*(np.pi**2)*Dm*d*np.sinh(dR*(1+1j)/delta))
        Zt = Zr + 1j*2*np.pi*dataf*Lt + (1+1j)*wLs
        
        #ataf.append(f)
        #data_wLs.append(wLs)
        #data_Lt.append(Lt)
        #data_Zr.append(np.abs(Zr))
        #data_Zt.append(np.abs(Zt))

        title = '%s Transfer Impedance'% model 
        plt.figure(0)
        if NORMALIZE:
            plt.loglog(dataf, np.abs(Zt)/np.abs(data_Zt[0]), label = 'Kley Transfer Impedance (normalized)')
        else:
            plt.loglog(dataf, np.abs(Zt), label = 'Kley Transfer Impedance')
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
    
        
        if __name__ != '__main__':
            return dataf, np.abs(Zt)
        
        rt20v = np.abs(Zr)
        mt20 = np.abs(wLs + Lt)
        
    elif model.lower() == 'beatrics':
        Dm = 2*a + 2*d
        G0 = c*n*d / (2*np.pi*Dm)  # minimal filling factor
        G = G0 / np.cos(x*np.pi/180)  # filling factor
        B = G*(2-G)  # optical coverage
        Rgs = 4 / (sigma*c*n*(d**2)*np.pi*np.cos(x*np.pi/180))  # DC resistance
        
        k1 = (np.pi/4) / ((2/3)*G0 + np.pi/10)
        k2 = (np.pi/4) / ((2/3)*G0 + 3/8)
        te = 12*G*np.cbrt((B**2)*d/Dm)
        th = 9.6*G*np.cbrt((B**2)*d/Dm)
        dR = 0.67*d/np.sqrt(np.cos(x*np.pi/180))
        
        dataf = []
        data_wLs = []
        data_Lt = []
        data_Zr = []
        data_Zt = []
        
        freq = np.logspace(np.log10(freq_low), np.log10(freq_high), 500)
        
            
        R0 = 4 / (np.pi * d**2 * n * c * sigma * np.cos(x*np.pi/180))
        delta = np.sqrt(2/(2*np.pi*freq*u*sigma))
        gamma = (1 + 1j) / delta
        Zd = R0 * gamma * d / np.sinh(gamma * d)
        
        F = n * c * d / (2 * np.pi * Dm * np.cos(x*np.pi/180))
        tau = 9.6 * F * np.cbrt(F**2 * (2 - F)**2 * d / Dm)
        
        #w = (1 - F) * n * d / (F * np.cos(x*np.pi/180)) # check definitions of w and l--they're reversed in Vance, but gave me w > l, which should not be true for alpha < 45
        #l = (1 - F) * n * d / (F * np.sin(x*np.pi/180))
        w = 2*np.pi*Dm / c - n*d / np.cos(x*np.pi/180)
        l = w / np.tan(x*np.pi/180)
        
        if x <= 45:
            e = np.sqrt(1 - (w / l)**2)
        else:
            e = np.sqrt(1 - (l / w)**2)
            #e = np.sqrt(1 - np.arctan(x * np.pi / 180)**2)
            
        Vmy = 2 / (3 * np.sqrt(np.pi)) * (l / w)**(3/2) * ((1 - e**2) * e**2 / (ellipe(e) - (1 - e**2) * ellipk(e)))
        Sr = l * w / 2
        my = Sr**(3/2) * Vmy
        # I can't get the my values to come out correctly--Vmy seems to be throwing it off. (Possibly the dimensionless l/w coefficient.)
        K = 2*F - F**2
        Mh_Vance = np.pi*u/(6*c) * (1-K)**(3/2) * (e**2 / (ellipe(e) - (1 - e**2) * ellipk(e)))
        Mh = 0.875 * 0.5079 * np.exp(-tau) * Mh_Vance        
        
        hhat = 2*d / (1 + n*(1-F) / F)
        Mb = -u * d * hhat * (1 - np.tan(x*np.pi/180)**2) / (4 * np.pi * Dm * d)
        
        #nu = 4 * np.pi * (Dm / 2) * np.sin(x*np.pi/180) * np.cos(x*np.pi/180) * F**2 / (n**2 + d**2)  # holes per unit length
        # I think holes per unit length is already integrated into the inductance formulae
        
        Zt = Zd + 1j * 2*np.pi*freq * (Mh + Mb)
        
        rt20v = np.abs(Zd)
        mt20 = np.abs(Mh + Mb)
        

        title = '%s Transfer Impedance'% model 
        plt.figure(0)
        if NORMALIZE:
            plt.loglog(freq, np.abs(Zt)/np.abs(data_Zt[0]), label = 'Beatrics Transfer Impedance (normalized)')
        else:
            plt.loglog(freq, np.abs(Zt), label = 'Beatrics Transfer Impedance')
        #plt.loglog(dataf,zt_mag20, label = '38.6 Degrees Simulated')
        #plt.loglog(dataf,zt_sim, label = '38.6 Degrees Simulated')
        plt.title(title,fontsize=16)
        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Transfer Impedance (Ohm/m)")
        plt.xlim(1e4,2e10)  
        plt.ylim(1e-6,1e2)  
        plt.grid(True)
        
        
    print('The transfer resistance is:', rt20v[0])
    print('The transfer Inductance is:',mt20)
    
    



if __name__ == '__main__':
    model = input('Which analytical model would you like to use to compute the transfer impedance for a braided shield? (Vance/Tyni/Demoulin/Kley/Beatrics):')

    if model.lower() in ['vance', 'tyni', 'demoulin', 'kley', 'beatrics']:
        calc_braid_impedance(model, corrections=False)
    else:
        print('Enter a valid analytical model.')