# -*- coding: utf-8 -*-
"""
Author: Jingying Hu
Created on Tue Oct 19 10:07:35 2018
CBE 5790 Midterm Project
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def elba(n0 = (31999, 0, 1, 0, 0, 0), nVadd = 0, tAdd = None, timeSpan = 120, nMax = 200000, nRun = 1):
    """
    This function is written for .
    Jingying Hu
    Created: 10/19/2018
    People who helped me: James F. Rathman, Bryan Hobocienski, Edward Hughes.
    """
    # deal with n0
    if type(n0) == str:
        raise TypeError
    # deal with nVadd
    if type(nVadd) == str:
        raise TypeError
    elif nVadd < 0:
        raise Exception('ValueError - additional doses of vaccine must be >= 0')
    # deal with timeSpan   
    if type(timeSpan) == str:
        raise TypeError
    if timeSpan < 0:
        raise Exception('Error - timeSpan must be positive')
    # deal with nMax
    if type(nMax) != int:
        raise TypeError 
    if nMax < 0:
        raise Exception('Error - maximum time steps must be postive')
    
    """================the continuous deterministic model=================="""
    
    p = [6.022e+23, 1.76e-05, 0.1, 0.01, 3.52e-06]
    def ode_model(n,t):
        H = n[0] # not immuned
        F = n[1] # free riders
        S = n[2] # sick
        I = n[3] # immuned
        D = n[4] # dead 
        V = n[5] # doses of vaccine
        # constant values of N_av, k_1, k_2, k_3, k_4                
        dHdt = -p[1]*H*S-p[4]*H*V
        dFdt = -p[1]*F*S
        dSdt = p[1]*(H+F)*S-p[2]*S-p[3]*S
        dIdt = p[2]*S+p[4]*H*V
        dDdt = p[3]*S
        dVdt = -p[4]*H*V
        return [dHdt,dFdt,dSdt,dIdt,dDdt,dVdt]    
    def Jacob(n,t): 
        H = n[0]
        F = n[1]
        S = n[2]
        I = n[3]
        D = n[4]
        V = n[5]
        return np.array([[-p[1]*S-p[4]*V, 0, -p[1]*H, 0, 0, -p[4]*H],
                         [0, -p[1]*S, -p[1]*F, 0, 0, 0],
                         [p[1]*S, p[1]*S, p[1]*H+p[1]*F-p[2]-p[3], 0, 0, 0],
                         [p[4]*V, 0, p[2], 0, 0, p[4]*H],
                         [0, 0, p[3], 0, 0, 0],
                         [-p[4]*V, 0, 0, 0, 0 -p[4]*H]])   

    # deal with tAdd
    if type(tAdd) == str:
        raise TypeError
    elif tAdd != None and tAdd < 0:
        raise Exception('ValueError - the day at which additional doses of vaccine provided can not be negative')
    elif tAdd != None and tAdd > timeSpan:
            raise Exception('Error - exceed  time span')
    elif tAdd != None and tAdd < timeSpan and tAdd > 0:
        # [continuous]two periods
        t1 = np.linspace(0,tAdd,nMax*(tAdd/timeSpan)+1) 
        t2 = np.linspace(tAdd, timeSpan, nMax*(1-tAdd/timeSpan)+2)
        tc = np.append(t1,t2)
        # [continuous]integrate first period 0 ~ t1
        nc1 = integrate.odeint(ode_model, n0, t1, Dfun = Jacob)
        nc1[np.where(t1 == tAdd),5] += nVadd #updata the number of vaccine
        # [continuous]integrate second period t1 ~ t2
        n1 = np.reshape(nc1[np.where(t1 == tAdd),:],(6,))
        nc2 = integrate.odeint(ode_model, n1, t2, Dfun = Jacob)
        nc = np.append(nc1,nc2,axis = 0)
    elif tAdd == None: 
    # [continuous]
        tc = np.linspace(0,timeSpan,nMax)
        nc = integrate.odeint(ode_model, n0, tc, Dfun = Jacob)
        
    
    """========================the stochastic model=========================""" 
    def Gillespie_model(n0, nVadd, tAdd, timeSpan, nMax):
        t = np.zeros((nMax+1,1)) # list of time steps
        n = np.zeros((nMax+1,6)) # list of H,F,S,I,D,V number variation
        n[0,:] = n0
        H = n[:,0] # each column for each variable
        F = n[:,1]
        S = n[:,2]
        I = n[:,3]
        D = n[:,4]
        V = n[:,5]
        # deal with tAdd
        if tAdd != None and tAdd < timeSpan and tAdd > 0:         
            V[tAdd] = nVadd  # [stochasitc]provide additional vaccine
        r = np.zeros((nMax+1,5)) # pre-allocate reaction rate values
        # initial rates
        r[0,:] = [p[1]*H[0]*S[0], p[1]*F[0]*S[0], p[2]*S[0], p[3]*S[0],p[4]*H[0]*V[0]]
        w = np.zeros((nMax+1,1)) # pre-allocate random number for each reactiont time
        # time length for first reaction
        w[0] = np.random.random_sample()   
        T = np.sum(r, axis = 1) # sum up of reaction rates at each time step
        prob = np.zeros((nMax+1,5)) #pre-allocate prob. for each reaction at each time step
        # initial prob. of each reaction at time zero 
        prob[0,:] = np.divide(r[0], T[0]) 
        csp = np.zeros((nMax+1,5)) # pre-allocate csp
        csp[0] = np.cumsum(prob[0]) 
        q = np.zeros((nMax+1,1)) # pre-allocate random prob. of reaction
        q[0] = np.random.random_sample() # generate prob.
        k = np.zeros((nMax+1,1)) #pre-allocate row index for stoichiometry
        k[0] = np.where(q[0] < csp[0])[0][0] # column index in a row
        sto = np.array([[-1,0,1,0,0,0],
                        [0,-1,1,0,0,0],
                        [0,0,-1,1,0,0],
                        [0,0,-1,0,1,0],
                        [-1,0,0,1,0,-1]])
        
        for i in range (1, nMax+1):
            t[i] = t[i-1] - np.log(w[i-1])/T[i-1]
            if t[i] >= timeSpan:
                break
            n[i] = n[i-1] + sto[int(k[i-1])]
            r[i,:] =  [p[1]*H[i]*S[i], p[1]*F[i]*S[i], p[2]*S[i], p[3]*S[i],p[4]*H[i]*V[i]]
            T = np.sum(r, axis = 1)
            if T[i] == 0:
                break
            prob[i] = np.divide(r[i],T[i])
            csp[i] = np.cumsum(prob[i])
            q[i] = np.random.random_sample()
            k[i] = np.where(q[i] < csp[i])[0][0]
            w[i] = np.random.random_sample()
            i = i + 1
        dt = D[i-1]
        #removed unused rows        
        Vlist = [n,t,r,prob,csp,k,w,q,dt]      
        for j in range(0,len(Vlist)-1):
            Vlist[j] = np.delete(Vlist[j], np.s_[i:nMax+2],axis = 0 )
        return Vlist

    """=======================subplot for continuous========================"""
    #sum over people who are willing to be vaccinated and free riders
    healthyc = nc[:,0] + nc[:,1]
    axc = plt.subplot(1,2,1)
    axc.plot(tc,healthyc,'b-',label = 'healthy')
    axc.plot(tc,nc[:,2],'r-', label = 'sick')
    axc.plot(tc,nc[:,3],'c-', label = 'immune')
    axc.plot(tc,nc[:,4],'m-', label = 'dead')
    axc.legend(loc='upper right')
    axc.set_xlabel('t - time (day)')
    axc.set_ylabel('n - number of individuals')
    axc.set_title('Demographic trends predicted by the continuous-variable deterministic model')
    
    #deal with nRun 
    #n0=n0,nVadd=nVadd,tAdd=tAdd,timeSpan=timeSpan,nMax=nMax
    if type(nRun) != int:
        raise TypeError
    if nRun < 1:
        raise Exception('Error - nRun must be >= 1')
    if nRun == 1:
        Vlist = Gillespie_model(n0=n0,nVadd=nVadd,tAdd=tAdd,timeSpan=timeSpan,nMax=nMax)
        healthys = Vlist[0][:,0] + Vlist[0][:,1]
        # subplot for stochastic
        axs = plt.subplot(1,2,2)
        axs.plot(Vlist[1],healthys,'b-',label = 'healthy')
        axs.plot(Vlist[1],Vlist[0][:,2],'r-', label = 'sick')
        axs.plot(Vlist[1],Vlist[0][:,3],'c-', label = 'immune')
        axs.plot(Vlist[1],Vlist[0][:,4],'m-', label = 'dead')
        axs.legend(loc='upper right')
        axs.set_xlabel('t - time (day)')
        axs.set_ylabel('n - number of individuals')
        axs.set_title('Demographic trends predicted by the stochastic model')
    if nRun > 1:
        death = np.zeros((nRun,2))
        death[:,0] = np.arange(1,nRun+1,1)        
        for run in range(0, nRun):
            death[run,1] = Gillespie_model(n0=n0,nVadd=nVadd,tAdd=tAdd,timeSpan=timeSpan,nMax=nMax)[8]
            run += 1
        axs = plt.subplot(1,2,2)
        axs.hist(death[:,1],20)
        axs.set_xlabel('final death toll')
        axs.set_ylabel('frequncy')
        axs.set_title('Histogram of the final death toll resutls obtained by stochastic model')
    return 
    
    
