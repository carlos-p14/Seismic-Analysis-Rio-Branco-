# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 20:36:16 2021

@author: Carlos

This is the main side function! It is responsible for GMPEs evaluation
"""
import numpy as np
import math

def GMPE_calc(GMPE, M, R_JB, IM2, x, y, x2, y2):
    #%% Fault input
    delta = 40 # Default value as per the USGS and Kaklamanos et al. (2011) for reverse faults
    Vs30 = 180 # Shear wave velocity in the site (760 = B/C soil boundary)
    
    #%% Determination of some preliminary parameters
    W = 10**(-1.61+0.41*M) # Reverse events
    
    Z_HYP = 11.24-0.2*M
    
    # Depth to top distance (Kaklamanos et al. (2011))
    Z_TOR = max(Z_HYP-0.6*W*np.sin(math.radians(delta)), 0) # (Kaklamanos et al. (2011))
    Z_TOR2 = (max(2.704-1.226*max(M-5.849,0),0))**2 # reverse and reverse-oblique faulting (Chiou and Youngs (2014))
    
    # Depth to Vs = 1.0 km/s (Kaklamanos et al. (2011))
    if Vs30 < 180:
        Z_1 = math.exp(6.745)
    
    elif Vs30 <= 500:
        Z_1 = math.exp(6.745-1.35*math.log(Vs30/180))
        
    else:
        Z_1 = math.exp(5.394-4.48*math.log(Vs30/500))
        
    Z_1_2 = math.exp(-7.15/4*math.log((Vs30**4 + 571**4)/(1360**4+571**4))) # for California and non-japan regions (Chiou and Youngs (2014))
    
    # Depth to Vs = 2.5 km/s
    Z_25 = (519+3.595*Z_1)/1000
    
    #%% Determination of source geometric parameters
    
    Rx = x - x2 # cell coordinate - source coordinate
    Ry = y - y2
    
    if Rx < 0:
        if y2 > y:
            alfa = -(np.arctan(np.absolute((x-x2)/(y-y2))))
        else:
            alfa = -(math.pi - np.arctan(np.absolute((x-x2)/(y-y2))))
    else:
        if y2 > y:
            alfa = np.arctan(np.absolute((x-x2)/(y-y2)))
        else:
            alfa = math.pi - np.arctan(np.absolute((x-x2)/(y-y2)))
            
    # Determine R_RUP
    if Rx < Z_TOR*np.tan(math.radians(delta)):
        R_RUP2 = np.sqrt(Rx**2 + Z_TOR**2)
        
    elif Rx <= Z_TOR*np.tan(math.radians(delta)) + W*np.sin(math.radians(delta)):
        R_RUP2 = Rx*np.sin(math.radians(delta)) + Z_TOR*np.cos(math.radians(delta))
        
    else:
        R_RUP2 = np.sqrt((Rx-W*np.cos(math.radians(delta)))**2 + (Z_TOR+W*np.sin(math.radians(delta)))**2)
        
    R_RUP = np.sqrt(R_RUP2**2 + Ry**2)
    
    #%% Ground Motion predictions equations     
    
    #%% Abrahamson and Silva 2014 GMPE
    def ASK13():  
        FRV = 1
        FN = 0
        
        coef = [
            ['PGA',100,6.75,660,-1.47,2.4,4.5,0.587,-0.79,0.275,-0.1,-0.41,2.154,0,-0.015,1.735,0,-0.1,0.6,-0.3,1.1,-0.0072,0.1,0.05,0,-0.05,-0.0015,0.0025,-0.0034,-0.1503,0.265,0.337,0.188,0,0.088,-0.196,0.044,0.754,0.52,0.47,0.36,0.741,0.501,0.54,0.63],
            ['PGV',1000,6.75,330,-2.02,2400,4.5,5.975,-0.919,0.275,-0.1,-0.41,2.366,0,-0.094,2.36,0,-0.1,0.25,0.22,0.3,-0.0005,0.28,0.15,0.09,0.07,-0.0001,0.0005,-0.0037,-0.1462,0.377,0.212,0.157,0,0.095,-0.038,0.065,0.662,0.51,0.38,0.38,0.66,0.51,0.58,0.53],
            ['SA(0.01)',0.01,6.75,660,-1.47,2.4,4.5,0.587,-0.79,0.275,-0.1,-0.41,2.154,0,-0.015,1.735,0,-0.1,0.6,-0.3,1.1,-0.0072,0.1,0.05,0,-0.05,-0.0015,0.0025,-0.0034,-0.1503,0.265,0.337,0.188,0,0.088,-0.196,0.044,0.754,0.52,0.47,0.36,0.741,0.501,0.54,0.63],
            ['SA(0.02)',0.02,6.75,680,-1.46,2.4,4.5,0.598,-0.79,0.275,-0.1,-0.41,2.146,0,-0.015,1.718,0,-0.1,0.6,-0.3,1.1,-0.0073,0.1,0.05,0,-0.05,-0.0015,0.0024,-0.0033,-0.1479,0.255,0.328,0.184,0,0.088,-0.194,0.061,0.76,0.52,0.47,0.36,0.747,0.501,0.54,0.63],
            ['SA(0.03)',0.03,6.75,770,-1.39,2.4,4.5,0.602,-0.79,0.275,-0.1,-0.41,2.157,0,-0.015,1.615,0,-0.1,0.6,-0.3,1.1,-0.0075,0.1,0.05,0,-0.05,-0.0016,0.0023,-0.0034,-0.1447,0.249,0.32,0.18,0,0.093,-0.175,0.162,0.781,0.52,0.47,0.36,0.769,0.501,0.55,0.63],
            ['SA(0.05)',0.05,6.75,915,-1.22,2.4,4.5,0.707,-0.79,0.275,-0.1,-0.41,2.085,0,-0.015,1.358,0,-0.1,0.6,-0.3,1.1,-0.008,0.1,0.05,0,-0.05,-0.002,0.0027,-0.0033,-0.1326,0.202,0.289,0.167,0,0.133,-0.09,0.451,0.81,0.53,0.47,0.36,0.798,0.512,0.56,0.65],
            ['SA(0.075)',0.075,6.75,960,-1.15,2.4,4.5,0.973,-0.79,0.275,-0.1,-0.41,2.029,0,-0.015,1.258,0,-0.1,0.6,-0.3,1.1,-0.0089,0.1,0.05,0,-0.05,-0.0027,0.0032,-0.0029,-0.1353,0.126,0.275,0.173,0,0.186,0.09,0.506,0.81,0.54,0.47,0.36,0.798,0.522,0.57,0.69],
            ['SA(0.1)',0.1,6.75,910,-1.23,2.4,4.5,1.169,-0.79,0.275,-0.1,-0.41,2.041,0,-0.015,1.31,0,-0.1,0.6,-0.3,1.1,-0.0095,0.1,0.05,0,-0.05,-0.0033,0.0036,-0.0025,-0.1128,0.022,0.256,0.189,0,0.16,0.006,0.335,0.81,0.55,0.47,0.36,0.795,0.527,0.57,0.7],
            ['SA(0.15)',0.15,6.75,740,-1.59,2.4,4.5,1.442,-0.79,0.275,-0.1,-0.41,2.121,0,-0.022,1.66,0,-0.1,0.6,-0.3,1.1,-0.0095,0.1,0.05,0,-0.05,-0.0035,0.0033,-0.0025,0.0383,-0.136,0.162,0.108,0,0.068,-0.156,-0.084,0.801,0.56,0.47,0.36,0.773,0.519,0.58,0.7],
            ['SA(0.2)',0.2,6.75,590,-2.01,2.4,4.5,1.637,-0.79,0.275,-0.1,-0.41,2.224,0,-0.03,2.22,0,-0.1,0.6,-0.3,1.1,-0.0086,0.1,0.05,0,-0.03,-0.0033,0.0027,-0.0031,0.0775,-0.078,0.224,0.115,0,0.048,-0.274,-0.178,0.789,0.565,0.47,0.36,0.753,0.514,0.59,0.7],
            ['SA(0.25)',0.25,6.75,495,-2.41,2.4,4.5,1.701,-0.79,0.275,-0.1,-0.41,2.312,0,-0.038,2.77,0,-0.1,0.6,-0.24,1.1,-0.0074,0.1,0.05,0,0,-0.0029,0.0024,-0.0036,0.0741,0.037,0.248,0.122,0,0.055,-0.248,-0.187,0.77,0.57,0.47,0.36,0.729,0.513,0.61,0.7],
            ['SA(0.3)',0.3,6.75,430,-2.76,2.4,4.5,1.712,-0.79,0.275,-0.1,-0.41,2.338,0,-0.045,3.25,0,-0.1,0.6,-0.19,1.03,-0.0064,0.1,0.05,0.03,0.03,-0.0027,0.002,-0.0039,0.2548,-0.091,0.203,0.096,0,0.073,-0.203,-0.159,0.74,0.58,0.47,0.36,0.693,0.519,0.63,0.7],
            ['SA(0.4)',0.4,6.75,360,-3.28,2.4,4.5,1.662,-0.79,0.275,-0.1,-0.41,2.469,0,-0.055,3.99,0,-0.1,0.58,-0.11,0.92,-0.0043,0.1,0.07,0.06,0.06,-0.0023,0.001,-0.0048,0.2136,0.129,0.232,0.123,0,0.143,-0.154,-0.023,0.699,0.59,0.47,0.36,0.644,0.524,0.66,0.7],
            ['SA(0.5)',0.5,6.75,340,-3.6,2.4,4.5,1.571,-0.79,0.275,-0.1,-0.41,2.559,0,-0.065,4.45,0,-0.1,0.56,-0.04,0.84,-0.0032,0.1,0.1,0.1,0.09,-0.002,0.0008,-0.005,0.1542,0.31,0.252,0.134,0,0.16,-0.159,-0.029,0.676,0.6,0.47,0.36,0.616,0.532,0.69,0.7],
            ['SA(0.75)',0.75,6.75,330,-3.8,2.4,4.5,1.299,-0.79,0.275,-0.1,-0.41,2.682,0,-0.095,4.75,0,-0.1,0.53,0.07,0.68,-0.0025,0.14,0.14,0.14,0.13,-0.001,0.0007,-0.0041,0.0787,0.505,0.208,0.129,0,0.158,-0.141,0.061,0.631,0.615,0.47,0.36,0.566,0.548,0.73,0.69],
            ['SA(1.0)',1,6.75,330,-3.5,2.4,4.5,1.043,-0.79,0.275,-0.1,-0.41,2.763,0,-0.11,4.3,0,-0.1,0.5,0.15,0.57,-0.0025,0.17,0.17,0.17,0.14,-0.0005,0.0007,-0.0032,0.0476,0.358,0.208,0.152,0,0.145,-0.144,0.062,0.609,0.63,0.47,0.36,0.541,0.565,0.77,0.68],
            ['SA(1.5)',1.5,6.75,330,-2.4,2.4,4.5,0.665,-0.79,0.275,-0.1,-0.41,2.836,0,-0.124,2.6,0,-0.1,0.42,0.27,0.42,-0.0022,0.22,0.21,0.2,0.16,-0.0004,0.0006,-0.002,-0.0163,0.131,0.108,0.118,0,0.131,-0.126,0.037,0.578,0.64,0.47,0.36,0.506,0.576,0.8,0.66],
            ['SA(2.0)',2,6.75,330,-1,2.4,4.5,0.329,-0.79,0.275,-0.1,-0.41,2.897,0,-0.138,0.55,0,-0.1,0.35,0.35,0.31,-0.0019,0.26,0.25,0.22,0.16,-0.0002,0.0003,-0.0017,-0.1203,0.123,0.068,0.119,0,0.083,-0.075,-0.143,0.555,0.65,0.47,0.36,0.48,0.587,0.8,0.62],
            ['SA(3.0)',3,6.82,330,0,2.4,4.5,-0.06,-0.79,0.275,-0.1,-0.41,2.906,0,-0.172,-0.95,0,-0.1,0.2,0.46,0.16,-0.0015,0.34,0.3,0.23,0.16,0,0,-0.002,-0.2719,0.109,-0.023,0.093,0,0.07,-0.021,-0.028,0.548,0.64,0.47,0.36,0.472,0.576,0.8,0.55],
            ['SA(4.0)',4,6.92,330,0,2.4,4.5,-0.299,-0.79,0.275,-0.1,-0.41,2.889,0,-0.197,-0.95,0,-0.1,0,0.54,0.05,-0.001,0.41,0.32,0.23,0.14,0,0,-0.002,-0.2958,0.135,0.028,0.084,0,0.101,0.072,-0.097,0.527,0.63,0.47,0.36,0.447,0.565,0.76,0.52],
            ['SA(5.0)',5,7,330,0,2.4,4.5,-0.562,-0.765,0.275,-0.1,-0.41,2.898,0,-0.218,-0.93,0,-0.1,0,0.61,-0.04,-0.001,0.51,0.32,0.22,0.13,0,0,-0.002,-0.2718,0.189,0.031,0.058,0,0.095,0.205,0.015,0.505,0.63,0.47,0.36,0.425,0.568,0.72,0.5],
            ['SA(6.0)',6,7.06,330,0,2.4,4.5,-0.875,-0.711,0.275,-0.1,-0.41,2.896,0,-0.235,-0.91,0,-0.2,0,0.65,-0.11,-0.001,0.55,0.32,0.2,0.1,0,0,-0.002,-0.2517,0.215,0.024,0.065,0,0.133,0.285,0.104,0.477,0.63,0.47,0.36,0.395,0.571,0.7,0.5],
            ['SA(7.5)',7.5,7.15,330,0,2.4,4.5,-1.303,-0.634,0.275,-0.1,-0.41,2.87,0,-0.255,-0.87,0,-0.2,0,0.72,-0.19,-0.001,0.49,0.28,0.17,0.09,0,0,-0.002,-0.14,0.15,-0.07,0,0,0.151,0.329,0.299,0.457,0.63,0.47,0.36,0.378,0.575,0.67,0.5],
            ['SA(10.0)',10,7.25,330,0,2.4,4.5,-1.928,-0.529,0.275,-0.1,-0.41,2.843,0,-0.285,-0.8,0,-0.2,0,0.8,-0.3,-0.001,0.42,0.22,0.14,0.08,0,0,-0.002,-0.0216,0.092,-0.159,-0.05,0,0.124,0.301,0.243,0.429,0.63,0.47,0.36,0.359,0.585,0.64,0.5]
            ]
   
        if Rx < 0: # Zero hanging wall
            FHW = 0 # 1 = sites over hanging walls 0 otherwise
        
        else:
            FHW = 1 
            
        for i in range(len(IM2)):
            
            k = 0
            
            while k <= len(coef):
                if IM2[i] == coef[k][0]:
                    Period = k
                    k = len(coef)
                
                k = k + 1
                
            M2 = 5 # See the paper
            n = 1.5 # See the paper
            a2HW = 0.2 # See the paper
            M1 = coef[Period][2]
            VLin = coef[Period][3]
            b = coef[Period][4]
            c = coef[Period][5]
            c4 = coef[Period][6]
            a1 = coef[Period][7]
            a2 = coef[Period][8]
            a3 = coef[Period][9]
            a4 = coef[Period][10]
            a5 = coef[Period][11]
            a6 = coef[Period][12]
            a7 = coef[Period][13]
            a8 = coef[Period][14]
            a10 = coef[Period][15]
            a11 = coef[Period][16]
            a12 = coef[Period][17]
            a13 = coef[Period][18]
            a14 = coef[Period][19]
            a15 = coef[Period][20]
            a17 = coef[Period][21]
            a43 = coef[Period][22]
            a44 = coef[Period][23]
            a45 = coef[Period][24]
            a46 = coef[Period][25]
            a25 = coef[Period][26]
            a28 = coef[Period][27]
            a29 = coef[Period][28]
            a31 = coef[Period][29]
            a36 = coef[Period][30]
            a37 = coef[Period][31]
            a38 = coef[Period][32]
            a39 = coef[Period][33]
            a40 = coef[Period][34]
            a41 = coef[Period][35]
            a42 = coef[Period][36]
            s1e = coef[Period][37]
            s2e = coef[Period][38]
            s3 = coef[Period][39]
            s4 = coef[Period][40]
            s1m = coef[Period][41]
            s2m = coef[Period][42]
            s5 = coef[Period][43]
            s6 = coef[Period][44]
            
            s1 = s1e
            s2 = s2e
            
            for j in range(1):
                
                # Magnitude term
                if M > 5:
                    c4m = c4
                    
                elif M > 4:
                    c4m = c4 - (c4-1)*(5-M)
                    
                else:
                    c4m = 1
                    
                R = np.sqrt(R_RUP**2 + c4m**2)
                
                if M > M1: # Magnitude term
                    f1 = a1 + a5*(M-M1) + a8*(8.5-M)**2 + (a2 + a3*(M - M1))*math.log(R) + a17*R_RUP
                
                elif M >= M2 and M < M1: 
                    f1 = a1 + a4*(M-M1) + a8*(8.5-M)**2 + (a2 + a3*(M - M1))*math.log(R) + a17*R_RUP
                    
                else:
                    f1 = a1 + a4*(M2 - M1) + a8*(8.5 - M2)**2 + a6*(M - M2) + a7*(M - M2)**2 + (a2 + a3*(M2 - M1))*math.log(R) + a17*R_RUP2
                    
                # Style of Faulting Term
                if M > 5:
                    f7 = a11
                    f8 = a12
                    
                elif M >= 4:
                    f7 = a11*(M - 4)
                    f8 = a12*(M - 4)
                    
                else:
                    f7 = 0
                    f8 = 0
                    
                # Hanging Wall Term (Model)
                R1 = W*np.cos(delta)
                R2 = 3*R1
                Ry1 = Rx*np.tan(math.radians(20))
                h1 = 0.25
                h2 = 1.5
                h3 = -0.75
                
                if delta > 30:
                    T1 = (90-delta)/45
                    
                else:
                    T1 = 60/45
                    
                if M > 6.5:
                    T2 = 1 + a2HW*(M-6.5)
                    
                elif M > 5.5:
                    T2 = 1 + a2HW*(M - 6.5) - (1 - a2HW)*(M - 6.5)**2
                    
                else:
                    T2 = 0
                    
                if Rx < R1:
                    T3 = h1 + h2*(Rx/R1) + h3*(Rx/R1)**2
                    
                elif Rx <= R2:
                    T3 = 1 - (Rx - R1)/(R2 - R1)
                    
                else:
                    T3 = 0
                    
                if Z_TOR <= 10:
                    T4 = 1 - Z_TOR**2/100
                    
                else:
                    T4 = 0
                
                if R_JB == 0:
                    T5 = 1
                    
                elif R_JB < 30:
                    T5 = 1 - R_JB/30
                    
                else:
                    T5 = 0
                    
                f4 = a13*T1*T2*T3*T4*T5
                
                # Depth to Rupture
                if Z_TOR < 20:
                    f6 = a15*Z_TOR/20
                
                else:
                    f6 = a15
                
                # Soil Depth Model
                Z_1Ref = 1/1000*math.exp(-7.67/4*math.log((Vs30**4+610**4)/(1360**4+610**4))) # For California
                                
                if Vs30 < 200:
                    f10 = a43*math.log((Z_1 + 0.01)/(Z_1Ref + 0.01))
                    
                elif Vs30 < 300:
                    f10 = a44*math.log((Z_1 + 0.01)/(Z_1Ref + 0.01))
                    
                elif Vs30 < 500:
                    f10 = a45*math.log((Z_1 + 0.01)/(Z_1Ref + 0.01))
                    
                else:
                    f10 = a46*math.log((Z_1 + 0.01)/(Z_1Ref + 0.01))
                    
                ln_Y2 = f1 + FRV*f7 + FN*f8 + FHW*f4 + f6 + f10 # Aftershock is not considered and regional attenuation teerm no included; see the paper
                Sa1180 = math.exp(ln_Y2)
                
                # Site response term (model)
                if M1 <= 0.5:
                    V1 = 1500
                
                elif M1 < 3:
                    V1 = math.exp(-0.35*math.log(M1/0.5)+math.log(1500))
                    
                else:
                    V1 = 800
                    
                if Vs30 < V1:
                    Vs30_2 = Vs30
                    
                else:
                    Vs30_2 = V1
                                
                if Vs30 >= VLin:
                    f5 = (a10 + b*n)*math.log(Vs30_2/VLin)
                
                else:
                    f5 = a10*math.log(Vs30_2/VLin) - b*math.log(Sa1180 + c) + b*math.log(Sa1180 + c*(Vs30_2/VLin)**n)
                    
                # Variability
                if M < 4: # Linear intra-event standard deviation for soil
                    PhiAL = s1
                    
                elif M <= 6:
                    PhiAL = s1 + (s2-s1)/2*(M-4)
                    
                else:
                    PhiAL = s2
                                   
                if M < 5: # Inter-event variability for the linear range
                    TauAL = s3
                
                elif M <= 7:
                    TauAL  = s3 + (s4 - s3)/2*(M-5)
                    
                else:
                    TauAL = s4
                    
                PhiAmp = 0.4 # Standard Deviation of the site amplification
                Phi_B = math.sqrt((PhiAL**2 - PhiAmp**2))
                
                if Vs30 >= VLin:
                    FACTOR = 0
                    
                else: 
                    FACTOR = (-b*Sa1180 + c) + (b*Sa1180)/(Sa1180 + c*(Vs30/VLin)**2)
                                   
                intraevent_var = math.sqrt((Phi_B**2*(1 + FACTOR)**2+PhiAmp**2))
                interevent_var = TauAL*(1+FACTOR)
                ln_Y = f1 + FRV*f7 + FN*f8 + f5 + FHW*f4 + f6 + f10 # Aftershock is not considered and regional attenuation teerm no included; see the paper
                                
                return(ln_Y, intraevent_var, interevent_var)
            
#%% Campbell and Bozorgnia 2014
    def CB13():
       
        FRV = 1
        FNM = 0
        SJ = 0
        
        coef = ['PGA', 100, -4.416,   0.984,    0.537,   -1.499,   -0.496,   -2.773,   0.248,   6.768,   -0.212,   0.720,    1.090,    2.186,   1.420,   -0.0064,   -0.202,   0.393,    0.0977,   0.0333,   0.00757,   -0.0055,   0.0000,   0.167,   0.241,   1.474,   -0.715,   -0.337,   -0.270,    865,   -1.186,   1.839,   0.734,   0.492,   0.409,   0.322,   0.166,    1.000]
        
        h4 = 1
        n = 1.18
        c = 1.88
        c8 = 0
        c0 = coef[2]
        c1 = coef[3]
        c2 = coef[4]
        c3 = coef[5]
        c4 = coef[6]
        c5 = coef[7]
        c6 = coef[8]
        c7 = coef[9]
        c9 = coef[10]
        c10 = coef[11]
        c11 = coef[12]
        c12 = coef[13]
        c13 = coef[14]
        c14 = coef[15]
        c15 = coef[16]
        c16 = coef[17]
        c17 = coef[18]
        c18 = coef[19]
        c19 = coef[20]
        c20 = coef[21]
        Dc20 = coef[22]
        a2 = coef[23]
        h1 = coef[24]
        h2 = coef[25]
        h3 = coef[26]
        h5 = coef[27]
        h6 = coef[28]
        k1 = coef[29]
        k2 = coef[30]
        k3 = coef[31]
        phi1 = coef[32]
        phi2 = coef[33]
        tau1 = coef[34]
        tau2 = coef[35]
        phiC = coef[36]
        rholny = coef[37]
        
        for i in range(1):
            
            # Magnitude Term
            if M <= 4.5:
                fmag = c0 + c1*M
            
            elif M <= 5.5:
                fmag = c0 + c1*M + c2*(M - 4.5)
                
            elif M <= 6.5:
                fmag = c0 + c1*M + c2*(M - 4.5) + c3*(M - 5.5)
            
            else:
                fmag = c0 + c1*M + c2*(M - 4.5) + c3*(M - 5.5) + c4*(M - 6.5)
                
            # Geometric attenuation term
            fdis = (c5 + c6*M)*math.log(math.sqrt(R_RUP**2 + c7**2))
            
            # Style-of-faulting term
            ff1tF = c8*FRV + c9*FNM
            
            if M <= 4.5:
                ff1tM = 0
                
            elif M <= 5.5:
                ff1tM = M - 4.5
                
            else:
                ff1tM = 1
                
            ff1t = ff1tF*ff1tM
            
            # Hanging-Wall term
            R1 = W*np.cos(delta/180*math.pi)
            R2 = 62*M - 350
            
            if Rx < 0:
                f_hng_Rx = 0
                
            elif Rx < R1:
                f1 = h1 + h2*(Rx/R1) + h3*(Rx/R1)**2
                f_hng_Rx = f1
                
            else:
                f2 = h4 + h5*(Rx-R1)/(R2-R1) + h6*((Rx-R1)/(R2-R1))**2
                f_hng_Rx = max(f2, 0)
                
            if R_RUP == 0:
                f_hng_RRup = 1
                
            elif R_RUP > 0:
                f_hng_RRup = (R_RUP - R_JB)/R_RUP
                
            if M <= 5.5:
                f_hng_M = 0
                
            elif M < 6.5:
                f_hng_M = (M - 5.5)*(1+a2*(M-6.5))
                
            else:
                f_hng_M = 1 + a2*(M-6.5)
                
            if Z_TOR > 16.66:
                f_hng_Z = 0
                
            else:
                f_hng_Z = 1 -0.06*Z_TOR
                
            f_hng_delta = (90-delta)/45
            fhng = c10*f_hng_Rx*f_hng_RRup*f_hng_M*f_hng_Z*f_hng_delta
            
            # Basin response term
            if Z_25 <= 1:
                fsed = (c14+c15*SJ)*(Z_25-1)
                
            elif Z_25 <= 3:
                fsed = 0
                
            else:
                fsed = c16*k3*math.exp(-0.75)*(1 - math.exp(-0.25*(Z_25-3)))

            # Hypocentral depth term
            if Z_HYP <=7:
                fhypH = 0
                
            elif Z_HYP <= 20:
                fhypH = Z_HYP - 7
                
            else:
                fhypH = 13
                
            if M <= 5.5:
                fhypM = c17
                
            elif M <= 6.5:
                fhypM = c17 + (c18-c17)*(M-5.5)
                
            else:
                fhypM = c18
                
            fhyp = fhypH*fhypM
            
            # Fault dip term
            if M <= 4.5:
                fdip = c19*delta
                
            elif M <= 5.5:
                fdip = c19*(5.5-M)*delta
                
            else:
                fdip = 0
                
            # Anelastic attenuation term
            if R_RUP > 80:
                fatn = (c20 + Dc20)*(R_RUP-80)
                
            else:
                fatn = 0
                
            ln_Y2 = fmag + fdis + ff1t + fhng + fsed + fhyp + fdip + fatn
            
            A1100 = math.exp(ln_Y2)
            
            # Shallow site response term
            if Vs30 <= k1:
                fsiteG = c11*math.log(Vs30/k1) + k2*(math.log(A1100 + c*(Vs30/k1)**n)-math.log(A1100+c))
              
            else: 
                fsiteG = (c11 + k2*n)*math.log(Vs30/k1)

            if Vs30 <= 200:
                fsiteJ = (c12 + k2*n)*(math.log(Vs30/k1)- math.log(200/k1))
                
            else:
                fsiteJ = (c13 + k2*n)*math.log(Vs30/k1)
                
            fsite = fsiteG + fsiteJ*SJ   
            
            # Variability
            if M <= 4.5:
                tau_lnY = tau1
                phi_lnY = phi1
                tau_lnPGAB = coef[34]
                phi_lnPGAB = coef[32]
                
            elif M < 5.5:
                tau_lnY = tau2 + (tau1 - tau2)*(5.5 - M)
                phi_lnY = phi2 + (phi1 - phi2)*(5.5 - M)
                tau_lnPGAB = coef[35] + (coef[34] - coef[35])*(5.5 - M)
                phi_lnPGAB = coef[33] + (coef[32] - coef[33])*(5.5-M)
                
            else:
                tau_lnY = tau2
                phi_lnY = phi2
                tau_lnPGAB = coef[35]
                phi_lnPGAB = coef[33]
                
            if Vs30 < k1:
                alfa = k2*A1100*((A1100 + c*(Vs30/k1)**n)**(-1) - (A1100+c)**(-1))
                
            else:
                alfa = 0
                
            phi_lnAF = 0.3
            tau_lnyb = tau_lnY
            phi_lnyb = math.sqrt(phi_lnY**2 - phi_lnAF**2)
                
            interevent_var = (tau_lnyb**2 + alfa**2*tau_lnPGAB**2 + 2*alfa*rholny*tau_lnyb*tau_lnPGAB)**0.5
            intraevent_var = math.sqrt(phi_lnyb**2 + phi_lnAF**2 + alfa**2*phi_lnPGAB**2 + 2*alfa*rholny*phi_lnyb*phi_lnPGAB)
                
            ln_Y = fmag + fdis + ff1t + fhng + fsite + fsed + fhyp + fdip + fatn
                
            return(ln_Y, intraevent_var, interevent_var)
            
#%% Boore et al. (2014)
    def BSSA13():
       
        U = 0
        SS = 0
        NS = 0
        RS = 1
        
        coef = ['PGA',       0,      0.447300,    0.485600,    0.245900,    0.453900,   1.431000,    0.050530,   -0.166200,   5.50,   -1.134000,   0.191700,   -0.008088,   4.50,   -0.0025500,   -0.600000,   1500.00,   -0.150000,   -0.007010,   -9.900,   -9.900,   110.000,   270.000,   0.100,   0.070,   0.695,   0.495,   0.398,   0.348]
        
        Mref = 4.5
        f_1 = 0
        V1 = 225
        V2 = 300
        Rref = 1
        Vref = 760
        f3 = 0.1
        e0 = coef[2]
        e1 = coef[3]
        e2 = coef[4] 
        e3 = coef[5]
        e4 = coef[6]
        e5 = coef[7]
        e6 = coef[8]
        Mh = coef[9]
        c1 = coef[10] 
        c2 = coef[11]
        c3 = coef[12]
        h = coef[13]
        Dc3 = coef[14]
        c = coef[15]
        Vc = coef[16]
        f4 = coef[17] 
        f5 = coef[18]
        f6 = coef[19]
        f7 = coef[20]
        R1 = coef[21]
        R2 = coef[22]
        DfR = coef[23]
        DfV = coef[24]
        f1 = coef[25]
        f2 = coef[26]
        t1 = coef[27]
        t2 = coef[28]
        R = math.sqrt(R_JB**2 + h**2)
        
        # PGAR
        if M <= Mh:
            FE = coef[2]*U + coef[3]*SS + coef[4]*NS + coef[5]*RS + coef[6]*(M-Mh) + coef[7]*(M-Mh)**2
            
        else:
            FE = coef[2]*U + coef[3]*SS + coef[4]*NS + coef[5]*RS + coef[8]*(M-Mh)
            
        FP = (coef[10] + coef[11]*(M-Mref))*math.log(R/Rref) + (coef[12] + Dc3)*(R-Rref)
        PGAr = math.exp(FE+FP)
        
        # Source Function
        if M <= Mh:
            FE = e0*U + e1*SS + e2*NS + e3*RS + e4*(M-Mh) + e5*(M-Mh)**2
            
        else:
            FE = e0*U + e1*SS + e2*NS + e3*RS + e6*(M-Mh) 
            
        # Path function
        FP = (c1 + c2*(M-Mref))*math.log(R/Rref) + (c3 + Dc3)*(R - Rref)
        
        # Site function
        if Vs30 <= Vc:
            ln_F1in = c*math.log(Vs30/Vref)
            
        else:
            ln_F1in = c*math.log(Vc/Vref)
            
        f_2 = f4*(math.exp(f5*(min(Vs30,760)-360)) - math.exp(f5*(760-360)))
        ln_Fn1 = f_1 + f_2*math.log((PGAr + f3)/f3)
        
        muz1 = math.exp(-7.15/4*math.log((Vs30**4+570.94**4)/(1360**4 + 570.94**4)) - math.log(1000))
        Dz1 = Z_1 - muz1
        
        if coef[1] < 0.65:
            FDz1 = 0
        
        elif coef[1] >= 0.65 and Dz1 <= f7/f6:
            FDz1 = f6*Dz1
            
        else:
            FDz1 = f7
            
        FS = ln_F1in + ln_Fn1 + FDz1
        
        # inter-event var
        if M <= 4.5:
            interevent_var = t1
            
        elif M < 5.5:
            interevent_var = t1 + (t2 - t1)*(M - 4.5)
        
        else: 
            interevent_var = t2
            
        #intra-event var
        if M <= 4.5:
            phi_M = f1
            
        elif M < 5.5:
            phi_M = f1 + (f2 - f1)*(M - 4.5)
            
        else:
            phi_M = f2
            
        if R_JB <= R1:
            phi_M_R = phi_M
            
        elif R_JB >= R1 and R_JB <= R2:
            phi_M_R = phi_M + DfR*math.log((R_JB/R1)/(R2/R1))
            
        elif R_JB > R2:
            phi_M_R = phi_M + DfR
            
        if Vs30 >= V2:
            phi_M_R_V = phi_M_R
            
        elif Vs30 >= V1:
            phi_M_R_V = phi_M_R - DfV*math.log((V2/Vs30)/(V2/V1))
            
        else:
            phi_M_R_V = phi_M_R - DfV
            
        intraevent_var = phi_M_R_V
            
        ln_Y = FE + FP + FS
    
        return(ln_Y, intraevent_var, interevent_var)    

#%% Chiou and Youngs (2014)
    def CY13():
        FRV = 1
        FNM = 0
        
        coef = ['PGA',       100,   -1.5065,  0.165,  -0.255,  -0.165,  0.255,  16.0875,  4.9993,  1.06,  1.9636,  -2.1,  -0.5,  50,  6.4551,  3.0956,  0.4908,  0.0352,   0.0462,  0.,     0.2695,  0.4833,  0.9228,  0.1202,  6.8607,  0.,      -0.4536,    -0.007146,  -0.006758,  4.2542,  -0.521,   -0.1417,   -0.00701,   0.102151,  0.,     300,  1.5817,  0.7594,  -0.6846,  0.459,    800.,        0.4,     0.26,    0.4912,  0.3762,  0.8,     0.4528]
        
        c2 = 1.06
        c4 = -2.1
        c4a = -0.5
        crb = 50
        c8a = 0.2695
        c11 = 0
        phi6 = 300
        phi6jp = 800
        DDPP = 0
        c1 = coef[2] 
        c1a = coef[3]
        c1b = coef[4]
        c1c = coef[5]
        c1d = coef[6]
        cn = coef[7]
        cm = coef[8]
        c2 = coef[9]
        c3 = coef[10] 
        c4 = coef[11]
        c4a = coef[12]
        crb = coef[13]
        c5 = coef[14]
        chm = coef[15]
        c6 = coef[16]
        c7 = coef[17]
        c7b = coef[18]
        c8 = coef[19]
        c8a = coef[20]
        c8b = coef[21]
        c9 = coef[22]
        c9a = coef[23]
        c9b = coef[24]
        c11 = coef[25]
        c11b = coef[26]
        cg1 = coef[27]
        cg2 = coef[28] 
        cg3 = coef[29]
        phi1 = coef[30] 
        phi2 = coef[31]
        phi3 = coef[32]
        phi4 = coef[33]
        phi5 = coef[34]
        phi6 = coef[35]
        gjpit = coef[36]
        gwn = coef[37]
        phi1jp = coef[38]
        phi5jp = coef[39]
        phi6jp = coef[40]
        tau1 = coef[41]
        tau2 = coef[42]
        sig1 = coef[43]
        sig2 = coef[44]
        sig3 = coef[45]
        sig2jp = coef[46]
        
        if Rx < 0:
            FHW = 0
        
        else:
            FHW = 1

        DZ_TOR = Z_TOR - Z_TOR2
        DZ_1 = Z_1 - Z_1_2
               
        ln_Y_ref = c1 + (c1a + c1c/np.cosh(2*max(M-4.5,0)))*FRV \
            + (c1b + c1d/np.cosh(2*max(M-4.5,0)))*FNM \
            + (c7 + c7b/np.cosh(2*max(M-4.5,0)))*DZ_TOR \
            + (c11 + c11b/np.cosh(2*max(M-4.5,0)))*np.cos(delta/180*math.pi)**2 \
            + c2*(M-6) + (c2-c3)/cn*math.log(1 + math.exp(cn*(cm-M)))   \
            + c4*math.log(R_RUP + c5*np.cosh(c6*max(M-chm,0)))   \
            + (c4a - c4)*math.log(math.sqrt(R_RUP**2 + crb**2)) \
            + (cg1 + cg2/np.cosh(max(M-cg3,0)))*R_RUP   \
            + c8*max(1-max(R_RUP-40,0)/30,0)*min(max(M-5.5,0)/0.8,1)*math.exp(-c8a*(M-c8b)**2)*DDPP \
            + c9*FHW*np.cos(delta/180*math.pi)*(c9a + (1-c9a)*np.tanh(Rx/c9b))*(1 - math.sqrt(R_JB**2 + Z_TOR**2)/(R_RUP + 1))

        Y_ref = math.exp(ln_Y_ref)
        
        ln_Y = ln_Y_ref \
            + phi1*min(math.log(Vs30/1130),0) \
            + phi2*(math.exp(phi3*(min(Vs30,1130)-360)) - math.exp(phi3*(1130-360)))*math.log((Y_ref + phi4)/phi4) \
            + phi5*(1-math.exp(-DZ_1/phi6))
            
        # Variability
        tau = tau1 + (tau2-tau1)/1.5*(min(max(M,5),6.5)-5)
        NL0 = phi2*(math.exp(phi3*(min(Vs30,1130)-360))-math.exp(phi3*(1130-360)))*(Y_ref/(Y_ref + phi4))
        sigmaNL0 = sig1 + (sig2 - sig1)/1.5*(min(max(M,5),6.5)-5)*math.sqrt(sig3+0.7+(1+NL0)**2)
        
        interevent_var = (1 + NL0)*tau
        intraevent_var = sigmaNL0   
        
        return(ln_Y, intraevent_var, interevent_var)
        
#%% Defining evaluation engine 
          
    if GMPE[3][0] == 1: # Abrahamson and Silva 2014 GMPE
        mean_ln_Y_ASK13, phi_ss_ASK13, tau_total_ASK13 = ASK13()
        
        return(mean_ln_Y_ASK13, phi_ss_ASK13, tau_total_ASK13)
    
    elif GMPE[4][0] == 1: # Campbell and Bozorgnia 2014
        mean_ln_Y_CB13, phi_ss_CB13, tau_total_CB13 = CB13() 
       
        return(mean_ln_Y_CB13, phi_ss_CB13, tau_total_CB13)
    
    elif GMPE[5][0] == 1:
        mean_ln_Y_BSSA13, phi_ss_BSSA13, tau_total_BSSA13 = BSSA13() 
       
        return(mean_ln_Y_BSSA13, phi_ss_BSSA13, tau_total_BSSA13)        
      
    elif GMPE[6][0] == 1:
        mean_ln_Y_CY13, phi_ss_CY13, tau_total_CY13 = CY13() 
       
        return(mean_ln_Y_CY13, phi_ss_CY13, tau_total_CY13)  

    elif GMPE[7][0] == 1:
        mean_ln_Y_ASK13, phi_ss_ASK13, tau_total_ASK13 = ASK13()
        
        mean_ln_Y_CB13, phi_ss_CB13, tau_total_CB13 = CB13()
        
        mean_ln_Y_BSSA13, phi_ss_BSSA13, tau_total_BSSA13 = BSSA13()
        
        mean_ln_Y_CY13, phi_ss_CY13, tau_total_CY13 = CY13()
        
        mean_ln_Y = 0.25*(mean_ln_Y_ASK13 + mean_ln_Y_CB13 + mean_ln_Y_BSSA13 + mean_ln_Y_CY13) 
        
        phi_ss = 0.25*(phi_ss_ASK13 + phi_ss_CB13 + phi_ss_BSSA13 + phi_ss_CY13)
        
        tau_total = 0.25*(tau_total_ASK13 + tau_total_CB13 + tau_total_BSSA13 + tau_total_CY13)
        
        return(mean_ln_Y, phi_ss, tau_total)
        
    else: 
        print('Im sorry. This GMPE has not been implemented yet') 
        
        
        
        
        
        
        
        
        
        
        
        