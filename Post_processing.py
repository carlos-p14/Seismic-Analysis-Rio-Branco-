# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 17:32:39 2021

@author: carlo

This script aim to summarize the results obtained by the main code
"""

from sys import exit as DebugStop 

# first the modules
import matplotlib as mpl
import matplotlib.pyplot as plt # to plot
import matplotlib.cm as cm 
import numpy as np # to be able to operate the variables
import pickle # to load the exposure layer
import os # to be able to create directories
import time

from pylab import flipud

from coordinates2_TFC import coord_X_norm, coord_Y_norm

from Budget import Bridge_budget

#%% The main plot function
def plot(decay, hist, d_curve, three_D, budget):
    
    if decay:
        damage_decay()
        
    if hist:
        damage_hist()
        
    if d_curve:
        damage_curve()
        
    if three_D:
        plot3d()
        
    if budget:
        TVE()

#%% Plotting the damage decay
def damage_decay():
    plt.rcParams.update({'font.size':14})
    s = 7.5
    plt.figure(101)
    plt.plot(ax, light_mean, marker = 'o', ls = '--', color='k', mfc = 'b', label = 'Dano Leve', lw=2, markersize=s)
    plt.plot(ax, moderate_mean, marker = 'o', ls = '--', color='k', mfc = 'r', label = 'Dano Moderado', lw=2, markersize=s)
    plt.plot(ax, extensive_mean, marker = 'o', ls = '--', color='k', mfc = 'y', label = 'Dano Extensivo', lw=2, markersize=s)
    plt.plot(ax, complete_mean, marker = 'o', ls = '--', color='k', mfc = 'g', label = 'Dano Completo', lw=2, markersize=s)
    #plt.xticks(ax)
    plt.xlabel('Distância ao Centro da Cidade [km]')
    plt.ylabel('Nº Médio de Construções Danificadas')
    plt.title(f'Dano em Função da Distância - M{M/10: .1f}')
    plt.yticks(ay)
    plt.legend()
    props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)
    plt.annotate('$N_{máx}$ de Construções:' +f' {n_max}', xy = (0.99,0.55), xycoords = 'axes fraction',bbox = props, horizontalalignment='right', verticalalignment='top')
    plt.grid()
    plt.savefig(os.path.join(directory, f'Decaimento do dano - M{M}'))

#%% the damaged buildings histograms
def damage_hist():
    plt.rcParams.update({'font.size':12})
    plt.figure(102)
    y, x, _ = plt.hist(damage[i, :, 0], bins=15, edgecolor = 'k', alpha = .5, color = 'b', label = 'Dano Leve')
    y1 =  max(y)
    
    y, x, _ = plt.hist(damage[i, :, 1], bins=15, edgecolor = 'k', alpha = .5, color = 'r', label = 'Dano Moderado')
    y2 = max(y)
    
    y, x, _ = plt.hist(damage[i, :, 2], bins=15, edgecolor = 'k', alpha = .5, color = 'y', label = 'Dano Extensivo')
    y3 = max(y)
    
    y, x, _ = plt.hist(damage[i, :, 3], bins=15, edgecolor = 'k', alpha = .5, color = 'g', label = 'Dano Completo')
    y4 = max(y)
    
    plt.vlines(light_mean[i], 0, y1, ls='--', color = 'b', lw = 2)
    plt.vlines(light_median[i], 0, y1+1000, ls=':', color = 'b', lw = 2)
    plt.text(light_mean[i], y1, f'{light_mean[i]: .0f}', ha = 'center', va = 'bottom', color = 'b')
    plt.text(light_median[i], y1+1000, f'{light_median[i]: .0f}', ha = 'center', va = 'bottom', color = 'b')
    
    plt.vlines(moderate_mean[i], 0, y2+500, ls='--', color = 'r', lw = 2)
    plt.vlines(moderate_median[i], 0, y2+1000, ls=':', color = 'r', lw = 2)
    plt.text(moderate_mean[i], y2+500, f'{moderate_mean[i]: .0f}', ha = 'center', va = 'bottom', color = 'r')
    plt.text(moderate_median[i], y2+1000, f'{moderate_median[i]: .0f}', ha = 'center', va = 'bottom', color = 'r')
    
    plt.vlines(extensive_mean[i], 0, y3+500, ls='--', color = 'y', lw = 2)
    plt.vlines(extensive_median[i], 0, y3+1000, ls=':', color = 'y', lw = 2)
    plt.text(extensive_mean[i], y3+500, f'{extensive_mean[i]: .0f}', ha = 'center', va = 'bottom', color = 'y')
    plt.text(extensive_median[i], y3+1000, f'{extensive_median[i]: .0f}', ha = 'center', va = 'bottom', color = 'y')
    
    plt.vlines(complete_mean[i], 0, y4+500, ls='--', color = 'g', lw = 2)
    plt.vlines(complete_median[i], 0, y4+500, ls=':', color = 'g', lw = 2)
    plt.text(complete_mean[i], y4+500, f'{complete_mean[i]: .0f}', ha = 'center', va = 'bottom', color = 'g')
    plt.text(complete_median[i], y4+500, f'{complete_median[i]: .0f}', ha = 'center', va = 'bottom', color = 'g')
    
    plt.vlines(1,0,1,color='k',ls='--',label='Médias')
    plt.vlines(1,0,1,color='k',ls=':',label='Medianas')
    
    plt.grid()
    plt.legend()
    plt.xlabel('Construções Danificadas')
    plt.ylabel('Nº Aparições')
    plt.title(f'Histogramas de Estados de Dano - M{M/10: .1f}')
    plt.savefig(os.path.join(directory, f'Distribuição do Dano - M{M}'))

#%% Damage Curves
def damage_curve(): 
    plt.rcParams.update({'font.size':12})
    prob = np.arange(1, n_sim+1)/n_sim
    
    LS = ['Leve', 'Moderado', 'Extensivo', 'Completo']
    
   # text = 'OBS: Números menores do que 1 \n para construções danificadas se dão \n pelo fato de haver baixos valores de PGA \n e, portanto, um pequeno dano à ponte'
    
    if M == 57:
        end = 12
    
    else:
        end = 21
    
    for k in  range(0,end,2):
        for j in range(4):
            eLC = np.flip(np.sort(damage[k, :, j]))
            
            plt.figure(201+j)
            plt.vlines(10**0, 0, 1, ls='--', color='k', lw=2)
            plt.plot(eLC, prob, linewidth = 2, label = f'D = {ax[k]} km')
            plt.xlabel('Estruturas Danificadas')
            plt.ylabel('Probabilidade')
            plt.title(f'Curva Dano {LS[j]} - M{M/10: .1f}')
            plt.yscale('log')
            plt.xscale('log')
            plt.grid(True, which='both', ls='--', alpha = .25)
            plt.legend()
            #props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)
            #plt.annotate(text, xy = (0.175,0.3), xycoords = 'axes fraction',bbox = props, horizontalalignment='left', verticalalignment='top')
            plt.savefig(os.path.join(directory, f'Curva de Dano {LS[j]} - M{M}'))
            plt.show()       

#%% Damage Building Map
def plot3d():
    plt.rcParams.update({'font.size':14})
    Mean_Damaged_Buildings_Grid = np.zeros((num_limit_states, y_grid, x_grid))
    Median_Damaged_Buildings_Grid = np.zeros((num_limit_states, y_grid, x_grid))
    
    LS = ['Leve', 'Moderado', 'Extensivo', 'Completo']
    
    x = coord_X_norm.ravel()
    y = coord_Y_norm.ravel()
    z_pos = np.zeros_like(x)
    width = depth = 1
    
    colormap_median = 'PuBuGn' 
    colormap_mean = 'YlOrRd'
    
    for k in range(num_limit_states):
        for j in range(x_grid):
            for i in range(y_grid):
                Mean_Damaged_Buildings_Grid[k,i,j] = np.mean(Damage_map[k,:,i,j])
                Median_Damaged_Buildings_Grid[k,i,j] = np.median(Damage_map[k,:,i,j])
                
        #Mean_Damaged_Buildings_Grid[Mean_Damaged_Buildings_Grid==0] = None
        aux = Mean_Damaged_Buildings_Grid[k,:,:].ravel()
        aux2 = Median_Damaged_Buildings_Grid[k,:,:].ravel()
        
        cmap = cm.get_cmap(colormap_mean) # Get desired colormap
        max_height = np.max(aux)   # get range of colorbars
        min_height = np.min(aux)
        
        # scale each z to [0,1], and get their rgb values
        rgba = [cmap((k-min_height)/max_height) for k in aux] 
        
        norm = mpl.colors.Normalize(vmin=np.min(aux), vmax=np.max(aux))
        
        fig = plt.figure(301+k)
        ax = fig.add_subplot(111, projection = '3d')
        ax.bar3d(x, y, z_pos, width, depth, aux, shade = True, color = rgba, zsort='average', alpha = 1)
        ax.set_zlim3d(0,10)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap = colormap_mean), ax = ax)
        plt.title(f'Média de Construções Danificadas \n M{M/10: .1f} Dano {LS[k]}')
        plt.savefig(os.path.join(directory, f'Média de Construções Danificadas - M{M} Dano {LS[k]}'))
        plt.show()  
        
        cmap2 = cm.get_cmap(colormap_median)
        max_height = np.max(aux2)   # get range of colorbars
        min_height = np.min(aux2)  
        
        rgba = [cmap2((k-min_height)/max_height) for k in aux2] 
        
        norm = mpl.colors.Normalize(vmin=np.min(aux2), vmax=np.max(aux2))
        
        fig = plt.figure(401+k)
        ax = fig.add_subplot(111, projection = '3d')
        ax.bar3d(x, y, z_pos, width, depth, aux2, shade = True, color = rgba, zsort='average', alpha = 1)
        ax.set_zlim3d(0,10)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap = colormap_median), ax = ax)
        plt.title(f'Mediana de Construções Danificadas \n M{M/10: .1f} Dano {LS[k]}')
        plt.savefig(os.path.join(directory, f'Mediana de Construções Danificadas - M{M} Dano {LS[k]}'))
        plt.show() 
        
#%% Total Exposed Value
def TVE():
    plt.rcParams.update({'font.size':12.5})
    deficit = np.zeros((len(last_name), y_grid, x_grid))
    ve = np.zeros(len(last_name))
    tve = np.zeros(len(last_name))
    
    tv = np.sum(Bridge_budget)
    
    for k in range(len(last_name)):
        deficit[k, :, :] = Bridge_budget*mean_damage_map[k]
        ve[k] = np.sum(deficit[k,:,:])
        tve[k] = ve[k]/tv
    
    angle = 45
    
    plt.figure(501)
    plt.bar(ax, ve/10**6, width = 14, edgecolor = 'k')
    plt.xlabel('Distância ao Centro da Cidade [km]')
    plt.ylabel('Impacto Econômico [milhões R$]')
    plt.title('Valor Exposto')
    plt.xlim(-20, max(ax)+30)
    plt.ylim(0, max(ve/10**6)+30)
    
    for k in range(len(last_name)):
        plt.text(ax[k]+4,ve[k]/10**6 , 'R${:,.2f}'.format(ve[k]/10**6), ha = 'center', va = 'bottom', color = 'k', rotation = angle)
 
    props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)
    plt.annotate('Valor total estimado: \n R${:,.2f}'.format(tv), xy = (0.75,0.9), xycoords = 'axes fraction',bbox = props, ha = 'center', va = 'center')        
    plt.grid(True, alpha = .3)
    plt.savefig(os.path.join(directory, f'Valor Exposto - M{M}'))
    plt.show()
    
    plt.figure(502)
    plt.bar(ax, tve*100, width = 15, edgecolor = 'k', color = 'r')
    plt.xlabel('Distância ao Centro da Cidade [km]')
    plt.ylabel('Fração em Relação ao Total Exposto [%]')
    plt.title('Valor Total Exposto')
    plt.xlim(-20, max(ax)+30)
    plt.ylim(0, 100)
    
    for k in range(len(last_name)):
        plt.text(ax[k]+4,tve[k]*100 , '{:.2f}%'.format(tve[k]*100), ha = 'center', va = 'bottom', color = 'k', rotation = angle)
        
    props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)
    plt.annotate('Valor total estimado: \n R${:,.2f}'.format(tv), xy = (0.75,0.9), xycoords = 'axes fraction',bbox = props, ha = 'center', va = 'center')        
    plt.grid(True, alpha = .3)
    plt.savefig(os.path.join(directory, f'Valor Total Exposto - M{M}'))
    plt.show()   
             
#%% Creating the folder to save the figures
cd = os.getcwd() # Getting current path

directory = os.path.join(cd, 'Figuras') # Creating directory for results to be saved

if os.path.exists(directory):
    print(f"{dir} already exists") # Cheking wheter the directory already exists
    
else:
    os.mkdir(directory)
    print(f"{dir} has been created") # Creating it otherwsise

#%% Loading the results
first_name = 'DamageTotal_'
first_name2 = 'DamageMap_'

x_grid = len(coord_X_norm[0]) # number of cells on the X direction
y_grid = len(coord_X_norm)

n_max = 42

num_limit_states = 4

# Set the magnitude value: 76 or 57
M = 57
# Distance to the city's downtown 
D = 0
# Distance converted to code index
i = int(D/20)

# Number of simulations done
n_sim = 10000

if M == 57:
    last_name = [
                'D0M57',
                'D20M57',
                'D40M57',
                'D60M57',
                'D80M57',
                'D100M57',
                'D120M57',
                'D140M57',
                'D160M57',
                'D180M57',
                'D200M57',
                'D220M57'
                # 'D240M57',
                # 'D260M57',
                # 'D280M57',
                # 'D300M57',
                # 'D320M57',
                # 'D340M57',
                # 'D360M57',
                # 'D380M57',
                # 'D400M57'
                ]

elif M == 76:
    last_name = [
                'D0M76',
                'D20M76',
                'D40M76',
                'D60M76',
                'D80M76',
                'D100M76',
                'D120M76',
                'D140M76',
                'D160M76',
                'D180M76',
                'D200M76',
                'D220M76',
                'D240M76',
                'D260M76',
                'D280M76',
                'D300M76',
                'D320M76',
                'D340M76',
                'D360M76',
                'D380M76',
                'D400M76'
                ] 
else:
    print('Magnitude not defined')
    DebugStop()

damage = np.zeros((len(last_name), 10000, 4))

dm = np.zeros((len(last_name), 10000, y_grid, x_grid))

for k in range(len(last_name)):

    
    file_name = first_name + last_name[k] + '.pkl'
    another_file = first_name2 + last_name[k] + '.pkl'
    
    print(f'Reading: \n {file_name} \n {another_file} \n \n')
    
    with open(file_name,'rb') as file:
     Damage_total = pickle.load(file)

    with open(another_file, 'rb') as second_file:
        Damage_map = pickle.load(second_file)
 
    damage[k, :, :] = Damage_total[:,:]
    dm[k,:,:,:] = Damage_map[3,:,:,:]
    
another_file = first_name2 + last_name[0] + '.pkl'

with open(another_file, 'rb') as second_file:
    Damage_map = pickle.load(second_file)

    
#%% Doing some previous set ups

# getting the path where the figures will be saved
cd = os.getcwd()

mean_damage_map = np.zeros((len(last_name), y_grid, x_grid))

# setting the variables up
light_mean = damage.mean(axis = 1)[:,0]

moderate_mean = damage.mean(axis = 1)[:,1]

extensive_mean = damage.mean(axis = 1)[:,2]

complete_mean = damage.mean(axis = 1)[:,3]

light_median = np.median(damage, axis = 1)[:,0]

moderate_median = np.median(damage, axis = 1)[:,1]

extensive_median = np.median(damage, axis = 1)[:,2]

complete_median = np.median(damage, axis = 1)[:,3]

for k in range(len(last_name)):
    mean_damage_map[k, :, :] = dm.mean(axis=1)[k]

if M == 57:
    ax = np.arange(0,221,20)
else:    
    ax = np.arange(0,401,20)

ay = np.arange(0,46,5)

#%% Calling the plot function
decay = 0
hist = 0
d_curve = 0
three_D = 0
budget = 1

plot(decay, hist, d_curve, three_D, budget) 