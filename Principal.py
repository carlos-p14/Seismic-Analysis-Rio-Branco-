# -*- coding: utf-8 -*-
"""
Created on Sat Jun 26 20:11:04 2021

@author: carlos

main code
"""
#%% Importing modules and variables
# Uncomment the next line to be able to stop the code wherever you want
from sys import exit as DebugStop 

# first the modules
import matplotlib.pyplot as plt # to plot
import numpy as np # to be able to operate the variables
import pickle # to load the exposure layer
import math # some useful math tools
import os # to be able to create directories
import time

from pylab import flipud
# importing the variables we're gonna need
# Exposure Model
# loading the coordinates
from coordinates2_TFC import  coord_X_norm, coord_Y_norm

# and the exposure_layer
with open('Exposure_layer_TFC.pkl','rb') as file:
    exposure_layer = pickle.load(file)

# Fragility Model
from FragilityFunctions_TFC import Fragility_Functions, IM2

# The function used to correlate the grid
from Corr_Mat import correlation_matrix

# The function used to interpolate the Fragility Functions values
from Find_interpolate import interpolate

# The function that will evaluate the GMPEs
from GMPEs import GMPE_calc 

#%% Creating the directory where the results will be stored
cd = os.getcwd() # Getting current path

directory = os.path.join(cd, 'Results_TFC') # Creating directory for results to be saved

if os.path.exists(directory):
    print(f"{dir} already exists") # Cheking wheter the directory already exists
    
else:
    os.mkdir(directory)
    print(f"{dir} has been created") # Creating it otherwsise
#%% Scenario-based earthquake assessment input
start = time.time()

# Earthquake input
M2 = 5.7 # Moment magnitude
x2 = -370# x coordinate
y2 = 10 # y coordinate

# Number of simulations
n_sim = 10000

# Ground Motion Prediction Equation (GMPE)
# 1st column: use the GMPE; 2nd col: the associated name
GMPE = [
        [0, 0], # Campbell and Bozorgnia (2008); NGA-West
        [0, 0], # Abrahamson and Silva (1997)
        [0, 0], # Akkar e Bommer (2010); for Europe (recommend by GEM)
        [0, 'ASK13'], # Abrahamson et al. (2014); NGA-West2
        [0, 'CB13'], # Campbell and Bozorgnia (2014); NGA-West2
        [0, 'BSSA13'], # Boore et al. (2014); NGA-West2
        [0, 'CY13'],  # Chiou and Youngs (2014); NGA-West2 (recommend by GEM)
        [1, 'TFC']
        ]

i = 0

while GMPE[i][0] == 0:
    i += 1
    
GMPE_name = GMPE[i][1]
#%% To plot Fragility Functions

# typology = 0 #  it's the typology of the fragility to be plotted (it goes from 0 to 4)

# plt.figure(1)

# label = ['Leve','Moderado','Extensivo','Completo']

# for j in range(1,5):
#     plt.plot(Fragility_Functions[typology][0], Fragility_Functions[typology][j], label = f'{label[j-1]}', linewidth=2)

# plt.grid()
# plt.xlabel('PGA')
# plt.ylabel('P[D$\geq$C|IM]')
# plt.legend()
# plt.xlim(0,max(Fragility_Functions[typology][0]))
# plt.ylim(0,1.01)

# plt.show()
# plt.savefig(os.path.join(directory, 'Curvas de Fragilidade'))

# DebugStop()
#%% Previous set ups

x_grid = len(exposure_layer[0][1][0]) # number of cells on the X direction
y_grid = len(exposure_layer[0][1]) # number of cells on the Y direction
num_typologies = len(exposure_layer)-2 # number of typologies

# # Generation of the spatial correlation matrix between the cells of the city
corr_L = correlation_matrix(x_grid,y_grid,IM2)

# # Number of limit states (equal 4: slight, moderate, extensive, complete)
num_limit_states = len(Fragility_Functions[0])-1

# Generation of the inter-event (or between-event) variability;
# This is an important variable!! It is part of the development of the
# GMPEs, and it is directly related to the uncertainties between
# earthquakes
between_ev_var = np.random.normal(0,1,size=(n_sim,1)) # Random normally distributed (mean = 0 and standard deviation = 1)

# # Pre-alocate space before simulations loop
intra_map = np.zeros((y_grid,x_grid))

Damaged_buildings = np.zeros((4,y_grid,x_grid)) # total damaged buildings in each cell

Damaged_buildings2 = np.zeros((4,num_typologies, y_grid, x_grid)) # total damaged buildings in each cell per typology

Damage_map = np.zeros((num_limit_states, n_sim,y_grid,x_grid))

Damage_map2 = np.zeros((num_limit_states, num_typologies, n_sim, y_grid,x_grid))

Damage_total = np.zeros((n_sim,num_limit_states))

Damage_total2 = np.zeros((num_limit_states, num_typologies, 1, n_sim))

PGA_map = np.zeros((n_sim,y_grid,x_grid))

#%% Perfirm the scenario-based damage assessment

for i in range(n_sim):
    print(f'Simulation {i}')
    
    # Intra-event Variability
    intra_var = np.random.normal(0, 1, size = (x_grid*y_grid,1))
    intra_var2 = np.matmul(np.transpose(intra_var),corr_L)
    
    k = 0   
    
    for i_y in range(y_grid):
        for i_x in range(x_grid):
            intra_map[i_y][i_x] = intra_var2[0][k]
            k = k+1
            
    Damaged_buildings = np.zeros((4,y_grid,x_grid))
    Damaged_buildings2 = np.zeros((4,num_typologies, y_grid, x_grid))
    
    for i_x in range(x_grid):
        for i_y in range(y_grid):
            if exposure_layer[len(exposure_layer)-2][1][i_y][i_x] != 0:
                # Previous calculation
                x = i_x + 0.5 # X coordinate of the center of the cell
                y = i_y + 0.5 # Y coordinate of the center of the cell
                
                dist = np.sqrt((x2 - x)**2 + (y2 - y)**2) # Distance between earthquake epicenter and the cell
                R_JB = dist # Joyne-Boore distance, assumed = dist

                within_ev_var = intra_map[i_y][i_x]
                
                mean_ln_Y, phi_ss, tau_total = GMPE_calc(GMPE, M2, R_JB, IM2, x, y, x2, y2)
                ln_Y = mean_ln_Y + within_ev_var*phi_ss + between_ev_var[i,0]*tau_total
                Y = math.exp(ln_Y)
                #DebugStop()
                
                # Calculating damage
                for k in range(num_typologies):
                    if exposure_layer[k][1][i_y][i_x] !=0:
                        for j in range(num_limit_states):
                            
                            Damage_vector = interpolate(Fragility_Functions[k,0,:],Fragility_Functions[k,j+1,:], Y)
                            Damaged_buildings[j, i_y, i_x] = exposure_layer[k][1][i_y][i_x]*Damage_vector + Damaged_buildings[j, i_y, i_x]
                            Damaged_buildings2[j, k, i_y, i_x] = exposure_layer[k][1][i_y][i_x]*Damage_vector
                            
                PGA_map[i,i_y,i_x] = Y
                
    for j in range(num_limit_states):
        Damage_total[i,j] = np.sum(Damaged_buildings[j, :, :])
        Damage_map[j,i,:,:] = Damaged_buildings[j,:,:]
        
        Damage_total2[j,:,0,i] = np.sum(Damaged_buildings2[j,:,:,:])
        
        for k in range(num_typologies):
            Damage_map2[j,k,i,:,:] = Damaged_buildings2[j,k,:,:]

#%% Post Processing 

#%% Saving the final variables
with open('PGA_D380M57.pkl','wb') as file:
    pickle.dump(PGA_map,file)

with open('DamageMap_D380M57.pkl','wb') as file:
    pickle.dump(Damage_map,file)
    
with open('DamageTotal_D380M57.pkl','wb') as file:
    pickle.dump(Damage_total,file)
    
#%% Loading the results when the simulations were done
# with open('PGA.pkl','rb') as file:
#     PGA_map = pickle.load(file)
    
# with open('DamageMap.pkl','rb') as file:
#     Damage_map = pickle.load(file)
    
# with open('DamageTotal.pkl','rb') as file:
#     Damage_total = pickle.load(file)
#%% Exposure Map
# exposure_map = np.zeros((y_grid, x_grid))

# for j in range(x_grid):
#     for i in range(y_grid):
#         exposure_map[i,j] = exposure_layer[1][1][i][j]
        
# exposure_map[exposure_map==0] = np.nan

# plt.figure(101)
# plt.pcolor(flipud(exposure_map), cmap = 'turbo', edgecolors = 'k', linewidth = 1)
# plt.tight_layout()
# plt.gca().set_aspect('equal', adjustable='box')
# plt.colorbar(label = 'Número de Construções')

# plt.scatter(x2, np.max(coord_Y_norm + (-y2)), s = 100, marker = '*', color = 'k')
# plt.savefig(os.path.join(directory, 'Mapa de Exposição D0M76'))
# plt.show()

# #%% PGA map        
# Mean_PGA_map = np.zeros((y_grid, x_grid)) 
# std_PGA_map = np.zeros((y_grid, x_grid)) 
# cv_PGA_map = np.zeros((y_grid, x_grid)) 

# for j in range(x_grid):
#     for i in range(y_grid):
#         Mean_PGA_map[i,j] = np.mean(PGA_map[:,i,j])
#         std_PGA_map[i,j] = np.std(PGA_map[:,i,j])
#         cv_PGA_map[i,j] = std_PGA_map[i,j]/Mean_PGA_map[i,j]
        
# Mean_PGA_map[Mean_PGA_map==0] = np.nan

# plt.figure(102)
# plt.pcolor(flipud(Mean_PGA_map), cmap = 'turbo', edgecolors = 'k', linewidth = 1)
# plt.tight_layout()
# plt.gca().set_aspect('equal', adjustable='box')
# plt.colorbar(label = 'Mean PGA')

# plt.scatter(x2, np.max(coord_Y_norm + (-y2)), s = 100, marker = '*', color = 'k')
# plt.savefig(os.path.join(directory, 'Mapa de PGA - D0M76'))
# plt.show()

# #%% Damaged Buildings Map
# Mean_Damaged_Buildings_Grid = np.zeros((num_limit_states, y_grid, x_grid))
# Std_Damaged_Buildings_Grid = np.zeros((num_limit_states, y_grid, x_grid))
# CV_Damaged_Buildings_Grid = np.zeros((num_limit_states, y_grid, x_grid))

LS = ['Leve', 'Moderado', 'Extensivo', 'Completo']

# for k in range(num_limit_states):
#     for j in range(x_grid):
#         for i in range(y_grid):
#             Mean_Damaged_Buildings_Grid[k,i,j] = np.mean(Damage_map[k,:,i,j])
#             Std_Damaged_Buildings_Grid[k,i,j] = np.std(Damage_map[k,:,i,j]) 
#             CV_Damaged_Buildings_Grid[k,i,j] = Std_Damaged_Buildings_Grid[k,i,j]/Mean_Damaged_Buildings_Grid[k,i,j]
            
#     Mean_Damaged_Buildings_Grid[Mean_Damaged_Buildings_Grid==0] = np.nan
#     aux = Mean_Damaged_Buildings_Grid[k,:,:]
    
#     plt.figure(103+k)
#     plt.pcolor(flipud(aux), cmap = 'turbo', edgecolors = 'k', linewidth = 1)
#     plt.tight_layout()
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.colorbar()
#     plt.title(f'Mean Damage Buildings - {LS[k]}')
#     plt.scatter(x2, np.max(coord_Y_norm + (-y2)), s = 100, marker = '*', color = 'k')
#     plt.savefig(os.path.join(directory, f'Mapa de Dano {LS[k]} (D0M76)'))
#     plt.show()      
   
# #%% Histogram and damage curves of damaged buildings

Mean_Damaged_Buildings = Damage_total.mean(axis=0)
Std_Damaged_Buildings = Damage_total.std(axis=0)
CV_Damaged_Buildings =  Std_Damaged_Buildings/Mean_Damaged_Buildings
Relative_Damaged_Buildings = Mean_Damaged_Buildings/exposure_layer[-1][1]

# for k in range(num_limit_states):
    
#     mu_l = Damage_total.mean(axis=0)[k]
#     med_l = np.median(Damage_total, axis=0)[k]
#     prob = np.arange(1, n_sim+1)/n_sim
#     eLC = np.flip(np.sort(Damage_total[:,k]))#/exposure_layer[2][1]
    
#     str1 = f'Média = {mu_l:.2f}'
#     str2 = f'Mediana = {med_l:.2f}'
        
#     plt.figure(107+k)
#     y, x, _ = plt.hist(Damage_total[:,k], bins=15, edgecolor = 'k', label = 'Histogram')
#     y1 = [0, max(y)]
#     #x1 = [0, max(Fragility_Functions[k][0])]
#     plt.vlines(mu_l, y1[0], y1[1], ls='--', color = 'red', linewidth = 2, label = 'Média')
#     plt.vlines(med_l, y1[0], y1[1], ls = '--', color = 'yellow', linewidth = 2, label = 'Mediana')
#     plt.legend()
#     plt.xlabel('Nº de Construções Danificadas')
#     plt.ylabel('Contagem')
#     plt.title(f'Histograma - {LS[k]}')
#     props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)
#     plt.annotate(f'{str1} \n{str2}', xy = (0.975,0.75), xycoords = 'axes fraction',bbox = props, horizontalalignment='right', verticalalignment='top')
#     plt.grid()
#     plt.savefig(os.path.join(directory, f'Histograma - {LS[k]} (D0M76)'))
#     plt.show()
    
#     plt.figure(201+k)
#     plt.plot(eLC, prob, color = 'k', linewidth = 2, label = 'Curva de Dano')
#     plt.xlabel('Estruturas Danificadas')
#     plt.ylabel('Probabilidade')
#     plt.title(f'Curva Teste de Dano {LS[k]}')
#     plt.yscale('log')
#     plt.xscale('log')
#     plt.xlim(max(10**0,min(eLC)), exposure_layer[2][1]+10)
#     plt.vlines(mu_l, y1[0], 1, ls='--', color = 'red', linewidth = 2, label = 'Média')
#     plt.vlines(med_l, y1[0], 1, ls = '--', color = 'yellow', linewidth = 2, label = 'Mediana')
#     plt.grid(True, which='both', ls='--', alpha = .25)
#     plt.legend()
#     plt.savefig(os.path.join(directory, f'Curva de Dano {LS[k]} (D0M76)'))
#     plt.show()
         
#%% Saving output
file_name = 'Resultados_D380M57.txt'

file_path = os.path.join(directory, file_name)

file_text = open(file_path, 'w')

file_text.write("--- RESULTADOS ---\n\n"
                f'GMPE: {GMPE_name}\n')

for i in range(num_limit_states):
    file_text.write(
                    f"## Estado Limite {LS[i]}\n\n"
                    f"Total Médio de Construções Danificadas: {Mean_Damaged_Buildings[i]: .4f}\n"
                    f"Desvio Padrão: {Std_Damaged_Buildings[i]: .4f}\n"
                    f"Coeficiente de Variação: {CV_Damaged_Buildings[i]: .4f}\n"
                    f"Construções Danificadas / Construções Existentes: {Relative_Damaged_Buildings[i]: .4f}\n\n")

file_text.close()

end = time.time()

print(f'\nTime elapsed: {end - start}\n')


            