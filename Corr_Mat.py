# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 19:35:32 2021

@author: carlos

This script is responsible to correlate the grid used in the main script
"""
import numpy as np

def correlation_matrix(x_grid, y_grid, IM2):

# This variable contains the kind of the Intensity Measure and its coeficientes 
# used to correlate the grid.
# coef = [IM, alphs, beta, gamma] 

    coef  = [
        ['PGA', 0.06, 0.283, 5],
        ['SA(0.1)', 0.062, 0.276, 5],
        ['SA(0.2)', 0.073, 0.248, 5],
        ['SA(0.3)', 0.086, 0.219, 5],
        ['SA(0.5)', 0.073, 0.248, 5],
        ['SA(1.0)', 0.051, 0.329, 5],
        ['SA(2.0)', 0.061, 0.421, 3.035],
        ['SA(3.0)', 0.092, 0.671, 1.189],
        ['SA(5.0)', 0.071, 0.741, 1.201],
        ['mean', 0.054, 0.319, 5]
        ]
    
    aux = []
    for j in range(len(coef)):
        if IM2[0] == coef[j][0]:
            aux.append(j)
     
    if aux == []:
        print('ERROR: IM not implemented')
        print('Please go check IM2 or the correlation_matrix function')
        
    point = np.zeros([y_grid*x_grid, 2])
    alpha = coef[aux[0]][1]
    beta = coef[aux[0]][2]
    gamma = coef[aux[0]][3]
    
    k = 0
    for i_y in range(y_grid):
        for i_x in range(x_grid):
            
            x = i_x + 0.5
            y = i_y + 0.5
            point[k, 0] = x
            point[k, 1] = y
            
            k += 1
                
    #print(len(point))
    corr_mat = np.zeros([len(point), len(point)])
    
    for i in range(len(point)):
        for j in range(i,len(point)):
            d = np.sqrt((point[i][0]-point[j][0])**2 + (point[i][1]-point[j][1])**2)
            corr_mat[i,j] = max(gamma*np.exp(-alpha*d**beta)-gamma+1, 0)
            corr_mat[j,i] = corr_mat[i,j]
            
    #print(corr_mat)
    L = np.linalg.cholesky(corr_mat)    
    L = L.transpose()
    
    return(L)
                
#x_grid = 44
#y_grid = 43
#IM2 = ['PGA']

#teste = correlation_matrix(x_grid, y_grid, IM2)