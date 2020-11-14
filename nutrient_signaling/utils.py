import PyDSTool as dst
import pandas as pd
import nutrient_signaling. modelreader as md
import ast
import os

def processliteral(string_specification):
    """
    Wraps the `literal_eval` function from python's 
    Abstract Syntax Tree. Used in reading dictionaries 
    storing mutant specifications.

    :param string_specification: String representing a python object
    :type string_specification: str
    :return  S: Object interpreted by the `literal_eval` function.
    """
    S = ast.literal_eval(string_specification)
    return(S)

def minmaxnorm_special(input_list, mi, mx):
    """
    Like minmaxnorm, but with specified min (mi) and max (mx)

    :param input_list: list of values
    :type input_list: list
    :param mi: user defined min
    :type mi: float
    :param mx: user defined max
    :type mx: float
    :return output_list: List of scaled values
    :type output_list: list
    """
    output_list = []
    for ele in input_list:
        output_list.append((ele - mi) / (mx - mi))
    return(output_list)
    
def minmaxnorm(input_list):
    """
    Min-Max scale the values in a list. If there is no sufficient 
    difference between min and max in the list, write to log in home directory.

    :param input_list: list of values
    :type input_list: list
    :return output_list: Min-Max scaled values from input_list
    :type output_list: list
    """
    mi = min(input_list)
    mx = max(input_list)
    output_list = []
    for ele in input_list:
        if (mx-mi) < 1e-5:
            denom = 1
        else:
            denom = mx - mi
        output_list.append((ele - mi) / (denom))
        
    return(output_list)

def get_protein_abundances():
    """
    Returns dictionary of protein abundance values
    """
    scale_abundance_dict = {'Cyr1':4000,
                            'Gln1':23700,
                            'Gcn2':4500,
                            'Sak':2000,
                            'PKA':20000,
                            'Trehalase':12000,
                            'Rtg13':2300,
                            'Ras':19000,
                            'Gcn4':4600,
                            'Gln3':1600,
                            'Dot6':8000,
                            'PDE':40000,
                            'Mig1':3000,
                            'Tps1':26000,
                            'TORC1':6000,
                            'Gis1':5100,
                            'Snf1':11300,
                            'EGO':2700,
                            'Sch9':12000,
                            'EGOGAP':300}    
    return(scale_abundance_dict)
'''
Author: Amogh Jalihal
Date: 2018-10-08

Description: 
    Python implementation of the Magnus-Neudecker Eliminator-Duplicator
    matrix vectorization method for symmetric matrices

Note:
    This is a reimplementation of the R package matrixcalc which has 
    utilities for the same. Matrixcalc is authored by Frederick Novomestky,
    and the package lives at 
    https://cran.r-project.org/web/packages/matrixcalc/index.html
    
Usage:
To generate the eliminator matrix L for an nxn symmetric matrix, call L(n). Likewise for the duplicator matrix D.

'''

import numpy as np
import sys
from tqdm import tqdm
import pandas as pd
def L(n):
    '''
    Arguments:
    n: integer >=2
    Returns:
    Numpy array of dimensions n(n+1)/2 x n**2
    '''
    if type(n) != int:
        print('Incorrect type')
        sys.exit()
    if n <=2 :
        print('n is less than 2')

    k, I = u_vectors(n)

    E =  E_matrices(n)
    p = int(n*(n+1)/2)
    nsquare = n**2
    
    L = np.zeros((p,nsquare))
    for j in tqdm(range(0,n)):
        for i in range(j,n):
            L = L + np.matmul(I[np.ix_(np.arange(0,len(I)),[int(k[i][j])])],
                              E[i][j].reshape((-1,1),order='F').transpose())

    return(L)


def D(n):
    '''
    Arguments:
    n: integer >=2
    Returns:
    Numpy array of dimensions n**2 x n(n+1)/2
    '''
    if type(n) != int:
        print('Incorrect type')
        sys.exit()
    if n <=2 :
        print('n is less than 2')
    
    p = int(n*(n+1)/2)
    nsquare = n**2
    Dt = np.zeros((p,nsquare))
    k, I = u_vectors(n)
    T = T_matrices(n)
    for j in tqdm(range(0,n)):
        for i in range(j,n):
            Dt = Dt + np.matmul(I[np.ix_(np.arange(0,len(I)),[int(k[i][j])])],
                              T[i][j].reshape((-1,1),order='F').transpose())
    D = Dt.transpose()
    return D

def T_matrices(n):
    E = E_matrices(n)
    T = list()
    for i in range(0,n):
        T.append(list())
        for j in range(0,n):
            if i==j:
                T[-1].append(E[i][j])
            else:
                T[-1].append(E[i][j] + E[j][i])
    return T
                
def u_vectors(n):
    p = int(n*(n+1)/2)
    I = np.eye(p)
    k = np.zeros((n,n))
    
    for j in range(1,n+1):
        for i in range(j,n+1):
            k[i-1][j-1] = int((j-1)*n + i -0.5*(j)*(j-1)) -1
    return(k, I)


def E_matrices(n):
    I = np.eye(n)
    #print(I)
    E = list()
    for i in range(0,n):
        E.append(list())
        for j in range(0,n):
            E[-1].append(np.outer(I[i],I[j]))
    return E
