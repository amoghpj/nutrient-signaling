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
