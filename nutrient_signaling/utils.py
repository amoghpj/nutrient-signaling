import PyDSTool as dst
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
            home = os.path.expanduser('~')
            with open(home + 'SIMULATION_LOGS','a') as outfile:
                today = datetime.datetime.today()
                outfile.write("\nDIVIDE BY ZERO ENCOUNTERED IN MINMAXNORM on " + str(today.year) + '-' + \
                              str(today.month) + '-' + str(today.day) + \
                              '-' + str(today.hour) + '-' +str(today.minute))
            denom = 1
        else:
            denom = mx - mi
        output_list.append((ele - mi) / (denom))
        
    return(output_list)


