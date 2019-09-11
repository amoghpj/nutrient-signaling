"""
Stores the rows in a csv file as entries in a yaml file
"""
import pandas as pd
import yaml
from optparse import OptionParser
import sys
import ast

def csv2yaml(csvpath, yamlpath, inputtype):
    if inputtype == 'perturb':
        data = read_perturb(csvpath)
    elif inputtype == 'time':
        data = read_timecourse(csvpath)
    elif inputtype == 'experiment':
        data = read_experiment(csvpath)
    else:
        print('Input type not recognized')
        sys.exit()
    print(data)
    with open(yamlpath, 'w') as outfile:
        outfile.write(yaml.dump(data, default_flow_style=False))

def read_time(csvpath):
    data = []
    skeleton = {}
    return data

def read_experiment(csvpath):
    data = []
    skeleton = {}
    return data

def read_perturb(csvpath):
    df = pd.read_csv(csvpath)
    data = []
    columns = df.columns
    print(columns)
    for i, row in df.iterrows():
        print(row)
        currdata = {'id':i}
        param = process_literal(row['parameter_change'])
        wt = {'pars':None,'ics':None}
        treat = {'pars':None,'ics':None}
        if type(param) == list:
            wt['pars'] = param[0]['parameters']
            wt['ics'] = param[0]['inconds']
            treat['pars'] = param[1]['parameters']
            treat['ics'] = param[1]['inconds']
        else:
            wt['pars'] = param['parameters']
            wt['ics'] = param['inconds']
            
        skeleton = {
            'description':row['description'],
            'readout':row['readout'],
            'type':row['Shift Type'],
            'source':row['source'],
            'doi':row['DOI'],
            'time':{'wt':row['WT Simulation end time'],
                    'mut':row['Mutant Simulation end time']},
            'units':row['Assay units'],
            'value':{'preshift':{'wt':row['Experimental WT pre shift level'],
                                 'mutant':row['Experimental Mutant pre shift level']},
                     'postshift':{'wt':row['Experimental WT post shift level'],
                                  'mutant':row['Experimental Mutant post shift level level']}},
            'simulations':{'preshift':{'pars':process_literal(row['Simulation pre shift specification'])},
                           'postshift':{'pars':process_literal(row['Simulation post shift specification'])},
                           'mutant':{'pars':mutant['pars'],
                                     'ics':mutant['ics']},
                           'treat':{'pars':treat['pars'], 
                                    'ics':treat['ics']},                           
            }}        
        data.append(skeleton)
    return data

def process_literal(string_specification):
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


def main():
    parser = OptionParser()
    parser.add_option('-i','--input',dest='csvpath',
                      help="Path to input csv file")
    parser.add_option('-o','--output',dest='yamlpath',
                      help="Path to output yaml file")
    parser.add_option('-t','--type',dest='inputtype',default='',
                      help="[perturb, time, experiment]")    
    (options, args) = parser.parse_args()
    csv2yaml(options.csvpath,
             options.yamlpath,
             options.inputtype)

    
if __name__ == '__main__':
    main()
    
