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
    df = pd.read_csv(csvpath,sep='\t')
    data = []
    columns = df.columns
    print(columns)
    for i, row in df.iterrows():
        print(row)
        currdata = {'id':i}
        param = process_literal(row['parameter_change'])
        perturbation = {'pars':None,'ics':None}
        treat = {'pars':None,'ics':None}
        if type(param) == list:
            perturbation['pars'] = param[0]['parameters']
            perturbation['ics'] = param[0]['inconds']
            treat['pars'] = param[1]['parameters']
            treat['ics'] = param[1]['inconds']
        else:
            perturbation['pars'] = param['parameters']
            perturbation['ics'] = param['inconds']
            
        skeleton = {
            'simid':row['sim_id'],
            'simulate':row['simulate'],
            'description':row['description'],
            'type':row['sim_type'],
            'whichcondition':row['whichcondition'],
            'readout':row['readout'],
            'source':row['source'],
            'time':{'wt':row['tend_wt'],
                    'perturb':row['tend_mut']},
            'units':row['units'],
            'value':{'preshift':{'wt':row['wt_pre'],
                                 'perturb':row['mut_pre']},
                     'postshift':{'wt':row['wt_post'],
                                  'perturb':row['mut_post']}},
            'spec':{'preshift':{'pars':process_literal(row['pre_spec'])},
                    'postshift':{'pars':process_literal(row['post_spec'])},
                    'perturbation':{'pars':perturbation['pars'],
                                    'ics':perturbation['ics']},
                    'treat':{'pars':treat['pars'], 
                             'ics':treat['ics']},
            },
            'citation':row['citation']
        }        
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
    
