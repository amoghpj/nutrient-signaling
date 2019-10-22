"""
Stores the rows in a csv file as entries in a yaml file
"""
import bibtexparser as bp
import pandas as pd
import yaml
from optparse import OptionParser
import sys
import ast
import numpy as np

def csv2yaml(csvpath, yamlpath, inputtype, bibpath):
    if bibpath != '':
        with open(bibpath, 'r') as infile:
            bib = bp.load(infile)
    if inputtype == 'perturb':
        data = read_perturb(csvpath, bib)
    elif inputtype == 'time':
        data = read_timecourse(csvpath, bib)
    elif inputtype == 'experiment':
        data = read_experiment(csvpath, bib)
    else:
        print('Input type not recognized')
        sys.exit()
    with open(yamlpath, 'w') as outfile:
        outfile.write(yaml.dump(data, default_flow_style=False))

def read_time(csvpath, bib):
    data = []
    skeleton = {}
    return data

def get_doi(key,bib):
    key = key.split('][')[1].replace(']]','')
    for e in bib.entries:
        if e['ID'] == key:
            entry = e
            break
    checkfields = ['url', 'doi']
    for f in checkfields:
        if f in entry.keys():
            return entry[f]
    return(entry['ID'])
            
    
def read_experiment(csvpath, bib):
    df = pd.read_csv(csvpath)
    data = []
    for i, row in df.iterrows():
        print(row['ID'])
        preshift = {}
        if not pd.isna(row['pre_pars']):
            preshift = process_literal(row['pre_pars'])
        postshift = {}
        if not pd.isna(row['post_pars']):
            postshift = process_literal(row['post_pars'])
        mutant = {'pars':{}, 'ics':{}}
        if not pd.isna(row['mutant']):
            print(row['mutant'])
            mutspec = process_literal(row['mutant'])
            if 'pars' in mutspec.keys():
                mutant['pars'] = mutspec['pars']
            if 'ics' in mutspec.keys():                
                mutant['ics'] = mutspec['ics']

        skeleton = {'id':row['ID'],
                    'strain':row['strain'],
                    'nutrientCondition':row['Nutrient input'],
                    'phenotype':row['Phenotype'],
                    'growth':row['Growth characteristic'],
                    'expReadout':row['Experimental Readout'],
                    'background':row['Background'],
                    'comments':row['Comments'],
                    'doi':get_doi(row['citation'], bib),
                    'preshift':{k:v for k, v in preshift.items() },
                    'postshift':{k:v for k, v in postshift.items() },
                    'mutant':{'pars':{k:v for k, v in mutant['pars'].items() },
                              'ics':{k:v for k, v in mutant['ics'].items() }},
                    'simReadout':row['Simulation Readout'],
        }
        data.append(skeleton)
    return data

def read_perturb(csvpath):
    df = pd.read_csv(csvpath,sep='\t')
    data = []
    columns = df.columns
    for i, row in df.iterrows():
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
            'vmax':row['vmax'],
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
    parser.add_option('-b','--bib',dest='bib',default='',
                      help="Path to bib file")        
    (options, args) = parser.parse_args()
    csv2yaml(options.csvpath,
             options.yamlpath,
             options.inputtype,
             options.bib)

    
if __name__ == '__main__':
    main()
    
