"""
Collection of functions to read csv tables containing experimental
data and converting them to yaml files.
Contains three "readers"
1. read_perturb()
2. read_timecourse()
3. read_experiment()

Each reader loads a csv file using pandas,
and parses out the required data by assigning values to keys.
Note that typically, the definition of a mutant or a shift
will involve specifying parameters and initial conditions. In the
csv files, they are stored a raw python dictionary strings. I use
the ast library to parse and store these in nested yaml entries.
The "skeletons" in each reader define the requisite entries for each
type of analysis.
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
        yaml.safe_dump(data, outfile,  default_flow_style=False)

def read_timecourse(csvpath, bib):
    """
    PROG: not yet complete
    """
    data = []
    df = pd.read_csv(csvpath,sep='\t')
    for i, row in df.iterrows():
        skeleton = {
            'description':row['description'],
            'title':row['title'],
            'source':row['source'],
            'citation':row['citation'],                        
            'ts_type':row['ts_type'],
            'tunits':row['tunits'],
            'readout':row['readout'],            
            'time':row['time'],
            'value':row['values'],
            'normalize':{'lower':row['lower'],
                         'upper':row['upper']},
            'spec':{'preshift':{'pars':process_literal(row['pre_shift'])['parameters'],
                                'ics':process_literal(row['pre_shift'])['inconds']},
                    'postshift':{'pars':process_literal(row['pre_shift'])['parameters'],
                                'ics':process_literal(row['pre_shift'])['inconds']}}
        }        
        data.append(skeleton)
    
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

def get_shortname(key,bib):
    key = key.split('][')[1].replace(']]','')
    for e in bib.entries:
        if e['ID'] == key:
            entry = e
            break
    first = entry['author'].split('and')[0]
    name = ''
    if ',' in first:
        name = first.split(',')[1].strip()
    elif '.' in first:
        name = ''.join([c for c in first.split(' ') if len(c) > 2])
    else:
        name = first.split(' ')[1]

    shortname = name + ' et al, '  + str(entry['year'])
    return(shortname)
            
    
def read_experiment(csvpath, bib):
    df = pd.read_csv(csvpath)
    data = []
    nancheck = lambda C : C if not pd.isna(C) else None
    nancheckliteral = lambda C : process_literal(C) if not pd.isna(C) else None    
    
    for i, row in df.iterrows():
        print(row['ID'])
        preshift = nancheckliteral(row['pre_pars'])
        postshift = nancheckliteral(row['post_pars'])
        postshiftics = nancheckliteral(row['post_ics'])
        mutant = {'pars':None, 'ics':None}
        mutspec = nancheckliteral(row['mutant'])
        if mutspec is not None:
            mutant = {p:mutspec.get(p, None) for p in ['pars','ics']}
        forGrowth = None
        if not pd.isna(row['For Growth']):
            terms = []
            forGrowth = {}
            if ',' in row['For Growth']:
                terms = row['For Growth'].split(',')
            else:
                terms = [row['For Growth']]
            forGrowth = {t.split(':')[0]:t.split(':')[1]\
                         for t in terms}
        phenoInter = None
        if not pd.isna(row['Phenotype Interpreted']):
            terms = []
            phenoInter = {}
            if ',' in row['Phenotype Interpreted']:
                terms = row['Phenotype Interpreted'].split(',')
            else:
                terms = [row['Phenotype Interpreted']]
            phenoInter = {t.split(':')[0]:t.split(':')[1]\
                         for t in terms}            
        skeleton = {'id':row['ID'],
                    'strain':row['strain'],
                    'nutrientCondition':row['Nutrient input'],
                    'phenotypeReported':row['Phenotype Reported'],
                    'phenotypeInterpreted': phenoInter,
                    'growth':row['Growth characteristic'],
                    'expReadout':row['Experimental Readout'],
                    'background':row['Background'],
                    'comments':nancheck(row['Comments']),
                    'doi':get_doi(row['citation'], bib),
                    'forGrowth':forGrowth,
                    'shortname': get_shortname(row['citation'], bib),
                    'expected':{tf.split(':')[0]:tf.split(':')[1] for tf in row['Expected State'].split(',')},
                    'preshift':preshift,
                    'postshift':postshift,
                    'postshiftics':postshiftics,
                    'mutant':{'pars':mutant['pars'],
                              'ics':mutant['ics']},
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
    
