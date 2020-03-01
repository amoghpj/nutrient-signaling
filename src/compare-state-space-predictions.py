import os
import sys
import yaml
import time
import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path
import multiprocessing as mp
from ctypes import c_wchar_p, c_char_p, c_char, c_wchar
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
# Local Imports
from nutrient_signaling.simulators import get_simulator

home = os.path.expanduser('~')
READOUTS = ['Gis1','Mig1','Dot6','Gcn4','Rtg13','Gln3']
SIMULATOR = 'cpp'
SIMARGS = {'executable':'robust.o'}
DESTINATION = 'tmp/'

Path(DESTINATION).mkdir(parents=True, exist_ok=True)

#####################################
##
def createNutrientInputSpace():
    """
    For each combination of high and low carbon and nitrogen
    input, create a dictionary whose keys are the nutrient input id (nid)
    and strings corresponding to the two letter abbreviation of a state
    	C	Gln	NH4	Pro	name
    0	0	0	0	0	LCLN
    1	0	1	0	0	LCHG
    2	0	0	1	0	LCHN
    3	0	0	0	1	LCHP
    4	1	0	0	0	HCLN
    5	1	1	0	0	HCHG
    6	1	0	1	0	HCHN
    7	1	0	0	1	HCHP
    """
    carbon_states = [{'C':0},{'C':1}]
    nitrogen_states = [
        {'G':0,'N':0,'P':0},
        {'G':1,'N':0,'P':0},
        {'G':0,'N':1,'P':0},
        {'G':0,'N':0,'P':1}]
    
    all_unique_states = []
    
    ## Create all combinations of low and high nutrient input states
    for c in carbon_states:
        for n in nitrogen_states:
            s = dict()
            s.update(c)
            s.update(n)
            all_unique_states.append(s)
    
    NutrientStates = []
    for nid, q in enumerate(all_unique_states):
        nstate = {}
        nstate = {
            'C':q['C'],
            'Gln':q['G'],
            'NH4':q['N'],
            'Pro':q['P']
        }
    
        C = 'LC'
        N = 'LN'
    
        if q['C'] == 1:
            C = 'HC'
        if q['G'] == 1:
            N = 'HG'
        if q['N'] == 1:
            N = 'HN'
        if q['P'] == 1:
            N = 'HP'
        nstate['name'] = C + N
        NutrientStates.append(nstate)
    return(NutrientStates)


def getMutantSpecs(modelpath):
    simargs = dict(SIMARGS)
    model = get_simulator(modelpath=modelpath,
                          simulator=SIMULATOR,**simargs)
    
    nodel = ['cAMP', 'Protein', 'Rib', \
             'Glutamine','Gis1','Rtg13',
             'Gcn4','Gln3','Mig1',
             'Gln1','Dot6','Trehalase','eIF',
             'Tps1']
    remove_interaction = [['w_gcn_torc','Sch9-inhibit-gcn2'],
                          ['w_gln_snf','Snf1-activates-nGln3'],
                          ['w_gln_sit','TORC1-inhibits-Gln3'],
                          ['w_pka_sch9','Sch9-inhibits-PKA'],
                          ['w_ras_pka','PKA-inhibits-Ras'],
    ]
    mutationSpecs = [{'name':'wt','pars':{},'ics':{}}]
    
    for v in model.get_variables():
        if v not in nodel:
            mutationSpecs.append({'pars':{v+'_T':0.0},'ics':{v:0.0},'name':v})
            
    for p in remove_interaction:
        if p[0] in model.get_pars():
            mutationSpecs.append({'pars':{p[0]:0.0},'ics':{},'name':p[1]})
    return(mutationSpecs)


def simulateGlobalStateSpace(args): 
    """
    PROG  Update shared dict of global state
    """
    NutrientStates = args['NutrientStates']
    strain = args['strain']
    print(strain['name'])
    modelpath = args['modelpath']
    fname = args['fname']
    cutoffs = args['cutoffs']
    paramdict = args['paramdict']
    
    simargs = dict(SIMARGS)

    if SIMULATOR == 'cpp':
        simargs['simfilename'] = cleanfname(strain['name']) + '.dat'
        
    model = get_simulator(modelpath=modelpath,
                          simulator=SIMULATOR,**simargs)
        
    collect = {}
    model.set_attr(pars=paramdict)
    
    for ninput in NutrientStates:
        shiftspec = {'pars':{'Carbon':ninput['C'],
                             'ATP':ninput['C'],
                             'Glutamine_ext':ninput['Gln'],
                             'Proline':ninput['Pro'],
                             'NH4':ninput['NH4']}}
        deletion_par = strain['pars']
        deletion_par.update(shiftspec['pars'])

        deletion_ics = strain['ics']
        model.set_attr(pars=deletion_par,
                       ics=deletion_ics)
        ss = model.get_ss()
        ss.update(deletion_ics)

        ss_b = tfStates(ss, READOUTS, cutoffs)
        outstr = ''.join([str(ss_b[r]['dec']) for r in READOUTS])
        collect[ninput['name']] = {strain['name']:outstr}
        
    collectdf = pd.DataFrame(collect)
    collectdf.to_csv(DESTINATION + fname)


def define_state_space(modelpath):
    statedict = {r:[] for r in READOUTS}
    
    nut = [{'Carbon':0.0,'ATP':0.0,'Glutamine_ext':0.0},
           {'Carbon':0.0,'ATP':0.0,'Glutamine_ext':1.0},
           {'Carbon':1.0,'ATP':1.0,'Glutamine_ext':0.0},
           {'Carbon':1.0,'ATP':1.0,'Glutamine_ext':1.0},
           {'Carbon':0.0,'ATP':0.0,'Glutamine_ext':0.0, 'Proline':0.0,'NH4':1.0},
           {'Carbon':1.0,'ATP':1.0,'Glutamine_ext':0.0, 'Proline':0.0,'NH4':1.0},
           {'Carbon':0.0,'ATP':0.0,'Glutamine_ext':0.0, 'Proline':1.0,'NH4':0.0},
           {'Carbon':1.0,'ATP':1.0,'Glutamine_ext':0.0, 'Proline':1.0,'NH4':0.0}
    ]
    model = get_simulator(modelpath=modelpath,
                          simulator=SIMULATOR,**SIMARGS)
    for n in nut:
        model.set_attr(pars=n, tdata=[0,200])
        P = model.get_ss()
        for k in statedict.keys():
            statedict[k].append(P[k])
            
    cutoffs = {}
    for k, v in statedict.items():
        if k in ['Mig1','Dot6']:
            lo = np.log10((min(v))/(1-min(v)))
            hi = np.log10((max(v))/(1-max(v)))
        else:
            lo = min(v)
            hi = max(v)
        mean = (lo+hi)/2.
        print(k, round(lo,3), round(hi,3), round(mean,3))
        cutoffs[k] = mean
    return statedict, cutoffs

def to_file(path, aggregate, colorder):
    """
    Write predictions to file
    - TODO Refactor to using pandas
    - NOTE This function is wholly unnecessary, a better
      approach would be to first create the data object, 
      then write it out to file using pandas.
    """
    if Path(path).exists():
        outstr = ''
    else:
        ## Write header
        outstr =  '\t'.join([col for col in colorder]) + '\n'
    for col in colorder:
        print(col)
        if col != 'cost':
            strain = col.split('_')[0]
            ninput = col.split('_')[1]
            outstr += aggregate[strain][ninput] + '\t'
        else:
            outstr += aggregate[col] + '\t'
    outstr += '\n'
    
    with open(path, 'a+') as outfile:
        outfile.write(outstr)
        
def simulateStrains(modelpath,
                    parameterSetPath='',
                    debug=False):
    
    NutrientStates = createNutrientInputSpace()
    mutationSpecs = getMutantSpecs(modelpath)
    aggregator = {'cost':''}
    colorder = ['cost']
    for strain in mutationSpecs:
        aggregator[strain['name']] = {}
        for ninput in NutrientStates:
            aggregator[strain['name']][ninput['name']] = ''
            colorder.append(strain['name'] + '_' + ninput['name'])
            
    if not debug:
        sd, cutoffs = define_state_space(modelpath)
    else:
        cutoffs = {'Gis1':0.498,
                   'Mig1':1.344,
                   'Dot6':1.058,
                   'Gcn4':0.489,
                   'Rtg13':0.46,
                   'Gln3':0.564}        
    parind = 0
    if parameterSetPath is not None:
        if Path(parameterSetPath).exists():
            psetdf = pd.read_csv(parameterSetPath)
        else:
            print('Invalid path')
            sys.exit()
        psets = psetdf.T.to_dict()
        psetlength = len(psets)
    else:
        psets = [None]
        psetlength = 1
    clock = time.time()
    outpath = 'output/global_space_{}.csv'.format(psetlength)
    while parind < psetlength:
        print(parind)

        if debug:
            for strain in mutationSpecs:
                expargs = {'strain':strain,
                           'NutrientStates':NutrientStates,
                           'fname':strain['name'],                           
                           'cutoffs': cutoffs,
                           'paramdict':psets[parind],
                           'modelpath':modelpath}                
                simulateGlobalStateSpace(expargs)
        else:
            workers = []
            counter = 0
            fnames = []
            jobs = []
            with mp.Pool() as pool:
                for strain in mutationSpecs:
                    fnames.append(strain['name'] + '.txt')
                    expargs = {'strain':strain,
                               'NutrientStates':NutrientStates,
                               'cutoffs': cutoffs,
                               'fname':fnames[-1],
                               'paramdict':psets[parind],
                               'modelpath':modelpath}
                    workers.append(mp.Process(target=simulateGlobalStateSpace, args=(expargs, )))
                    
                for p in workers:
                    p.start()
                for p in workers:
                    p.join()                    

            for fname in fnames:
                df = pd.read_csv(DESTINATION + fname, index_col=0, dtype='str')
                for col in df.columns:
                    aggregator[fname.split('.')[0]][col] = df.loc[fname.split('.')[0], col]
            

        aggregator['cost'] = str(psets[parind].get('cost',0))
        to_file(outpath, aggregator, colorder)
        parind += 1        
            
    print('Done!')
    print('This took ' + str(time.time() - clock) + 's')
    
#####################################

def cleanfname(fname):
    specialchars = ['\\', '{', '}', ' ']
    for sc in specialchars:
        fname = fname.replace(sc,'')
    return(fname)

def tfStates(ss,READOUTS,cutoffs):
    state = {}
    for rd in READOUTS:
        if rd in ['Mig1','Dot6']:
            val = np.log10((ss[rd])/(1.001-ss[rd]))
        else:
            val =  ss[rd]

        if val > cutoffs[rd]:
            dec = 1
        else:
            dec = 0
        state[rd] = {'val':val, 'dec':dec}
    return(state)


def load_data(path):
    with open(path, 'r') as infile:
        expdat = yaml.safe_load(infile)
    return(expdat)

def main():
    parser = OptionParser()
    parser.add_option('-d','--debug',action='store_true',default=False,
                      help='will not calculate thresholds')
    parser.add_option('-m','--model-path',type='str',default='',dest='modelpath',
                      help='Path to model file')    
    parser.add_option('-p','--parameter-sets',type='str', default='',
                      help='Path to parameter set file.')
    opt,args = parser.parse_args()
    
    simulateStrains(opt.modelpath,
                    debug=opt.debug,
                    parameterSetPath=opt.parameter_sets)
    
if __name__ == '__main__':
    main()



