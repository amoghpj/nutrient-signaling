import os
import sys
import yaml
import time
import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path
import multiprocessing as mp
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optparse import OptionParser
# Local Imports
from nutrient_signaling.simulators import get_simulator

home = os.path.expanduser('~')
READOUTS = ['Gis1','Mig1','Dot6','Gcn4','Rtg13','Gln3']
SIMULATOR = 'cpp'
SIMARGS = {'executable':'robust.o'}

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
            dec = 'On'
        else:
            dec = 'Off'
        state[rd] = {'val':val, 'dec':dec}
    return(state)

def plotTfs(P, READOUTS, whatsimulation, cutoffs, fname):
    plt.close()
    f = plt.figure(figsize=(6,4))
    axmap = {'Gln3':0,'Rtg13':1,'Mig1':2,'Gcn4':3,'Gis1':4,'Dot6':5}
    for i, rd in enumerate(READOUTS):
        ax = f.add_subplot(2,3,axmap[rd] + 1)
        if rd == 'Mig1' or rd == 'Dot6':
            ax.plot(P['t'], [np.log10((p)/(1-p)) for p in P[rd]], label=rd)
            if rd == 'Mig1':
                ax.set_ylim([0.5,1.6])
            elif rd == 'Dot6':
                ax.set_ylim([0, 2.8])
        else:
            ax.plot(P['t'], P[rd], label=rd)
            ax.set_ylim([0.0,1.0])
        ax.legend(fancybox=False, framealpha=0.0)
        ax.axhline(cutoffs[rd],c='k',ls='--',lw=1.0)
    plt.suptitle(whatsimulation)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])    
    plt.savefig('img/' + fname + '.png', dpi=200)

def load_data(path):
    with open(path, 'r') as infile:
        expdat = yaml.safe_load(infile)
    return(expdat)

def do_experiment(expargs):
    """
    Carry out an in silico shift experiment with the given experiment
    specification.
    PROG: incomplete
    """
    experiment = expargs['experiment']
    shareddict = expargs['shareddict']
    cutoffs = expargs['cutoffs']
    paramdict = expargs['paramdict']
    if paramdict is None:
        paramdict = {}
        
    modelpath = expargs['modelpath']    
    uid = str(experiment['id']) + '-' + experiment['strain']
    fname = str(uid)
    if SIMULATOR == 'cpp':
        SIMARGS['simfilename'] = cleanfname(uid) + '.dat'
    # Initialize model object
    model = get_simulator(modelpath=modelpath,
                          simulator=SIMULATOR,
                          **SIMARGS)
    model.set_attr(pars=paramdict)
    
    ## Do mutant simulation
    ## 1. If parameters are defined, set parameters
    ## 2. If initial conditions are specified, set initial conditions
    ## Else, leave at default, resulting in wt simulation
    if experiment['mutant']['pars'] is not None:
        mutpars = experiment['mutant']['pars']
        for k,v in mutpars.items():
            ## Simple check to implement a "X times" type of specification
            ## Update mutpars inplace
            if type(v) is str and 'x' in v:
                mutpars[k] = float(v.replace('x',''))*model.pars[k]
        model.set_attr(pars = mutpars)
    if experiment['mutant']['ics'] is not None:
        model.set_attr(ics = experiment['mutant']['ics'])
    
    prepars = experiment['preshift']
    if prepars is None:
        prepars = {}
    postpars = experiment['postshift']
    if postpars is None:
        postpars = {}
    postics = experiment['postshiftics']
    if postics is None:
        postics = {}    
    model.set_attr(pars=prepars, tdata=[0,90])
    newics = model.get_ss()
    newics.update(postics)
    model.set_attr(pars=postpars, tdata=[0,90], ics=model.get_ss())
    
    P = model.simulate_and_get_points()
    newss = model.get_ss()
    state = tfStates(newss, READOUTS, cutoffs)

    for rep in [' ', '\\']:
        if rep in fname:
            fname = fname.replace(rep,'')
    #plotTfs(P, READOUTS, uid, cutoffs, fname)
    interpret_experiment(uid, state, expargs, fname)    

def interpret_experiment(uid, state, expargs, fname):
    """
    Interpret each in silico simulation.
    Currently records the global predicted cell state,
    model-experiment match, and predicted viability of the state.
    """
    # TODO: Maybe classes will be more useful here?
    experiment = expargs['experiment']
    shareddict = expargs['shareddict']
    cutoffs = expargs['cutoffs']    

    ## The experiment specification can carry the field
    ## =phenotypeInterpreted= which is how I represent
    ## the strain's phenotype in terms of the model READOUTS.
    ## If this field is present, use this to interpret the
    ## simulation results
    phenotypeMatches = []
    if experiment['phenotypeInterpreted'] is not None:
        for k, v in experiment['phenotypeInterpreted'].items():
            if experiment['phenotypeInterpreted'][k].strip() == state[k]['dec'].strip():
                phenotypeMatches.append(True)
            else:
                phenotypeMatches.append(False)                
                # Simple string comparison to decide match.
                # Takes care of potential white spaces, but not
                # case of text. Printing below to help debug.
            # print(uid + '\t' + k +'\t' + experiment['phenotypeInterpreted'][k]\
            #       + '\t' + state[k]['dec'])
    else:
        ## Not desirable, but fall back on the full list of READOUTS.
        ## This is quite unreliable.
        for r in READOUTS:
            if experiment['expected'][k] == state[k]['dec'].strip():
                phenotypeMatches.append(True)
            else:
                phenotypeMatches.append(False)                
                # Simple string comparison to decide match.
                # Takes care of potential white spaces, but not
                # case of text. Printing below to help debug.
                print(uid + '\t' + k +'\t' + experiment['phenotypeInterpreted'][k]\
                      + '\t' + state[k]['dec'])            
            
    ## Next, irrespective of the strain, we can guesstimate the state
    ## the cell _should be in_ in order to grow in a given nutrient
    ## environment, in order to infer viability. This typically means
    ## a particular transcription factor should be expressed in a
    ## given nutrient condition. The field =forGrowth= records
    ## this(ese) state(s).
    viability = []
    if experiment['forGrowth'] is not None:
        for tf, tfstate in experiment['forGrowth'].items():
            if tfstate == state[tf]['dec']:
                viability.append(True)
            else:
                viability.append(False)
    
    ## Store the global TF state. Useful for analysis in the future
    s = '|*TF*|*Interpreted*|*Simulated*|*Simulation*|\n'

    for rd in READOUTS:
        # Format the value of each readout in individual rows
        emph = '/'
        expected = '-'
        simulated = '-'
        if rd in experiment['phenotypeInterpreted']:
            emph = '*'
            expected = experiment['phenotypeInterpreted'][rd]
            simulated = state[rd]['dec']
        s += "|%s|%s|%s|%0.3f|\n" %(emph + rd + emph,
                                    expected,
                                    simulated,
                                    state[rd]['val'])
    strainGrowthStatus = "growth"
    for v in viability:
        if not v:
            strainGrowthStatus = "no growth"
            break
        
    simulationAgreementStatus = "agrees"
    for p in phenotypeMatches:
        if not p:
            simulationAgreementStatus = "does not agree"

    modelMatchesFlag = False
    if simulationAgreementStatus == "agrees":
        shareddict['counter'] += 1
        #shareddict[uid]['counter'] +=1
        shareddict[uid] += 1        
        modelMatchesFlag = True
        
    template = "* {}\n"\
        ":PROPERTIES:\n:CUSTOM_ID: sec:{}\n:END:\n"\
        "{}\n\n"\
        "*Description*: {} studied a /{}/ strain ({}) grown in {}.\n\n"\
        "*Representation*:\n\n/Preshift Parameters/\n{}\n\n/Postshift Parameters/\n{}\n\n/Postshift Initial Conditions/\n{}\n\n"\
        "/Mutant/\n\nParameters:\n{}\n\nInitial Conditions:\n{}\n\n"\
        "*Growth*: In this medium, the strain showed {}. The model predicts that the strain will show {}.\n\n"\
        "*Model {} with experiment*.\n\n"\
        "#+ATTR_LATEX: :height 0.25\\textheight\n"\
        "[[./img/{}.png]]\n\n"

    report = template.format(uid,
                             fname,                             
                             s,
                             experiment['shortname'],
                             experiment['strain'],
                             experiment['background'],
                             experiment['nutrientCondition'],
                             stringify(experiment['preshift']),
                             stringify(experiment['postshift']),
                             stringify(experiment['postshiftics']),                             
                             stringify(experiment['mutant']['pars']),
                             stringify(experiment['mutant']['ics']),                             
                             experiment['phenotypeReported'],
                             strainGrowthStatus,
                             simulationAgreementStatus,
                             fname)

    #shareddict[uid]['report'] += report
    #print(shareddict)    

def stringify(datadict):
    if datadict is None:
        return('')
    else:
        return('\n'.join(['|' + str(k) + '|' + str(v)+'|' for k,v in datadict.items()]))
    
def define_state_space(datapath, modelpath):
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


def compare_model_predictions(experimentpath,
                              modelpath,
                              debug=False,
                              parameterSetPath=None):
    print('Initializing...')
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
    if not debug:
        sd, cutoffs = define_state_space(experimentpath, modelpath)
    else:
        cutoffs = {'Gis1':0.498,
                   'Mig1':1.344,
                   'Dot6':1.058,
                   'Gcn4':0.489,
                   'Rtg13':0.46,
                   'Gln3':0.564}    
    print('Loading data...')
    expdat = load_data(experimentpath)
    report = ''
    clock = time.time()
    kl = []
    start = 0
    end = len(expdat)
    print('Starting experiments...')
    aggregator = {'cost':[]}
    uidlist = []
    for i, experiment in enumerate(expdat[start:min(end, len(expdat))]):
        uid = str(experiment['id']) + '-' + experiment['strain']
        uidlist.append(uid)
        aggregator[uid] = []

    metacounter = {'counter':[], 'cost':[]}
    #print('uid\ttf\texpec\tsimul')
    parind = 0
    
    while parind < psetlength:
        manager = mp.Manager()
        shareddict = manager.dict()
        #shareddict.update({uid:{'report':'','counter':0} for uid in uidlist})
        for uid in uidlist:
            shareddict[uid] = 0
        shareddict['counter'] = 0

        if debug:
            for i, experiment in enumerate(expdat[start:min(end, len(expdat))]):
                expargs = {'experiment':experiment,
                           'shareddict': shareddict,
                           'cutoffs': cutoffs,
                           'paramdict':psets[parind],
                           'modelpath':modelpath}
                do_experiment(expargs)
        
        else:
            with mp.Pool() as pool:
                jobs = []
                for i, experiment in enumerate(expdat[start:min(end, len(expdat))]):
                    expargs = {'experiment':experiment,
                               'shareddict': shareddict,
                               'cutoffs': cutoffs,
                               'paramdict':psets[parind],
                               'modelpath':modelpath}
                    job = pool.apply_async(do_experiment, args=(expargs, ))
                    jobs.append(job)
                for job in jobs:
                    job.wait()

        aggregator['cost'].append(psets[parind].get('cost',0))
        counter = 0
        cumsum = 0
        for uid in uidlist:
            cumsum += shareddict[uid]
            if shareddict[uid] == 1:
                aggregator[uid].append(1)
                counter += 1
            else:
                aggregator[uid].append(0)                
        #metacounter['counter'].append(shareddict['counter'])
        #print(parind,'\t', shareddict['counter'], aggregator['cost'][-1])
        print(parind,'\t', cumsum, '\t', shareddict['counter'],'\t',counter, '\t', aggregator['cost'][-1])
        
        del manager
        del shareddict
        parind += 1        
            
    print('Done!')
    print('This took ' + str(time.time() - clock) + 's')
    
    # Write report to file
    if psetlength == 1 and psets[0] is None:
        preamble = '#+OPTIONS: toc:nil\n\n#+LATEX_HEADER: \\usepackage[margin=0.5in]{geometry}\n\n'
        # summary  = '{} out of {} experiments were correctly predicted\n\n'.format(shareddict['counter'],min(end, len(expdat))-start)
        # with open('report.org','w') as outfile:
        #     report = preamble + summary + ''.join([shareddict[k] for k in kl])
        #     outfile.write(report)
    else:
        agdf = pd.DataFrame(aggregator)
        agdf.to_csv('summary_' + str(len(psets.keys())) + '_test.csv')


def main():
    parser = OptionParser()
    parser.add_option('-d','--debug',action='store_true',default=False,
                      help='will not calculate thresholds')
    parser.add_option('-m','--model-path',type='str',default='',dest='modelpath',
                      help='Path to model file')    
    parser.add_option('-p','--parameter-sets',type='str', default='',
                      help='Path to parameter set file.')
    parser.add_option('-e','--experiments',default='',type='str',
                      help='Path to yaml file containing experimental results')
    opt,args = parser.parse_args()
    
    compare_model_predictions(opt.experiments,
                              opt.modelpath,
                              debug=opt.debug,
                              parameterSetPath=opt.parameter_sets)
    
if __name__ == '__main__':
    main()



