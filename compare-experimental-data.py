import os
import yaml
import time
import numpy as np
import multiprocessing as mp
import matplotlib
import sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nutrient_signaling.simulators import get_simulator

home = os.path.expanduser('~')
modelpath = home + '/jalihal_projects/Research/nutrient-signaling/data/2019-06-21'
readouts = ['Gis1','Mig1','Dot6','Gcn4','Rtg13','Gln3']

# cutoffs = {'Gis1':0.5,
#            'Gcn4':0.28,
#            'Rtg13':0.3,
#            'Gln3':0.5,
#            'Mig1':1.325,
#            'Dot6':1.85}

def tfStates(ss,readouts,cutoffs):
    state = {}
    for rd in readouts:
        if rd in ['Mig1','Dot6']:
            val = np.log10((ss[rd])/(1-ss[rd]))
        else:
            val =  ss[rd]

        if val > cutoffs[rd]:
            dec = 'On'
        else:
            dec = 'Off'
        state[rd] = {'val':val, 'dec':dec}
    return(state)

def plotTfs(P, readouts, whatsimulation):
    plt.close()
    f,ax = plt.subplots(1,3,figsize=(6,2))
    for rd in readouts:
        if rd == 'Mig1':
            ax[1].plot(P['t'], [np.log10((p)/(1-p)) for p in P[rd]], label=rd)
        elif rd == 'Dot6':
            ax[2].plot(P['t'], [np.log10((p)/(1-p)) for p in P[rd]], label=rd)
        else:
            ax[0].plot(P['t'], P[rd], label=rd)
    ax[1].set_ylim([1.1,1.6])
    ax[2].set_ylim([1.2,2.8])
    ax[0].set_ylim([0.0,1.0])
    ax[0].legend(fancybox=False, framealpha=0.0)
    ax[1].legend(fancybox=False, framealpha=0.0)
    ax[2].legend(fancybox=False, framealpha=0.0)
    plt.tight_layout()
    plt.suptitle(whatsimulation)
    plt.savefig('img/' + whatsimulation.replace(' ','') + '.png')

def load_data(path):
    with open(path, 'r') as infile:
        #lines = infile.readlines()
        expdat = yaml.safe_load(infile)
    return(expdat)

def do_experiment(expargs):
    experiment = expargs['experiment']
    shareddict = expargs['shareddict']
    cutoffs = expargs['cutoffs']
    model = get_simulator(modelpath=modelpath,
                       simulator='py')
    ismutant = False
    # Do mutant simulation
    if experiment['mutant']['pars'] is not None:
        mutpars = experiment['mutant']['pars']
        for k,v in mutpars.items():
            if type(v) is str and 'x' in v:
                mutpars[k] = float(v.replace('x',''))*model.model.pars[k]
        model.set_attr(pars = mutpars)
        ismutant = True
    if experiment['mutant']['ics'] is not None:
        model.set_attr(ics = experiment['mutant']['ics'])
        ismutant = True        
    
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
    state = tfStates(newss, readouts, cutoffs)
    uid = str(experiment['id']) + '-' + experiment['strain']        
    plotTfs(P, readouts, uid)
    interpret_experiment(uid, state, expargs)    

def interpret_experiment(uid, state, expargs):
    """
    TODO: Fix Dot 6 interpretation
    """
    experiment = expargs['experiment']
    shareddict = expargs['shareddict']
    cutoffs = expargs['cutoffs']
    s = '|TF|Expected|Simulation | Boolean|\n'
    matches = []
    for rd in readouts:
        s += "|%s|%s|%0.3f|%s|\n" %(rd, experiment['expected'][rd], state[rd]['val'],state[rd]['dec'])
        if experiment['expected'][rd] == state[rd]['dec']:
            matches.append(rd)
            
    if state['Dot6']['val'] < cutoffs['Dot6']:
        growthdec = 'Strain exhibits growth in this condition'
    elif state['Dot6']['val'] > cutoffs['Dot6'] and state['Dot6']['val'] < 1.3*cutoffs['Dot6']:
        growthdec = 'Strain exhibits slow growth'            
    else:
        growthdec = 'Strain does not grow in this condition'
    template = "* {}\n{}\n\n{} of 6 tf states are identical.\n\n*Experiment*: {} studied a /{}/ strain ({}) grown in {}. {} ({}). The strain is {}.\n\nPredicted: {}.\n\n[[./img/{}.png]]\n\n"
    report = template.format(uid,
                             s,
                             len(matches),
                             experiment['shortname'],
                             experiment['strain'],
                             experiment['background'],
                             experiment['nutrientCondition'],
                             experiment['phenotypeReported'],
                             experiment['expReadout'],
                             experiment['growth'],
                             growthdec,
                             uid.replace(' ',''))
    if len(matches) >=5 :
        shareddict['counter'] += 1
    shareddict[uid] += report

def define_state_space(datapath):
    statedict = {r:[] for r in readouts}
    
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
                       simulator='py')
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
    
def compare_model_predictions():
    datapath = 'data/experimental-data-full.yaml'
    print('Initializing...')
    sd, cutoffs = define_state_space(datapath)
    print('Loading data...')
    expdat = load_data(datapath)
    report = ''
    clock = time.time()
    manager = mp.Manager()
    shareddict = manager.dict()
    shareddict['counter'] = 0
    kl = []
    start = 0
    end = 62
    print('Starting experiments...')
    for i, experiment in enumerate(expdat[start:end]):
        uid = str(experiment['id']) + '-' + experiment['strain']
        shareddict[uid] = ''
        kl.append(uid)
        
    with mp.Pool() as pool:
        jobs = []
        for i, experiment in enumerate(expdat[start:end]):
            expargs = {'experiment':experiment,
                       'shareddict': shareddict,
                       'cutoffs': cutoffs}
            job = pool.apply_async(do_experiment, args=(expargs, ))
            jobs.append(job)
        for job in jobs:
            job.wait()
            
    print('Done!')
    print('This took ' + str(time.time() - clock) + 's')
    preamble = '#+OPTIONS: toc:nil\n\n#+LATEX_HEADER: \\usepackage[bottom=0.5in,margin=1in]{geometry}\n\n'
    summary  = '{} out of {} experiments were correctly predicted\n\n'.format(shareddict['counter'],end-start)
    with open('report.org','w') as outfile:
        report = preamble + summary + ''.join([shareddict[k] for k in kl])
        outfile.write(report)
    
compare_model_predictions()


