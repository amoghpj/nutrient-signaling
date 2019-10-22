import numpy as np
import time
import matplotlib
import sys
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nutrient_signaling.simulators import get_simulator
import os
import yaml

home = os.path.expanduser('~')
modelpath = home + '/jalihal_projects/Research/nutrient-signaling/data/2019-06-21'
readouts = ['Gis1','Mig1','Dot6','Gcn4','Rtg13','Gln3']

cutoffs = {'Gis1':0.5,
           'Gcn4':0.28,
           'Rtg13':0.3,
           'Gln3':0.5,
           'Mig1':1.325,
           'Dot6':1.85}

def tfStates(ss,readouts,cutoffs):
    state = {}
    for rd in readouts:
        if rd in ['Mig1']:#,'Dot6']:
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
    plt.savefig('img/' + whatsimulation + '.png')

def load_data(path):
    with open(path, 'r') as infile:
        #lines = infile.readlines()
        expdat = yaml.load(infile)
    return(expdat)

def do_experiment(experiment):
    model = get_simulator(modelpath=modelpath,
                       simulator='py')
    ismutant = False
    # Do mutant simulation
    if any(experiment['mutant']['pars']):
        mutpars = experiment['mutant']['pars']
        for k,v in mutpars.items():
            if type(v) is str and 'x' in v:
                mutpars[k] = float(v.replace('x',''))*model.model.pars[k]
        model.set_attr(pars = mutpars)
        ismutant = True
    if any(experiment['mutant']['ics']):
        model.set_attr(ics = experiment['mutant']['ics'])
        ismutant = True        
    
    prepars = experiment['preshift']
    postpars = experiment['postshift']    
    model.set_attr(pars=prepars, tdata=[0,90])
    model.set_attr(pars=prepars, tdata=[0,90], ics=model.get_ss())
    P = model.simulate_and_get_points()
    newss = model.get_ss()
    s = tfStates(newss, readouts, cutoffs)
    whatsimulation = str(experiment['id']) + '-' + experiment['strain']
    plotTfs(P, readouts, whatsimulation)
    return(s, whatsimulation)


def interpret_experiment(state, experiment, path):
    """
    TODO Fix Dot 6 interpretation
    """
    s = ''
    for rd in readouts:
        s += "%s\t%0.3f\t%s\n" %(rd, state[rd]['val'],state[rd]['dec'])
    if state['Dot6']['val'] < cutoffs['Dot6']\
       and state['Dot6']['val'] > 0.3*cutoffs['Dot6']:
        growthdec = 'Strain exhibits slow growth'
    elif state['Dot6']['val'] < 0.3*cutoffs['Dot6']:
        growthdec = 'Strain is inviable'            
    else:
        growthdec = 'wt growth/viable'
    report = "* %s-%s\n%s\nNutrient medium: %s\nObserved: %s\n"\
        "Predicted: %s\n[[./img/%s.png]]\n\n"\
        % (str(experiment['id']),
           experiment['strain'],
           s,
           experiment['nutrientCondition'],
           experiment['growth'],
           growthdec,
           path)    
    return(report)

def compare_model_predictions():
    datapath = 'data/experimental-data-full.yaml'
    expdat = load_data(datapath)
    report = ''
    start = time.time()
    for i, experiment in enumerate(expdat):
        print(experiment['id'])
        state, path = do_experiment(experiment)
        report += interpret_experiment(state, experiment, path)        
        # if i > 4:
        #     state, path = do_experiment(experiment)
        #     report += interpret_experiment(state, experiment, path)
        # if i > 10:
        #     break
    
    print('This took ' + str(time.time() - start) + 's')
    
    with open('report.org','w') as outfile:
        outfile.write(report)
    
compare_model_predictions()


