import nutrient_signaling.simulators as sim
import pandas as pd
import os
import sys
import importlib
utils = importlib.import_module("compare_state_space_predictions")
import confidence_state_space as css
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from tqdm import tqdm
from optparse import OptionParser

MODELPATH = "data/2018-9-26-12-3-no-sigma"
SIMARGS = {'executable': 'dynamics-1.o'}
SIMULATOR = 'cpp'
MAJORVARIABLES = [
    "Snf1",
    "Sch9",
    "PKA",
    "cAMP",
    "Gcn2",
    "Gis1",
    "Gln3",
    "Rtg13"]


def read_k_alternate_params(k=10):
    """
    Read file containing alternate parameter sets and sample /k/ sets from top 5000
    """
    path = './data/parameter-set-ensembles/expansion_iter4.csv'
    df = pd.read_csv(path)
    # Same as in TimeCourse.setNumberOfParameterSets()
    ParameterSetsToUse = df.sort_values(by='cost')[0:5000]
    ParameterSetsToUse = ParameterSetsToUse.sample(k)            
    return(ParameterSetsToUse)

def compute_mse(referencetraj, newtraj, is_camp=False):
    """
    Take two lists storing time coureses. Compute and return SSE 
    """
    if is_camp:
        mse = sum([np.power((r - n), 2.) for r,n in zip(referencetraj[0:5000],newtraj[0:5000])])/5000
    else:
        mse = sum([np.power((r - n), 2.) for r,n in zip(referencetraj,newtraj)])/len(newtraj)
    return mse

def simulate_strain_and_nutrient_condition(ns, pset={},mut={}):
    """
    Simulates a single strain using a parameter set if specified.
    Requires a nutrient condition specification
    """
    model = sim.get_simulator(modelpath = MODELPATH,
                              simulator="cpp",
                              executable=SIMARGS['executable'])
    tmax = 90
    # Set preshift parameters
    if ns['name'] == 'LCLN':
        preshift = {'Carbon':1.0,"Glutamine_ext":1.0,"ATP":1.0}
    else:
        preshift = {'Carbon':0.0,"Glutamine_ext":1.0,"ATP":0.0}
    # initialize with pset
    deletion_par = {}
    if type(pset) != dict:
        deletion_par.update(pset.to_dict())
    deletion_par.update(mut.get('pars', {}))
    deletion_par.update(preshift)
    deletion_ics = mut.get('ics', {})
    
    model.set_attr(pars=deletion_par,
                   ics=deletion_ics,
                   tdata = [0,90])
    # Reinitialize the parameter specification
    # NOTE This can be cleaner
    deletion_par = {}
    if type(pset) != dict:
        deletion_par.update(pset.to_dict())
    deletion_par.update(mut.get('pars', {}))
    postshift = {'Carbon':ns['C'],
                 "ATP":ns['C'],
                 "Glutamine_ext":ns['Gln'],
                 "NH4":ns['NH4'],
                 "Proline":ns['Pro']}

    deletion_par.update(postshift)
    model.set_attr(ics=model.simulate_and_get_ss(),
                   pars=deletion_par,
                   tdata=[0,tmax])

    points = model.simulate_and_get_points()
    return(points)

def dynamic_responses(settings,numpsets=10):
    """
    For each variable in defined in MAJORVARIABLES,
    make a tableu of time courses where columns are
    nutrient conditions and rows are mutants
    Preshift will always be starvation, except in the case of
    LCLN
    """

    psets = read_k_alternate_params(k=numpsets)
    
    nutrientStates = utils.createNutrientInputSpace()

    # order of nutrient shifts
    norder = ['HCHG','HCHN','HCHP', 'HCLN',
              'LCHG','LCHN','LCHP', 'LCLN']
    
    mutantSpecs = utils.getMutantSpecs(MODELPATH)
    referencetraj = {no:{mv:None for mv in MAJORVARIABLES}
                          for no in norder }
    for ns in nutrientStates:
        points = simulate_strain_and_nutrient_condition(ns)
        for mv in MAJORVARIABLES:
            referencetraj[ns['name']][mv] = points[mv]

    store = {mut['name']:{no:{mv:[] for mv in MAJORVARIABLES}
                          for no in norder }
             for mut in mutantSpecs}    
    # loop over the strains
    for mut in tqdm(mutantSpecs):
        # loop over the nutrient conditions
        for no in norder:
            # silly loop to get the correct specification 
            ns = [n for n in nutrientStates if n['name'] == no][0]
            # loop over parameter sets
            for i, pset in psets.iterrows():
                # Initialize model
                points = simulate_strain_and_nutrient_condition(ns, mut=mut, pset=pset)
                # Main loop, over the variables to plot
                for mv in MAJORVARIABLES:
                    if mv == 'cAMP':
                        is_camp = True
                    else:
                        is_camp = False
                    # Compute MSE with reference pset
                    v = compute_mse(referencetraj[ns['name']][mv],
                                    points[mv], is_camp=is_camp)                        
                    store[mut['name']][no][mv].append(v)
                    
    with open('sse_values_'+str(numpsets)+'.csv','w') as outfile:
        for mut in mutantSpecs:
            for no in norder:
                for mv in MAJORVARIABLES:
                    l = [settings.mapper[mut['name']],no,mv]
                    l.extend([str(s) for s in store[mut['name']][no][mv]])
                    outfile.write(','.join(l) + '\n')

def process_results(numpsets, settings):
    path = 'sse_values_' + str(numpsets) + '.csv'
    if not os.path.exists(path):
        print('please run dynamic_responses() by setting do_simulations = True.')
        sys.exit()
    df = pd.read_csv(path, header=None, index_col=False)

    # # NOTE Calculating a cutoff to consider a timecourse to be a good fit
    # # There are two ways of thinking about this
    # # 1. We can calculate the median MSE across all simulations.
    # #    This assumes that a majority of the MSEs are "good"
    # vals = df[[i for i in range(3,numpsets+3)]].values.flatten()
    # cutoff = np.median(vals)
    # # 2. Select a maximum MSE based on visual examination of time course
    vals = df[(df[0]=='wt')
              & (df[1]== 'HCHN')
              & (df[2] == 'Sch9')][[i for i in range(3,numpsets+3)]].values.flatten()
    cutoff = np.max(vals)
    
    revmapper = {v:k for k,v in settings.mapper.items()}
    # Adapted from confidence_state_space.py
    # make a nested dictionary to hold predictions
    confidencedf = pd.DataFrame(columns=['strain','variable','nutrientState','robustFraction'],
                              index=pd.Index([i for i in range(df.shape[1])]))
    confidence = [{'name':revmapper[strain],
                   'states':{nstate:[] for nstate in
                             settings.nutrientStates}}
                  for straindesc, strain in settings.mapper.items()]    

          
    for i, row in df.iterrows():
        strain = row[0]
        nutrientState = row[1]
        variable = row[2]
        vals = list(row[3:numpsets+3])
        robustFraction = float(sum([1 for i in range(len(vals))
                                    if vals[i] < cutoff]))/float(numpsets)
        confidencedf.loc[i] = {'strain':strain,
                           'variable':variable,
                           'nutrientState':nutrientState,
                           'robustFraction':robustFraction}
    strains = df[0].unique()
    for strain in strains:
        for curind, conf in enumerate(confidence):
            if conf['name'] == revmapper[strain]:
                break
        for nstate in settings.nutrientStates:
            strain_nut_dyn = confidencedf[(confidencedf.strain == strain) &
                                          (confidencedf.nutrientState == nstate)]
            for mv in MAJORVARIABLES[0:4]:
                # This is a problem. 'on' and 'off' don't carry the same meanings as the TF plot
                # 'on' is fraction that is robust, 'off' is fraction that is not
                # So if 'on' > 0.9 then the simulations are robust
                rob = strain_nut_dyn[strain_nut_dyn.variable == mv].robustFraction.values[0]
                #print(rob)
                confidence[curind]['states'][nstate].append({'off':1.0 - rob, 'on':rob})
    evidencedict = {}
    css.visualize(confidence, numpsets, evidencedict, settings,
                  plotname='dynamics-'+ str(numpsets)+'.pdf',
                  inputtype='dynamics')
        
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-s','--simulate',action='store_true',default=False,dest='do_simulations',
                      help='Do simulations')
    parser.add_option('-n','--numpsets',type='int',default=None,dest='numpsets',
                      help='Number of parameter sets to use')
    opt,args = parser.parse_args()
    
    numpsets = opt.numpsets
    do_simulations = opt.do_simulations
    
    settings = css.Settings()

    if do_simulations:
        dynamic_responses(settings, numpsets=numpsets)
    process_results(numpsets, settings)

