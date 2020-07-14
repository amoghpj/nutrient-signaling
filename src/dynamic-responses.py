import nutrient_signaling.simulators as sim
import pandas as pd
import sys
import importlib
utils = importlib.import_module("compare_state_space_predictions")
import matplotlib.pyplot as plt
import matplotlib
from tqdm import tqdm
font = {'family' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)
SIMARGS = {'executable': 'dynamics.o'}
SIMULATOR = 'cpp'
mapper = {'Sch9':'$\Delta$sch9',
          'TORC1':'$\Delta$tor1',
          'Snf1':'$\Delta$snf1',
          'EGO':'$\Delta$gtr1/2',
          'PDE':'$\Delta$pde1/2',
          'EGOGAP':'$\Delta$lst4/7',
          'Ras':'$\Delta$ras2',
          'Sak':'$\Delta$sak1',
          'Gcn2':'$\Delta$gcn2',
          'Cyr1':'$\Delta$cyr1',
          'PKA':'$\Delta$tpk1/2/3',
          'Sch9-inhibit-gcn2':'GCN2-S557',
          'Snf1-activates-nGln3':'GLN3 $\Delta$ST',
          'TORC1-inhibits-Gln3':'GLN3 $\Delta$TT',
          'Sch9-inhibits-PKA':'$\Delta$bcy1',
          'PKA-inhibits-Ras':'$\Delta$ira1/2',
          'wt':'wt'}

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
    
def plot_dynamic_responses():
    """
    For each variable in defined in MAJORVARIABLES,
    make a tableu of time courses where columns are
    nutrient conditions and rows are mutants
    Preshift will always be starvation, except in the case of
    LCLN
    """
    modelpath = "data/2018-9-26-12-3-no-sigma"
    majorvariables = [
        "Snf1",
        "Sch9",
        "PKA",
        "cAMP"
    ]
    psets = read_k_alternate_params(k=100)
    
    nutrientStates = utils.createNutrientInputSpace()

    # order of nutrient shifts
    norder = ['HCHG','HCHN','HCHP', 'HCLN',
              'LCHG','LCHN','LCHP', 'LCLN']
    
    mutantSpecs = utils.getMutantSpecs(modelpath)

    # Main loop, over the variables to plot
    for mv in majorvariables:
        print(mv)
        plt.close()
        f = plt.figure(figsize=(20,30))
        imgid = 1
        # loop over the strains
        for mut in tqdm(mutantSpecs):
            if mv == 'cAMP':
                tmax = 4.0
            else:
                tmax = 90.0
                
            # loop over the nutrient conditions
            for no in norder:
                # silly loop to get the correct specification 
                ns = [n for n in nutrientStates if n['name'] == no][0]
                ax = f.add_subplot(len(mutantSpecs),len(nutrientStates),imgid)
                # loop over parameter sets
                for i, pset in psets.iterrows():
                    # Initialize model
                    model = sim.get_simulator(modelpath = modelpath,
                                              simulator="cpp",
                                              executable=SIMARGS['executable'])
                    # Set preshift parameters
                    if ns['name'] == 'LCLN':
                        preshift = {'Carbon':1.0,"Glutamine_ext":1.0,"ATP":1.0}
                    else:
                        preshift = {'Carbon':0.0,"Glutamine_ext":1.0,"ATP":0.0}
                    # initialize with pset
                    deletion_par = pset.to_dict()
                    deletion_par.update(mut['pars'])
                    deletion_par.update(preshift)
                    deletion_ics = mut['ics']
                    model.set_attr(pars=deletion_par,
                                   ics=deletion_ics,
                                   tdata = [0,90])
                    # Reinitialize.
                    # NOTE This can be cleaner
                    deletion_par = pset.to_dict()
                    deletion_par.update(mut['pars'])
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
                    ax.plot(points['t'],points[mv],'k',alpha=0.4)
                    # End loop
                # Configure plot
                if imgid <= len(nutrientStates):
                    ax.set_title(ns['name'])

                if (imgid-1) % len(nutrientStates) == 0:
                    ax.set_ylabel(mapper[mut['name']])
                if mv != 'cAMP':
                    ax.set_ylim([0,1.1])
                else:
                    ax.set_ylim([0,3.0])
                ax.set_xticklabels([])
                ax.set_yticklabels([])                
                imgid += 1
        # plt.suptitle(mv)
        plt.tight_layout()
        plt.savefig('img/dynamic-responses-' + mv.lower() + '-confidence.png', dpi=100)
        
if __name__ == "__main__":
    plot_dynamic_responses()



