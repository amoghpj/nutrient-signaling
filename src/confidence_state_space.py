"""
Author: Amogh Jalihal
This script is used to generate Figure 4(c).
"""
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
font = {'family' : 'Sans',
        'size'   : 17}

matplotlib.rc('font', **font)

def compare_states(states, settings):
    '''
    states : list of strings, the predicted cellular state (boolean of 6 readouts)
             for every parameter set in df.
    '''
    # tuple for OFF vs ON confidence for each of the 6 TFs
    numTFs= len(settings.readouts)
    confidence = [{'off':0,'on':0} for _ in range(numTFs)]
    
    # States are of form 010010. Inner loop splits each character 'c'    
    splitstates = [[c for c in s] for s in states]
    total = float(len(states))

    ## For each prediction by a parameter set
    for s in splitstates:
        ## For each TF
        for i, c in enumerate(s):
            ## Record the number of times the TF is OFF or ON
            if c == '0':
                confidence[i]['off'] += 1.
            elif c == '1':
                confidence[i]['on'] += 1.
                
    ## Store as fraction
    for i in range(len(confidence)):
        confidence[i]['off'] = confidence[i]['off']/total
        confidence[i]['on'] = confidence[i]['on']/total
    return(confidence)


def draw_grey_background(odd, ax, ypos, nstates, width, height):
    bggrey = '#e6e6e6ff'#(0.95, 0.95, 0.95)
    if odd == 1:
        r = Rectangle((- width/2 - 0.1,ypos - height/2 - 0.3),
                      height=height + 0.4,
                      width = len(nstates),
                      facecolor=bggrey)
        ax.add_artist(r)
        odd = 0
    else:
        odd = 1
    return odd
        
def visualize(confidence, numpsets, evidencedict, settings,
              inputtype='tfstate', # 'dynamics'
              plotname='tf-predictions.pdf'):
    '''
    confidence: list of nested dictionaries with the following strucuture:
    [{
      'name':STRAIN-NAME, 
      'states':{
               NUTRIENT-CONDITION: {'off':FRACTION OFF, 'on':FRACTION ON},
               (7 more)
               },
    (16 more)
    }]
    There should be 17 strains, 8 nutrient conditions.
    This function can also handle predictions of dynamic responses. see dynamic-responses.py
    '''
    
    ## Color Settings. Leaving older colors in here just in case.
    robustgreen = '#d4efdf'#'k''#a9dfbf'##d4efdf'#
    robustred = '#f5b7b1' # '#ffffff'
    
    fragilered = 'r' # '#e74c3c'#
    fragilegreen = 'g'#'#27ae60'#'#196f3d'#'g'
    

    strains = [conf['name'] for conf in confidence]
    print(strains)
    nstates = settings.nutrientStates
    nested = {}

    y_positions = [i for i in range(len(strains))]
    x_positions = [i for i in range(len(nstates))]
    
    print(y_positions)
    
    plotwidth = 10
    plotheight = 10
    f, ax = plt.subplots(1,1,figsize=(plotwidth, plotheight))
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.xaxis.tick_top()

    ## NOTE I have hand tuned these values to get the bars to look nice
    bgheight = 0.5
    bgwidth = 0.8
    height = 0.1
    width = 0.05
    padding  = 0.95
    spacing = 0.06
    odd = 1
    for ypos, strain in enumerate(strains):
        # print(ypos)
        ## lay down the background grey bars
        odd = draw_grey_background(odd, ax, ypos, nstates, bgwidth, bgheight)
        for xpos, nstate in enumerate(nstates):
            readout_counter = 0
            for conf in confidence[ypos]['states'][nstate]:
                off = conf['off']
                on = conf['on']
                red = robustred
                green = robustgreen
                if inputtype == 'tfstates':
                    if min(on, off) > 0.1:
                        red = fragilered
                        green = fragilegreen
                elif inputtype == 'dynamics':
                    if on < 0.8:
                        red = fragilered
                        green = fragilegreen                        
                    
                ystart = ypos - 0.25
                yend_off = off/2.
                yend_on = on/2.                
                X = xpos + (readout_counter -3.)*(width + spacing)
                annotX = X + width/2.
                annotY = ystart - 0.2
                
                # First draw a rectangle that will serve as the border
                border = 0.02
                r = Rectangle((X- border/2. , ystart - border/2.),
                              height= 0.5 + border,
                              width = width + border,
                              facecolor='k')
                ax.add_artist(r)
                
                ## Draw the 'off' bar, in red
                r = Rectangle((X , ystart ),
                              height= yend_off,
                              width = width,
                              facecolor=red)
                ax.add_artist(r)

                ## Draw the 'on' bar, in green
                r = Rectangle((X , ystart + yend_off ),
                              height= yend_on,
                              width = width,
                              facecolor=green)
                ax.add_artist(r)

                # annotation for experiment agreement
                strain_nstate_readout = strain + '_' + nstate + '_' + settings.readouts[readout_counter]
                if strain_nstate_readout in evidencedict.keys():
                    gtstate = evidencedict[strain_nstate_readout]
                    color = 'r'
                    if gtstate == 'ON':
                        color = 'g'
                    predstate = 'ON'
                    ## NOTE Model prediction matches the experiment if the
                    ## majority prediction matches.
                    if off > on:
                        predstate = 'OFF'
                    if gtstate == predstate:
                        ax.plot(annotX,annotY, color+'o',ms=5.)
                    else:
                        ax.plot(annotX,annotY, color+'o',ms=5.,mfc='#ffffff')
                        
                readout_counter += 1
    ax.set_yticks(y_positions)
    strainnames = []
    for strain in strains:
        if strain == 'wt':
            strainnames.append('wt')
        else:
            strainnames.append(settings.mapper[strain])
    ax.set_yticklabels(strainnames)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(nstates)    
    #plt.xlabel('Gis1, Mig1, Dot6, Gcn4, Rtg13, Gln3')
    plt.xlim([min(x_positions) - width/2 - padding,max(x_positions) + width/2 + padding])
    plt.ylim([min(y_positions) - height/2 - 4*padding,max(y_positions) + height/2 + padding])    
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)        
    plt.tight_layout()
    # plt.savefig('state-space-comparison-robust-representative.pdf', dpi=500)
    plt.savefig(plotname, dpi=300)


class Settings():
    """
    Store some global settings that can be passed around 
    easily.
    1. self.readouts stores the TF names
    2. self.mapper maps mutant definitions in the input data file
       to interpretable, prettified LaTeX strings
    """
    def __init__(self):
        self.readouts = ['Gis1','Mig1','Dot6','Gcn4','Rtg13','Gln3']
        self.nutrientStates = ['HCHG','HCHN','HCHP','HCLN','LCHG','LCHN','LCHP','LCLN']
        self.mapper = {'Sch9':'$\Delta$sch9',
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
                       'wt':'wt'
        }        
        
def main():
    settings = Settings()
    
    ## State Predictions:
    ## This tab-separated file is generated using 
    ## src/compare-state-space-predictions.py
    predictionpath = 'output/global_space_24066.csv'
    predictiondf = pd.read_csv(predictionpath,
                                sep='\t', index_col=False,
                                dtype='str')
    ## Filter by cost, retain anything less than twice the minimum
    predictiondf = predictiondf.astype({'cost':'float'})
    cmin = predictiondf['cost'].min()
    predictiondf = predictiondf[predictiondf['cost'] <= 2.*cmin]

    # Curated evidence for states:
    ## This comma-separated file is manually curated,
    ## with references to the original publication
    evidencepath = 'data/csv/tf-state-experimental-evidence.csv'
    evidencedf = pd.read_csv(evidencepath, index_col=None)
    evidencedict = {row['strain']:row['state']                     # Create dict {strain:state}
                       for ind, row in                                # Loop over key, value of 
                       evidencedf.dropna(how='any').T.to_dict().items()} # evidence table converted to dict
    # Put wt on top
    colorder = [c for c in reversed(list(predictiondf.columns))]

    colorder.remove('cost')

    seen = set()
    allstrains = []
    for c in colorder:
        strain = c.split('_')[0]
        if strain not in seen:
            allstrains.append(strain)
            seen.add(strain)
            

    # make a nested dictionary to hold predictions
    confidence = [{'name':strain,
                   'states':{nstate:[] for nstate in
                             settings.nutrientStates}}
                  for strain in allstrains]
    strainMapper = {conf['name']:i for i, conf in enumerate(confidence)}
    print(strainMapper)
    # Record predictions
    ## Each entry in predictiondf is a string of 0s and 1s
    ## which records the state of the TFs in the order given
    ## by settings.readouts.  The columns in this file take
    ## the form STRAIN-NAME_NUTRIENT-CONDITION.
    ## The function compare_states() reads the list containing
    ## predictions from each parameter set, and for each TF,
    ## counts the number of time the parameter sets predicted
    ## an ON vs an OFF. This is reported as a fraction of the total
    ## number of predictions, and is returned as a list. (See docs for visualize())
    
    for col in colorder:
        strain, nstate = col.split('_')
        strainid = strainMapper[strain]
        confidence[strainid]['states'][nstate] = compare_states(list(predictiondf[col]), settings)
    print(colorder)
    numpsets = predictiondf.shape[0]

    # Plot the confidences as a tableu
    visualize(confidence, numpsets, evidencedict, settings)        

if __name__ == '__main__':
    main()
