import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patheffects as pe
font = {'family' : 'normal',
        'size'   : 17}

matplotlib.rc('font', **font)

READOUTS = ['Gis1','Mig1','Dot6','Gcn4','Rtg13','Gln3']

predictionpath = 'output/global_space_24066.csv'
groundtruthpath = 'data/ground_truth.csv'

df = pd.read_csv(predictionpath, sep='\t', index_col=False, dtype='str')
df = df.astype({'cost':'float'})
cmin = df['cost'].min()
df = df[df['cost'] <= 2.*cmin]
groundtruthdict = {row['strain']:row['state'] for ind, row in pd.read_csv(groundtruthpath, index_col=None).dropna(how='any').T.to_dict().items()}
print(groundtruthdict)

def compare_states(states):
    '''
    states : list of strings, the predicted cellular state (boolean of 6 readouts)
             for every parameter set in df.
    '''
    # tuple for OFF vs ON confidence for each position
    confidence = [[0,0] for _ in range(6)]
    
    # States are of form 010010. Inner loop splits each character 'c'    
    splitstates = [[c for c in s] for s in states]
    total = float(len(states))
    for s in splitstates:
        for i, c in enumerate(s):
            if c == '0':
                confidence[i][0] += 1.
            else:
                confidence[i][1] += 1.
    for i in range(len(confidence)):
        confidence[i][0] = confidence[i][0]/total
        confidence[i][1] = confidence[i][1]/total
    return(confidence)

def visualize(confidence, size, groundtruthdict):
    '''
    confidence: list of lists of size ((numstrains x numnutrConditions = 114?)
        x 2). First item is strain_nutrCondition, second item is list of
        lists of size (numReadouts (6) x 2)

    '''
    lightgreen = '#d4efdf'#'#a9dfbf'
    lightred = '#f5b7b1'
    darkred = '#e74c3c'#'r'
    darkgreen = '#27ae60'#'#196f3d'#'g'
    
    strains = []
    nstates = ['HCHG','HCHN','HCHP','HCLN','LCHG','LCHN','LCHP','LCLN']
    nested = {}
    # Retain non null
    for confitems in confidence:
        straindef, conf = confitems
        strain, nstate = straindef.split('_')
        if strain not in strains:
            strains.append(strain)
        if strains[-1] in nested.keys():
            nested[strain][nstate] = conf
        else:
            nested[strain] = {nstate:conf}
        
    y_positions = [i for i in range(len(strains))]
    x_positions = [i for i in range(len(nstates))]
    print(y_positions)
    plotwidth = 10
    plotheight = 10
    f, ax = plt.subplots(1,1,figsize=(plotwidth, plotheight))
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.xaxis.tick_top()

    height = 0.5
    width = 0.8
    padding  = 0.1
    odd = 1
    for ypos, strain in enumerate(strains):
        print(ypos)
        
        if odd == 1:
            r = Rectangle((- width/2 - 0.1,ypos - height/2 - 0.3),
                          height=height + 0.4,
                          width = len(nstates),
                          facecolor=(0.9,0.9,0.9)#'#d5d8dc'#'#eaecee')
                          )
            ax.add_artist(r)
            odd = 0
        else:
            odd = 1
            
        for xpos, nstate in enumerate(nstates):
            readout_counter = 0
            for conf in nested[strain][nstate]:

                off = conf[0]
                on = conf[1]
                red = lightred#'#faacbd'#'#ff9696'
                green = lightgreen#'#acfacd'#'#9fff96'
                
                if min(on, off) > 0.1:
                    red = darkred#'r'
                    green = darkgreen#'g'
                    
                ystart = ypos - height/2.0
                yend = ypos - height/2.0 + height*off
                X = xpos-width/2.0 + readout_counter*width/6.0
                annotX = X
                annotY = ystart - 0.2
                x = [X for _ in range(2)]
                LW = 4
                BORDER = 0.8
                ax.plot(x, [ystart, yend],
                        color=red,
                        lw=LW,
                        path_effects=[pe.Stroke(linewidth=LW+BORDER, foreground='k'), pe.Normal()])
                ystart = ypos - height/2.0 + height*off
                yend = ypos + height/2.0
                ax.plot(x, [ystart, yend],
                        color=green,
                        lw=LW,
                        path_effects=[pe.Stroke(linewidth=LW+BORDER, foreground='k'), pe.Normal()])
                
                # annotation for experiment agreement
                strain_nstate_readout = strain + '_' + nstate + '_' + READOUTS[readout_counter]
                if strain_nstate_readout in groundtruthdict.keys():
                    gtstate = groundtruthdict[strain_nstate_readout]
                    color = 'r'
                    if gtstate == 'ON':
                        color = 'g'
                    predstate = 'ON'
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
            strainnames.append(mapper[strain])
    ax.set_yticklabels(strainnames)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(nstates)    
    plt.xlabel('Gis1, Mig1, Dot6, Gcn4, Rtg13, Gln3')
    plt.xlim([min(x_positions) - width/2 - padding,max(x_positions) + width/2 + padding])
    plt.ylim([min(y_positions) - height/2 - 4*padding,max(y_positions) + height/2 + padding])    
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)        
    plt.tight_layout()
    plt.savefig('state-space-comparison-robust-representative.pdf', dpi=500)
    
mapper = {'Sch9':'$\Delta$sch9',
          'TORC1':'$\Delta$tor1',
          'Snf1':'$\Delta$snf1',
          'EGO':'$\Delta$gtr1/2',
          'PDE':'$\Delta$pde1/2/3',
          'EGOGAP':'$\Delta$egogap',
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
}        
collect = []
wt = []
notwt = []

# Put wt on top
colorder = reversed(list(df.columns))

for col in colorder:
    if col != 'cost':
        collect.append([col, compare_states(list(df[col]))])

size = df.shape[0]        
visualize(collect, size, groundtruthdict)        
