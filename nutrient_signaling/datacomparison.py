import os
import yaml
import sys
import matplotlib.pyplot as plt
from nutrient_signaling.perturbation import Perturb
from nutrient_signaling.simulators import get_simulator
from nutrient_signaling.utils import processliteral,\
    minmaxnorm

class Comparison:
    def __init__(self, datapath, modelpath, comparisontype, simulator='py',debug=False):
        self.datapath = datapath
        self.modelpath = modelpath
        self.comparisontype = comparisontype
        if self.comparisontype not in ['time','perturb','qualtitative']:
            print("comparisontype must be one of ['time','perturb','qualtitative']")
        self.simulator = simulator
        self.checkpaths()
        self.docomparison()
        self.debug = debug
        
    def checkpaths(self):
        validpath = True
        if not os.path.exists(self.datapath) and not os.path.isfile(self.datapath):
            print('Invalid datapath')
            validpath = False
        if not os.path.exists(self.modelpath) and not os.path.isdir(self.datapath):
            print('Invalid modelpath')
            validpath = False
        if validpath:
            print('Paths Validated.')
        else:
            sys.exit()

    def docomparison(self):
        if self.comparisontype == 'perturb':
            simobj = get_simulator(self.modelpath, self.simulator)
            compobj = Perturb(simobj)
            compobj.toggledebug()
            print('loading data')
            
            compobj.read_data(self.datapath)
            print('simulating model')
            compobj.comparison()
            print(compobj.predictions)
            self.plot_perturb(compobj.predictions)
            
    def plot_perturb(self, predictions):
        """
        Specifies visualization of perturbation data
        """
        print(len(predictions))
        f,ax = plt.subplots(1,2,
                          figsize=(6, len(predictions)
                          ))
        ypos = 0
        ax[0].set_xlim(-0.1,1.1)
        ax[1].set_xlim(-0.1,1.1)
        leftylabels = []
        for ypos, prediction in enumerate(predictions):
            print(prediction['value'])
            leftylabels.append(prediction['readout'])
            expvalues = {'preshift':{'wt':0,'perturb':0},
                         'postshift':{'wt':0,'perturb':0}}
            expvalues['preshift']['wt'],\
                expvalues['postshift']['wt'],\
                expvalues['preshift']['perturb'],\
                expvalues['postshift']['perturb'] = minmaxnorm([prediction['value']['preshift']['wt'],
                                                            prediction['value']['postshift']['wt'],
                                                            prediction['value']['preshift']['perturb'],
                                                            prediction['value']['postshift']['perturb']])


            linestyle = {'exp':{'wt':'k-','perturb':'r-','disp':0.15},
                         'pred':{'wt':'k--','perturb':'r--','disp':-0.15}}
            
            for axid, experimenttype in enumerate(['wt', 'perturb']):
                ax[axid].plot([expvalues['preshift'][experimenttype],
                               expvalues['postshift'][experimenttype]],
                              [ypos+linestyle['exp']['disp'],
                               ypos+linestyle['exp']['disp']],
                              linestyle['exp'][experimenttype])
                ax[axid].plot([prediction[experimenttype]['pre'],
                               prediction[experimenttype]['post']],
                              [ypos+linestyle['pred']['disp'],
                               ypos+linestyle['pred']['disp']],
                              linestyle['pred'][experimenttype])
        ax[0].set_yticks([i for i in range(len(leftylabels))], )
        ax[0].set_yticklabels(leftylabels)
        ax[1].set_yticklabels([])
        plt.tight_layout()
        plt.savefig('perturb-vis.png')            
