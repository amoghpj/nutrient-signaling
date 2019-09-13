import os
import yaml
import sys
import matplotlib.pyplot as plt
from nutrient_signaling.perturbation import Perturb
from nutrient_signaling.simulators import get_simulator
from nutrient_signaling.utils import processliteral,\
    minmaxnorm, minmaxnorm_special

class Comparison:
    def __init__(self, datapath, modelpath, comparisontype, simulator='py'):
        self.datapath = datapath
        self.modelpath = modelpath
        self.comparisontype = comparisontype
        self.debug = False                
        self.simulator = simulator 
        self.order = [12,5,None,13,9,10,4,None,6,8,11,7,3]       
        self.checkpaths()
        if self.comparisontype not in ['time','perturb','qualtitative']:
            print("comparisontype must be one of ['time','perturb','qualtitative']")
            sys.exit()

    def toggledebug(self):
        self.debug = True

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
            if self.debug:
                compobj.toggledebug()
            print('loading data')
            compobj.read_data(self.datapath)
            print('simulating model')
            compobj.comparison()
            self.visualize_perturb(compobj.predictions)

    def perturb_plot_row(self, values_to_plot, ax, ypos):
        style = {'exp':{'wtc':'k','ls':'-','perturbc':'r','disp':0.15},
                 'pred':{'wtc':'k','ls':'--','perturbc':'r','disp':-0.15}}
            
        for axid, etype in enumerate(['wt', 'perturb']):
            for dtype in ['exp', 'pred']:
                lw = 2.0                    
                ax[axid].plot([values_to_plot[etype][dtype]['preshift'],
                               values_to_plot[etype][dtype]['postshift']],
                              [ypos + style[dtype]['disp'],
                               ypos + style[dtype]['disp']],
                              c=style[dtype][etype + 'c'],
                              ls=style[dtype]['ls'],
                              lw=lw)
                # draw arrowheads
                epre, epost = values_to_plot[etype][dtype]['preshift'],\
                    values_to_plot[etype][dtype]['postshift']
                leftmarker = '.'
                rightmarker = '>'
                msize = 10
                left, right = float(epre), float(epost)
                mcolor = style[dtype][etype + 'c']
                if epost < epre:
                    leftmarker = '<'
                    rightmarker = '.'
                    left, right = float(epost), float(epre)
                if dtype == 'pred':
                    mcolor = 'white'
                for v, shape in [(left, leftmarker), (right, rightmarker)]:
                    ax[axid].plot(v, ypos + style[dtype]['disp'],
                                  marker=shape,
                                  c=style[dtype][etype +'c'],
                                  markeredgecolor=style[dtype][etype +'c'],
                                  markerfacecolor=mcolor,
                                  markersize=msize)

    def visualize_perturb(self, predictions):
        """
        Specifies visualization of perturbation data
        """
        scale = 0.8
        f,ax = plt.subplots(1,2,
                          figsize=(5, 0.5*len(predictions)))

        orderedpredictions = {'C':[],'N':[]}        
        leftylabels = []
        descriptions = []
        ypos = 0        
        for oi in self.order[::-1]:
            if oi is None:
                leftylabels.append("")
                descriptions.append("")
                ypos += 1
            else:
                for p in predictions:
                    if p['simid'] == oi:
                        prediction = p
                        break
                cond = prediction['whichcondition']
                leftylabels.append('('+cond+') ' + prediction['readout'])
                if prediction['vmax'] != 0:
                    w_e_pre, w_e_post, p_e_pre, p_e_post = minmaxnorm_special([prediction['value']['preshift']['wt'],
                                                                               prediction['value']['postshift']['wt'],
                                                                               prediction['value']['preshift']['perturb'],
                                                                               prediction['value']['postshift']['perturb']],0,
                                                                              prediction['vmax'])
                else:
                    w_e_pre, w_e_post, p_e_pre, p_e_post = minmaxnorm([prediction['value']['preshift']['wt'],
                                                                       prediction['value']['postshift']['wt'],
                                                                       prediction['value']['preshift']['perturb'],
                                                                       prediction['value']['postshift']['perturb']])                    
                w_p_pre, w_p_post ,p_p_pre, p_p_post = prediction['wt']['pre'],\
                    prediction['wt']['post'],\
                    prediction['perturb']['pre'],\
                    prediction['perturb']['post']
                
                if prediction['readout'] in ['cAMP','Rib']:
                    w_p_pre, w_p_post ,p_p_pre, p_p_post = minmaxnorm([w_p_pre,
                                                                       w_p_post,
                                                                       p_p_pre,
                                                                       p_p_post])
                
                values_to_plot = {'wt':{'exp':{'preshift':w_e_pre,
                                               'postshift':w_e_post},
                                        'pred':{'preshift':w_p_pre, 'postshift':w_p_post}},
                                  'perturb':{'exp':{'preshift':p_e_pre,
                                                    'postshift':p_e_post},
                                             'pred':{'preshift':p_p_pre, 'postshift':p_p_post}}}
                
    
                descriptions.append(prediction['description'])
                self.perturb_plot_row(values_to_plot, ax, ypos)
                ypos += 1

                
        ax[0].set_xlim(-0.1,1.1)
        ax[1].set_xlim(-0.1,1.1)
        ax[0].set_title("Control")        
        ax[1].set_title("Perturbation")
        ax[0].set_xticks([0., 0.5, 1.0])
        ax[0].set_xticklabels([0., 0.5, 1.0])
        ax[1].set_xticks([0., 0.5, 1.0])
        ax[1].set_xticklabels([0., 0.5, 1.0])        
        ax[0].set_yticks([i for i in range(len(leftylabels))], )        
        ax[0].set_yticklabels(leftylabels)
        ax[1].set_yticks([i for i in range(len(leftylabels))], )        
        ax[1].set_yticklabels(['' for i in range(len(leftylabels))])
        ax_right_twin=ax[1].twinx()
        ax_right_twin.set_ylim(ax[1].get_ylim())
        ax_right_twin.set_yticks(ax[1].get_yticks())
        ax_right_twin.set_yticklabels(descriptions)
        ax[1].set_yticks([])
        ax[1].set_yticklabels([])
        plt.tight_layout()
        plt.savefig('perturb-vis.pdf')            
