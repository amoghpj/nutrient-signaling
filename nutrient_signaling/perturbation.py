import yaml
import copy
import sys
import PyDSTool as dst
import matplotlib.pyplot as plt
import os
# Local imports
from nutrient_signaling.utils import minmaxnorm, minmaxnorm_special

class Perturb:
    """All perturbations are some sort of a shift experiment.  The
    simplest case is the wt compared to a single perturbation, be it
    a genetic perturbation or a chemical treatment. Thus,  the control 
    is ALWAYS the wt undergoing the same nutrient shift.
    """
    def __init__(self, simulator, data_path='../data/yaml/perturbation-data.yaml'):        
        #self.order = [12,5,None,13,9,10,4,None,6,8,11,7,3]
        self.order = [12]                       
        self.predictions = []
        self.debugflag = False
        self.simobj = simulator
        self.data_path = data_path
        self.data = []
        self.read_data()
        
    def toggledebug(self):
        self.debugflag = not self.debugflag
        
    def debug(self, prnt):
        if self.debugflag:
            print(prnt)
        
    def read_data(self):
        """
        From path, read perturbation data
        """
        # cwd = os.path.dirname(os.path.realpath(__file__))
        # print('I am in ' + cwd)
        # TODO what is the point of the above line?
        with open(self.data_path,'r') as infile:
            self.data = yaml.safe_load(infile)
            
    def comparison(self):
        """
        Call simulate on each item in yaml file
        """
        print("Starting Comparison")
        store_attributes = ['simid','description', 'units','readout','source','citation','value','whichcondition','type','vmax']
        for simid, experiment in enumerate(self.data):

            # TODO revert
            #if experiment['type'] != 'treat' and experiment['shouldsimulate'] == 'y':
            if experiment['shouldsimulate'] == 'y':                
                predictions = {}
                for att in store_attributes:
                    predictions[att] = experiment[att]
                self.debug(experiment['readout'])                
                self.debug(simid)                
                for experimenttype in ['wt', 'perturb']:
                    self.debug(experimenttype)
                    predictions[experimenttype] = self.simulate(experiment,
                                                                experimenttype=experimenttype)
                    print(predictions)
                self.predictions.append(predictions)
            break

    
    def simulate(self, experiment, experimenttype='wt'):
        """
        For each element in `data`, simulate experimenttype
        and store the steady state values
        """
        readout = experiment['readout']
        pars_perturb = experiment['spec']['perturbation']['pars']
        self.debug(pars_perturb)
        ics_perturb = experiment['spec']['perturbation']['ics']
        self.debug(ics_perturb)
        
        result = {'pre':0, 'post':0}

        model = copy.copy(self.simobj)
        self.debug('\t\tinitialized')
        
        # Initialize all variables with their ss values
        # Preshift
        preics = model.get_ss()
        print(preics[readout])
        prepars = experiment['spec']['preshift']['pars']
        model.set_attr(pars=prepars)

        #self.debug('\t\tset ics to ss')
        
        if experimenttype == 'perturb':
            # Parameters and initial conditions representing the experimental
            preics.update(ics_perturb)
            prepars.update(pars_perturb)

        model.set_attr(ics=preics,
                       pars=prepars,                       
                       tdata=[0, experiment['time'][experimenttype]])
        
        result['pre'] = model.get_ss()[experiment['readout']]
        print(result['pre'])
        self.debug('pre done')
        # Postshift
        postics = model.get_ss()
        postpars = experiment['spec']['postshift']['pars']
        
        if experimenttype == 'perturb':
            postics.update(ics_perturb)
            postpars.update(pars_perturb)
            
        model.set_attr(ics=postics,
                       pars=postpars,
                       tdata=[0, experiment['time'][experimenttype]])
        
        result['post'] = model.get_ss()[experiment['readout']]
        print(result['post'])
        self.debug('post done')        
        del(model)
        self.debug(result)
        return(result)

    def simulate_old(self,experiment, experimenttype='wt'):
        """
        For each element in `data`, simulate experimenttype
        and store the steady state values
        """
        so = SimulatorPython()
        result = {'pre':0, 'post':0}        
        modelds = so.createModelObject(copy.deepcopy(self.model))

        self.debug('\t\tinitialized')
        # Initialize all variables with their ss values        
        modelds.set(ics=so.get_ss(modelds))
        self.debug('\t\tset ics to ss')
        if experimenttype == 'perturb':
            # Parameters and initial conditions representing the experimental
            pars_perturb = experiment['spec']['perturbation']['pars']
            self.debug(pars_perturb)
            ics_perturb = experiment['spec']['perturbation']['ics']
            self.debug(ics_perturb)            
            modelds.set(pars=pars_perturb,
                        ics=ics_perturb)

            self.debug('perturb done')

        # Preshift 
        modelds.set(ics=so.get_ss(modelds),
                    tdata=[0, experiment['time'][experimenttype]],
                    pars=experiment['spec']['preshift']['pars'])
        
        result['pre'] = so.get_ss(modelds)[experiment['readout']]
        self.debug('pre done')
        # Postshift
        modelds.set(ics=so.get_ss(modelds),
                    tdata=[0, experiment['time'][experimenttype]],
                    pars=experiment['spec']['postshift']['pars'])
        result['post'] = so.get_ss(modelds)[experiment['readout']]
        self.debug('post done')        
        del(modelds)
        self.debug(result)
        return(result)

    ## Plotting utilies
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

    def visualize_perturb(self):
        """
        Specifies visualization of perturbation data
        """

        scale = 0.8
        f,ax = plt.subplots(1,2,
                          figsize=(5, 0.5*len(self.predictions)))

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
                prediction = self.predictions[oi]
                # for p in self.predictions:
                #     if p['simid'] == oi:
                #         prediction = p
                #         break
                cond = prediction['whichcondition']
                leftylabels.append('('+cond+') ' + predictions['readout'])
                # If the maxium value from an experiment is specified,
                # use it for minmax normaliation, otherwise compute min and
                # max from the list of values
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
        ax[0].set_title("Wild-type")        
        ax[1].set_title("Mutant")
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
