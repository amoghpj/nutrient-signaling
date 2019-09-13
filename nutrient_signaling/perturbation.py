import yaml
import copy
import sys
import PyDSTool as dst
import matplotlib.pyplot as plt
import os
from nutrient_signaling.simulators import SimulatorPython

class Perturb:
    """All perturbations are some sort of a shift experiment.  The
    simplest case is the wt compared to a single perturbation, be it
    a genetic perturbation or a chemical treatment. Thus,  the control 
    is ALWAYS the wt undergoing the same nutrient shift.
    """
    #def __init__(self, model, simulator):
        #self.model = model
    def __init__(self, simulator):        
        self.data = []
        self.predictions = []
        self.debugflag = False
        self.simobj = simulator
        
    def toggledebug(self):
        self.debugflag = not self.debugflag
        
    def debug(self, prnt):
        if self.debugflag:
            print(prnt)
        
    def read_data(self, path_to_yaml):
        """
        From path, read perturbation data
        """
        cwd = os.path.dirname(os.path.realpath(__file__))
        print('I am in ' + cwd)
        with open(path_to_yaml,'r') as infile:
            self.data = yaml.safe_load(infile)
            
    def comparison(self):
        """
        Call simulate on each item in yaml file
        """
        store_attributes = ['simid','description', 'units','readout','source','citation','value','whichcondition','type','vmax']
        for simid, experiment in enumerate(self.data):
            predictions = {}
            self.debug(simid)
            if experiment['type'] != 'treat' and experiment['simulate'] == 'y':
                for att in store_attributes:
                    predictions[att] = experiment[att]
                self.debug(experiment['readout'])
                for experimenttype in ['wt', 'perturb']:
                    self.debug(experimenttype)
                    predictions[experimenttype] = self.simulate(experiment,
                                                                experimenttype=experimenttype)
                        
                self.predictions.append(predictions)

    
    def simulate(self,experiment, experimenttype='wt'):
        """
        For each element in `data`, simulate experimenttype
        and store the steady state values
        """
        result = {'pre':0, 'post':0}

        model = copy.deepcopy(self.simobj)
        # print('init')
        # model.get_attr()
        self.debug('\t\tinitialized')
        # Initialize all variables with their ss values        
        model.set_attr(ics=model.get_ss())
        # print('ics')
        # model.get_attr()        
        self.debug('\t\tset ics to ss')
        if experimenttype == 'perturb':
            # Parameters and initial conditions representing the experimental
            pars_perturb = experiment['spec']['perturbation']['pars']
            self.debug(pars_perturb)
            ics_perturb = experiment['spec']['perturbation']['ics']
            self.debug(ics_perturb)            
            model.set_attr(pars=pars_perturb,
                           ics=ics_perturb)
            # print('ics')
            # model.get_attr()        
            self.debug('perturb done')

        # Preshift 
        model.set_attr(ics=model.get_ss(),
                       tdata=[0, experiment['time'][experimenttype]],
                       pars=experiment['spec']['preshift']['pars'])
        
        # if experiment['description'] == 'rapamycin' and experimenttype == 'perturb':
        #     print('preshift')
        #     model.get_attr()
        #     sys.exit()
        result['pre'] = model.get_ss()[experiment['readout']]
        self.debug('pre done')
        # Postshift
        model.set_attr(ics=model.get_ss(),
                       tdata=[0, experiment['time'][experimenttype]],
                       pars=experiment['spec']['postshift']['pars'])
        result['post'] = model.get_ss()[experiment['readout']]
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
    
