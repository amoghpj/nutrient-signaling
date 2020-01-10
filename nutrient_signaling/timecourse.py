import matplotlib.pyplot as plt
import yaml
import nutrient_signaling.utils as utils
from nutrient_signaling.simulators import get_simulator
import numpy as np
import copy
import sys

class TimeCourse:
    """
    
    """
    def __init__(self):
        self.data = []
        self.predictions = []
        self.debugflag = False
        self.data_path = '../data/yaml/time-course-data.yaml'
        self.customFuncDef = {'mig1': (lambda mig_ts: [np.log10((m)/(1e-3+1-m)) for m in mig_ts]), # ensures denom is not 0
                 'rib': (lambda rib_ts: [rib_ts[i]/(1e-2+rib_ts[0]) for i in range(0,len(rib_ts))])} #ensures denom is not 0

        self.customylabel = {'mig1' : 'log(nMig1/cMig1)',
                    'rib':'Rel Rib'}
        
    
    def toggledebug(self):
        self.debugflag = not self.debugflag
        
    def setSimulator(self, modelpath, simulatortype='py'):
        self.modelpath = modelpath
        self.simulatortype = simulatortype

    def makeSimulator(self):
        simobj = get_simulator(self.modelpath, self.simulatortype)
        return simobj
        
    def debug(self, prnt):
        if self.debugflag:
            print(prnt)
            
    def read_data(self, data_path='../data/yaml/perturbation-data.yaml'):
        self.data_path = data_path
        with open(self.data_path,'r') as infile:
            self.data = yaml.safe_load(infile)

    def comparison(self):
        """
        Call simulate on each item in yaml file
        TODO: This is getting messy again. 
        1. Need to handle plotting of Mig1 on log scale
        2. Need to handle non-normalized y axis for mig1 cleanly
        3. cAMP should always have 0-1 range.
        """
        print("Starting Comparison")
        store_attributes = ['value','time','readout']
        for simid, experiment in enumerate(self.data):
            self.debug(simid)
            traj = self.simulate(experiment)
            if experiment['tunits'] == 's':
                traj['t'] = [t*60. for t in traj['t']]
            if experiment['readout'].lower() in self.customFuncDef.keys():
                traj['v'] = self.customFuncDef[experiment['readout'].lower()](traj['v'])
            self.predictions.append(traj)
            f, ax = plt.subplots(1,1)
            self.plotdata(ax, experiment)
            self.plotpred(ax, self.predictions[-1])
            plt.savefig(f'img/{simid}.png')
            plt.close()

    def plotall(self):
        for simid, (experiment, pred) in enumerate(zip(self.data, self.predictions)):
            f, ax = plt.subplots(1,1)
            self.plotdata(ax, experiment)
            self.plotpred(ax, pred)
            plt.savefig(f'img/{simid}.png')
            plt.close()
        
    def plotdata(self, ax, experiment):
        self.debug(experiment['title'])
        lo = experiment['normalize']['lower']
        hi = experiment['normalize']['upper']
        if np.isnan(lo):
            lo = min(experiment['value'])
        if np.isnan(hi):
            hi = max(experiment['value'])

        normalizedData = utils.minmaxnorm_special(experiment['value'], lo, hi)
        ax.plot(experiment['time'], normalizedData, 'k.')
        ax.set_title(experiment['title'])
        
    def plotpred(self, ax, traj):
        print(len(traj['t']))
        ax.plot(traj['t'], traj['v'],'k--')
            
    def simulate(self,experiment):
        """
        Simulate a time course using the preshift and post shift specifications defined in
        the experiment
        """
        print(experiment['readout'])
        model = self.makeSimulator()
        preshift = experiment['spec']['preshift']
        postshift = experiment['spec']['postshift']
        print(preshift)
        print(postshift)
        tmax = 90
        if len(experiment['time']) >0:
            tmax = max(experiment['time'])
        
        if experiment['tunits'] == 's':
            tmax = tmax/60.

        # Set up preshift
        preics = model.get_ss()
        #preics.update(preshift['ics'])
        model.set_attr(ics=preics, pars=preshift['pars'])
        # set up post shift
        postics = model.get_ss()
        #postics.update(postshift['ics'])

        model.set_attr(ics=postics, pars=postshift['pars'], tdata=[0,tmax])
        y = model.simulate_and_get_points()
        traj = {'t':y['t'],
                'v':y[experiment['readout']]}

        return traj
            
    def plot_predicted(self, ax, traj, plotargs=None):
        """
        traj should be dictionary with keys 't' and 'v' with values being lists
        containing the time points and the values respectively
        """
        ax.plot(traj['t'], traj['v'])
