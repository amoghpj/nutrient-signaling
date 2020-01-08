import matplotlib.pyplot as plt
import yaml
import nutrient_signaling.utils as utils
import numpy as np
import copy
import sys

class TimeCourse:
    """
    
    """
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
            
    def read_data(self, path_to_yaml='../data/time-course-data.yaml'):
        with open(path_to_yaml,'r') as infile:
            self.data = yaml.safe_load(infile)

    def comparison(self):
        """
        Call simulate on each item in yaml file
        """
        print("Starting Comparison")
        store_attributes = ['value','time','readout']
        for simid, experiment in enumerate(self.data):
            self.debug(simid)
            self.predictions.append(self.simulate(experiment))
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
        model = copy.copy(self.simobj)
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
