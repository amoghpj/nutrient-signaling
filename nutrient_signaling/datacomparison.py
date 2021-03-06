import os
import yaml
import sys
import matplotlib.pyplot as plt
from nutrient_signaling.perturbation import Perturb
from nutrient_signaling.timecourse import TimeCourse
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
        self.checkpaths()
        if self.comparisontype not in ['time', 'perturb', 'qualtitative']:
            print("comparisontype must be one of ['time','perturb','qualtitative']")
            raise Warning
            #sys.exit()

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
        """Compare model predictions to either perturbation or time course data"""
        if self.comparisontype == 'perturb':
            compobj = Perturb()
            if self.debug:
                compobj.toggledebug()
            compobj.read_data(self.datapath)
            compobj.setSimulator(self.modelpath, self.simulator)
            compobj.comparison()
            compobj.visualize_perturb()

        if self.comparisontype == 'time':
            compobj = TimeCourse()
            compobj.read_data(self.datapath)
            compobj.setSimulator(self.modelpath, self.simulator)
            if self.debug:
                compobj.toggledebug()
            compobj.comparison()
            # compobj.plotall()

