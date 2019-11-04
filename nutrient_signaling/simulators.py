import os
import sys
import pandas as pd
import PyDSTool as dst
import nutrient_signaling.modelreader as md
from nutrient_signaling.cpputils.pydstool2cpp import PyDSTool2CPP

class SimulatorPython:
    def __init__(self, modeldict):
        self.pars = {}
        self.ics = {}
        self.createModelObject(modeldict)
        
    def simulateModel(self):
        """
        Takes as input PyDSTool Model object
        and returns PyDSTool Pointset with
        default dt=0.01.
    
        :param model: PyDSTool Model object
        :return pts: Solution of Model
        :type pts: dict
        """
        pts = self.model.compute('test').sample()
        return(pts)
    
    def simulate_and_get_points(self):
        return(self.simulateModel())
    
    def get_ss(self):
        """
        Returns steady states of all variables in the model
    
        :return SSPoints: Dictionary containing steady state values
        """
        Points = self.simulateModel()
        SSPoints={}
        for k in Points.keys():
            SSPoints[k]=Points[k][-1]
        return(SSPoints)

    def set_attr(self, pars={}, ics={}, tdata=[0, 90]):
        self.pars.update(pars)
        self.ics.update(ics)
        self.model.set(pars=self.pars, ics=ics, tdata=tdata)

    def get_attr(self):
        print("pars ", self.pars )
        print("ics ", self.ics)
        print("tend ", self.tend )
        
        
    def createModelObject(self, modeldict):
        """
        Takes dictionary object as input, and
        returns a PyDSTool Model object
    
        :param model: Model stored as dictionary. Should contain the keys `variables`, `parameters`, and `initiaconditions`
        :type model: dict
        :return ModelDS: a PyDSTool object that can be simulated to obtain the simulated trajectory
        :type ModelDS: PyDSTool Object
        """
        self.pars = modeldict['parameters']

        self.ics = modeldict['initialconditions']                
        ModelArgs = dst.args(
            name='test',
            varspecs=modeldict['variables'],
            pars=self.pars,
            ics=self.ics,
            tdata=[0, 60],
            ## Define inline functions that are specific to the NutSig modelx
            ## tRNA() computes min(tRNA_total, Amino acid)
            ## pRib() computes min(Rib, eIF)
            ## shs() is the soft-heaviside function
            ## shsm() is the soft heaviside function with a maximum specified.
            fnspecs={'tRNA': (['tRNA_tot', 'AmAc'], 'min(tRNA_tot, AmAc)'),
                     'pRib': (['rib_comp', 'init_factor'],
                              'min(rib_comp,init_factor)'),
                     'shs': (['sig', 'summation'],
                             '1/(1+e^(-sig*summation))'),
                     'shsm': (['sig', 'm', 'summation'],
                              'm/(1+e^(-sig*summation))')}
        )
        self.model = dst.Vode_ODEsystem(ModelArgs)


class SimulatorCPP:
    def __init__(self, modeldict,
                 execpath='./src/',
                 executable='main.o',
                 simfilename='values.dat'):
        self.execpath = execpath
        self.pars = modeldict['parameters']
        self.ics = modeldict['initialconditions']
        self.tdata = [0, 90]
        self.tend = self.tdata[1]
        self.step = 0.001
        self.plot = False
        self.simfilename = simfilename
        self.executable = executable
        
    def set_attr(self, pars={}, ics={},tdata=[0,90]):
        self.pars.update(pars)
        self.ics.update(ics)
        self.tend = tdata[1]
        
    def get_attr(self):
        print("pars ", self.pars )
        print("ics ", self.ics)
        print("tend ", self.tend )
        
    def construct_call(self):
        """
        Returns:
        --------
        cmd : str
            Shell command to call cpp executable
        """
        cmd = self.execpath + self.executable +" --ode --tend " + str(self.tend) + " "\
              + "--step " + str(self.step) + " "
        
        if self.plot:
            cmd += "--plot "
    
        if self.simfilename is not None:
            cmd += "--out " + self.simfilename + " "
            
        if self.pars is not None:
            cmd += "--pars "
            for k,v in self.pars.items():
                cmd += k + " " + str(v)  + " "
                
        if self.ics is not None:
            cmd += "--ics "
            for k,v in self.ics.items():
                cmd += k + " " + str(v)  + " "
        return cmd
    
    def read_sim(self):
        simpoints = pd.read_csv(self.simfilename, sep="\t", index_col = False)
        return simpoints
    
    def get_ss_as_dict(self, simpoints):
        SS = simpoints.tail(1)
        SSdict = SS.to_dict(orient='record')[0]
        SSdict = {k:float(v) for k,v in SSdict.items()}
        return SSdict
    
        
    def simulate(self, cmd):
        so = os.popen(cmd).read()
    
    #########
    ## Helpers
    def simulate_and_get_ss(self):
        cmd = self.construct_call()
        self.simulate(cmd)
        D = self.read_sim()
        SS = self.get_ss_as_dict(D)
        return(SS)
    # wrapper
    def get_ss(self):
        return(self.simulate_and_get_ss())
    
    def simulate_and_get_points(self):
        cmd = self.construct_call()
        self.simulate(cmd)
        D = self.read_sim()
        return(D.to_dict(orient='list'))
    
    def simulateModel(self):
        return(self.simulate_and_get_points())    


def get_simulator(modelpath='./', simulator='py',**kwargs):
    model = md.readinput(modelpath) ## change
    
    if simulator == 'py':
        simobj = SimulatorPython(model)
        return(simobj)
    
    elif simulator == 'cpp':
        validpath = False
        execpath = kwargs.get('execpath','./src/')
        executable = kwargs.get('executable','main.o')
        simfilename = kwargs.get('simfilename','values.dat')
        
        if not os.path.isdir(execpath):
            print('Path to executable does not exist')
            validpath = False
            
        if not os.path.isfile(execpath +  executable):
            print('Executable not found')
            validpath = False
            
        if not validpath:
            if not os.path.exists(execpath + executable):
                print(execpath + executable +' does not exist. Creating model file...')
                cwd = os.path.dirname(os.path.realpath(__file__))
                outpath = execpath
                if not os.path.exists(outpath):
                    os.mkdir(outpath)
                p2c = PyDSTool2CPP(modelpath)
                p2c.setwritepath(cwd + '/cpputils/')
                print(p2c.getwritepath())
                p2c.writecpp()
                headerpath = cwd + '/cpputils/'
                cmd = "g++ -std=c++11 {}main.cpp {}model.cpp -o {}/{}".format(headerpath,
                                                                              headerpath,
                                                                              execpath,
                                                                              executable)
                print(cmd)
                so = os.popen(cmd).read()
                so = os.popen('cp ' + headerpath + executable+ ' ' + outpath)
        simobj = SimulatorCPP(model,
                              execpath=execpath,
                              executable=executable,
                              simfilename=simfilename)
        return(simobj)
            
            
