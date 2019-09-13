import os
import pandas as pd
import PyDSTool as dst
import nutrient_signaling.modelreader as md
import sys
class SimulatorPython:
    def __init__(self, modeldict):
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
        self.model.set(pars=pars, ics=ics, tdata=tdata)
        
    def createModelObject(self, modeldict):
        """
        Takes dictionary object as input, and
        returns a PyDSTool Model object
    
        :param model: Model stored as dictionary. Should contain the keys `variables`, `parameters`, and `initiaconditions`
        :type model: dict
        :return ModelDS: a PyDSTool object that can be simulated to obtain the simulated trajectory
        :type ModelDS: PyDSTool Object
        """
        ModelArgs = dst.args(
            name='test',
            varspecs=modeldict['variables'],
            pars=modeldict['parameters'],
            ics=modeldict['initialconditions'],
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
    def __init__(self, execpath='./src/',executable='main.o'):
        self.execpath = execpath
        self.pars = {}
        self.ics = {}
        self.tdata = [0, 90]
        self.tend = self.tdata[1]
        self.step = 0.001
        self.plot = False
        self.simfilename = 'values.dat'
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
    if simulator == 'py':
        model = md.readinput(modelpath) ## change
        simobj = SimulatorPython(model)
        return(simobj)
    elif simulator == 'cpp':
        validpath = False
        execpath = './src/'
        executable = 'main.o'
        if 'execpath' in kwargs.keys():
            execpath = kwargs['execpath']
            if not os.path.isdir(kwargs['execpath']):
                print('Path to executable does not exist')
                validpath = False
        if 'executable' in kwargs.keys():
            executable = kwargs['executable']
            if not os.path.isfile(kwargs['execpath'] +  kwargs['executable']):
                print('Executable not found')
                validpath = False
        if not validpath:
            if not os.path.exists(execpath + executable):
                print(execpath + executable +' does not exist. Creating model file...')
                cwd = os.path.dirname(os.path.realpath(__file__))
                from nutrient_signaling.cpputils.pydstool2cpp import PyDSTool2CPP
                outpath = execpath#cwd + '/../' + execpath
                if not os.path.exists(outpath):
                    os.mkdir(outpath)
                p2c = PyDSTool2CPP(modelpath)
                p2c.setwritepath(cwd + '/cpputils/')
                print(p2c.getwritepath())
                p2c.writecpp()
                os.chdir(cwd + '/cpputils')
                cmd = "g++ -std=c++11 main.cpp model.cpp -o " + executable
                so = os.popen(cmd).read()
                so = os.popen('cp ' +executable+ ' ' + outpath)
                os.chdir(cwd + '/../')
        simobj = SimulatorCPP(execpath, executable )
        return(simobj)
            
            
