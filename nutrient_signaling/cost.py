from nutrient_signaling.perturbation import Perturb
from nutrient_signaling.timecourse import TimeCourse
import nutrient_signaling.simulators as sim
import nutrient_signaling.utils as utils
import numpy as np
import sys
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

class ComputeCost:
    """
    Compute cost of fit to time course and perturbation data.
    """
    def __init__(self, modelpath, simulator, psetpath, perturbpath = None,numpsets=10, executable="main.o"):
        self.pert = Perturb()
        self.perturbpath = None
        self.pert.readData(data_path=perturbpath)
        self.modelpath = modelpath
        self.simulator = simulator
        self.executable = executable

        self.tc = TimeCourse()
        self.tc.readData()
        # self.perturb_data = Perturb
        self.psetpath = psetpath
        self.numpsets = numpsets
        self.timeCourseSpecList = [
            'Snf1:glucose_repression',
            'wt:cAMP',
            'sch9Delta:cAMP',
            'Sch9P:100_glutamine',
            'Sch9P:1_glutamine',
            'Sch9:gtr1Delta',
            'Sch9:glucose_starve',
            'Sch9:glucose_relief',
            'Mig1:glucose_relief',
            'Rib:rap',
            'Gln3:rap'
        ]
        
    def read_pset(self):
        self.psetdf = pd.read_csv(self.psetpath, header=0,sep='\t').iloc[0:self.numpsets] #.sample(self.numpsets, axis='index')

    def compute_costs(self, outname="out"):
        self.read_pset()
        store_costs = {'cost':[], 'c_t': [], 'c_p':[],
                       'original_cost':[], 'original_c_t':[],
        'original_c_p':[]}
        for i, pset in tqdm(self.psetdf.iterrows()):
            pdict = pset.to_dict()
            c, ct, cp = self.compute(pdict)
            store_costs['cost'].append(c)
            store_costs['c_t'].append(ct)
            store_costs['c_p'].append(cp)
            store_costs['original_cost'].append(pset['cost'])
            store_costs['original_c_t'].append(pset['c_t'])
            store_costs['original_c_p'].append(pset['c_p'])
        storedf = pd.DataFrame(store_costs)
        storedf.to_csv(outname + '.csv')
        self.visualize_costs(storedf)

    def visualize_costs(self, storedf):
        fig = plt.figure(figsize=(8,12))
        ax = fig.add_subplot(3,2,1)
        ax.hist(storedf.c_t, bins=10,color='g', alpha=0.2, label='c_t')
        ax.legend()        
        ax = fig.add_subplot(3,2,2)        
        ax.hist(storedf.original_c_t, bins=10,color='k', alpha=0.2, label='o_c_t')
        ax.legend()
        # pert cost
        ax = fig.add_subplot(3,2,3)
        ax.hist(storedf.c_p, color='g', bins=10,alpha=0.2, label='c_p')
        ax.legend()
        ax = fig.add_subplot(3,2,4)        
        ax.hist(storedf.original_c_p, bins=10,color='k', alpha=0.2, label='o_c_p')
        ax.legend()
        # total cost
        ax = fig.add_subplot(3,2,5)
        ax.hist(storedf.cost, color='g', bins=10,alpha=0.2, label='c')
        ax.legend()
        ax = fig.add_subplot(3,2,6)
        ax.hist(storedf.original_cost, bins=10,color='k', alpha=0.2, label='o_c')
        ax.legend()
        plt.tight_layout()
        plt.savefig('cost-dists.png')
        
    def compute(self, paramset={}):
        """
        Computes the total time-course cost and perturbation
        cost for the given parameter set. If the argument is 
        empty, uses default parameter set.
        """
        cost_ts = self.compute_timecourse_cost(paramset)
        cost_pert = self.compute_perturb_cost(paramset)
        cost_full = cost_ts + cost_pert
        return(cost_full, cost_ts, cost_pert)
    
    def _get_list_item(self, st, li, key):
        """ returns index"""
        for i, item in enumerate(li):
            if item[key] == st:
                return i
            
    def _findindex(self, t, expt,tol=1e-1):
        """FIXME! briefly describe function
        :param t: 
        :param expt: 
        :param tol: 
        :returns: 
        :rtype: 
        """
        idx = 0
        for i, _t in enumerate(t):
            if abs(_t - expt) < tol:
                return i
            elif _t > expt:
                return i-1
            else:
                continue
        
    def compute_perturb_cost(self, paramset):
        self.pert.setSimulator(modelpath=self.modelpath,
                               simulatortype=self.simulator,
                               executable=self.executable)
        cost = 0.
        perturb_count = 0.
        for simid in [sid for sid in self.pert.order if sid is not None] :
            exp_data = self.pert.data[self._get_list_item(simid, self.pert.data, "simid")]
            # This is an old check, no simid in order should correspond to 'treat'
            if exp_data['type'] != 'treat' and exp_data['shouldsimulate'] == 'y':
                e_wt_pre, e_wt_post, e_mut_pre, e_mut_post = utils.minmaxnorm([
                    exp_data["value"]["preshift"]["wt"],
                    exp_data["value"]["postshift"]["wt"],                    
                    exp_data["value"]["preshift"]["perturb"],
                    exp_data["value"]["postshift"]["perturb"]])
                
                predictions = {}
                for experimenttype in ['wt', 'perturb']:
                    predictions[experimenttype] = self.pert.simulate(exp_data,
                                                                     experimenttype=experimenttype,
                                                                     paramset=paramset)
                if exp_data['readout'] in ['cAMP', 'Rib']:
                    p_wt_pre, p_wt_post, p_mut_pre, p_mut_post  = utils.minmaxnorm([
                        predictions["wt"]["pre"],
                        predictions["wt"]["post"],
                        predictions["perturb"]["pre"],
                        predictions["perturb"]["post"]])
                else:
                    p_wt_pre, p_wt_post, p_mut_pre, p_mut_post = predictions["wt"]["pre"],\
                        predictions["wt"]["post"],\
                        predictions["perturb"]["pre"],\
                        predictions["perturb"]["post"]
                    # if exp_data['vmax'] > 0:
                    #     p_wt_pre, p_wt_post, p_mut_pre, p_mut_post  = utils.minmaxnorm_special([
                    #         predictions["wt"]["pre"],
                    #         predictions["wt"]["post"],
                    #         predictions["perturb"]["pre"],
                    #         predictions["perturb"]["post"]], 0.0, exp_data['vmax'])
                    # else:            
                
                # compute cost
                # print(exp_data['description'], exp_data['readout'])
                # print(e_wt_pre, p_wt_pre)
                # print(e_wt_post, p_wt_post)
                # print(e_mut_pre, p_mut_pre)
                # print(e_mut_post, p_mut_post)
                # print('---')
                cost_wt = 0.5*((e_wt_pre - p_wt_pre)**2. + (e_wt_post - p_wt_post)**2.)
                cost_mut = 0.5*((e_mut_pre - p_mut_pre)**2. + (e_mut_post - p_mut_post)**2.)
                cost += 0.5*(cost_wt + cost_mut)
                perturb_count += 1.
        return(cost/float(perturb_count))
        
    def compute_timecourse_cost(self, paramset):
        cost = 0.
        timepoint_count = 0.
        
        for spec in self.timeCourseSpecList:
            exp_data_idx = self._get_list_item(spec,
                                           self.tc.data,
                                           'description')
            exp_data = self.tc.data[exp_data_idx]
            simspec = exp_data["spec"]
            expT = exp_data['time']
            expV = exp_data['value']
            if exp_data['tunits'] == 's':
                expT = [t / 60.0 for t in expT]
                
            if exp_data["ts_type"] == "ts":
                if np.isnan(exp_data['normalize']['lower']) and np.isnan(exp_data['normalize']['upper']):
                    expV = utils.minmaxnorm(expV)
                else:
                    expV = utils.minmaxnorm_special(expV,
                                                    exp_data['normalize']['lower'],
                                                    exp_data['normalize']['upper'])

            varname = exp_data['readout']
            
            parameters = paramset.copy()
            parameters.update(simspec['preshift']['pars'])
            ## Simulation details
            model = sim.get_simulator(self.modelpath, self.simulator, executable=self.executable)
            preshift_ics = model.simulate_and_get_ss()
            preshift_ics.update(simspec["preshift"]["ics"])
            model.set_attr(tdata=[0,90], pars=parameters, ics=preshift_ics)
            
            postshift_ics = model.simulate_and_get_ss()
            
            # post shift
            postshift_ics.update(simspec['postshift']['ics'])
            parameters.update(simspec['postshift']['pars'])
            model.set_attr(tdata=[0, expT[-1] + 10.],
                           pars=parameters,
                           ics=postshift_ics)
            
            pts = model.simulate_and_get_points()

            if 'func' in exp_data['ts_type']:
                predicted = self.tc.customFuncDef[exp_data['ts_type'].split(':')[1]](pts[varname])
            else:
                predicted = list(pts[varname])
            pred_time = pts['t']
            
            # for every time point in the ground truth, find
            # the corresponding index in the predicted time list,
            # and then compute the deviation for the predicted value
            #
            cost += sum([np.power(predicted[self._findindex(pred_time, et)] - ev, 2) \
                        for et, ev in zip(expT, expV)])
            timepoint_count += len(expT)
        return(cost/timepoint_count)
