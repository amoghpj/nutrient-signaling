import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# Local Imports
from nutrient_signaling.simulators import get_simulator
font = {'size'   : 16}
matplotlib.rc('font', **font)

f, ax = plt.subplots(1,2, figsize=(10,5))
LW = 4

## Prouteau, 2017
ax[0].plot([0.0,2.5,5.0,15.0,30.0],
           [1.0,0.34,0.13,0.15,0.20],
           'ko',markersize=8)
ax[1].plot([0.0,2.0,5.0,15.0,30.0],
           [0.20,0.63,0.92,0.91,1.0],
           'ko',
           markersize=8)

# Plot time course data
experiments = [{'pre':{'parameters':{'Carbon':1.0,
                                     'ATP':1.0,
                                     'Glutamine_ext':1.0},
                       'inconds':{}},
                'post':{'parameters':{'Carbon':0.0,
                                      'ATP':0.0,
                                      'Glutamine_ext':1.0},
                        'inconds':{}}},
               {'pre':{'parameters':{'Carbon':0.0,
                                     'ATP':0.0,
                                     'Glutamine_ext':1.0},
                       'inconds':{}},
                'post':{'parameters':{'Carbon':1.0,
                                      'ATP':1.0,
                                      'Glutamine_ext':1.0},
                        'inconds':{}}}
]

plotsettings = [{'style':'k--','label':'Model 1'}, # 'Snf1 -| TORC1'
                {'style':'#bbbbbb','label':'Model 2'}] #'no interaction'

MODELPATH = 'data/2018-9-26-12-3-no-sigma/'

for exp, axis in zip(experiments, ax):
    # For each experiment, create a wt  and mutant object,
    # where the mutant strain lacks the TORC1-Snf1 interaction.
    # Simulate the shift defined in the experiment, plot
    # the time courses
    wt = get_simulator(modelpath=MODELPATH,
                      SIMULATOR = 'cpp',
                      SIMARGS = {'executable':'robust.o'})

    mut = get_simulator(modelpath=MODELPATH,
                        SIMULATOR = 'cpp',
                        SIMARGS = {'executable':'robust.o'})
    mut.set_attr(pars={'w_torc_snf':0.0})
    strains = [wt, mut]
    for model, setting in zip(strains, plotsettings):
        # For both models, set the preshift parameters and ics
        model.set_attr(pars=exp['pre']['parameters'],
                       ics=exp['pre']['inconds'])
        ss = model.get_ss()
        model.set_attr(pars=exp['post']['parameters'],
                       ics=model.get_ss(),
                       tdata=[0,30])
        traj = model.simulate_and_get_points()
        axis.plot(traj['t'], traj['Sch9'],
                  setting['style'],
                  label=setting['label'],
                  lw=LW)
        


ax[0].set_ylabel('Sch9-P')
ax[0].legend()
ax[0].set_xlabel('time (min)')
ax[0].set_title('Acute Glucose Starvation')
ax[1].set_xlabel('time (min)')
ax[1].set_yticks([])
ax[1].set_title('Relief from Glucose Starvation')

#plt.suptitle('Snf1 mediated TORC1 inhibition')
outpath = 'img/Snf1-TORC1.png'
plt.tight_layout()
plt.savefig(outpath, dpi=300)

