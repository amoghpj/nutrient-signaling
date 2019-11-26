import os
import matplotlib.pyplot as plt
HOME = os.path.expanduser('~')
os.chdir(HOME + '/jalihal_projects/Research/data/ModelAnalysis/')
import modelreader as md
import simulator as sim

Model = md.main('./ParameterSets/2018-9-26-12-3-no-sigma/')
f, ax = plt.subplots(1,2,figsize=(8,4))
## Prouteau, 2017
pre = {'parameters':{'Carbon':1.0,'ATP':1.0,'Glutamine_ext':1.0},'inconds':{}}
post = {'parameters':{'Carbon':0.0,'ATP':0.0,'Glutamine_ext':1.0},'inconds':{}}
DS = sim.createModelObject(Model)
DS_mut = sim.createModelObject(Model)
DS_mut.set(pars={'w_torc_snf':0.0})

DS.set(pars=pre['parameters'],
       ics=pre['inconds'])
P = sim.get_ss(DS)

DS.set(pars=post['parameters'],
       ics=P,
       tdata=[0,30])
P = sim.simulateModel(DS)
ax[0].plot(P['t'], P['Sch9'],'k--',lw=2,label='Snf1 -| TORC1')

#####
DS_mut.set(pars=pre['parameters'],
       ics=pre['inconds'])
P = sim.get_ss(DS_mut)
DS_mut.set(pars=post['parameters'],
       ics=P,
       tdata=[0,30])
P = sim.simulateModel(DS_mut)
ax[0].plot(P['t'], P['Sch9'],'#777777',lw=2,label='no interaction')
ax[0].plot([0.0,2.5,5.0,15.0,30.0], [1.0,0.34,0.13,0.15,0.20],'k.',markersize=8)

####################
del DS_mut
del DS
pre = {'parameters':{'Carbon':0.0,'ATP':0.0,'Glutamine_ext':1.0},'inconds':{}}
post = {'parameters':{'Carbon':1.0,'ATP':1.0,'Glutamine_ext':1.0},'inconds':{}}
DS = sim.createModelObject(Model)
DS_mut = sim.createModelObject(Model)
DS_mut.set(pars={'w_torc_snf':0.0})

DS.set(pars=pre['parameters'],
       ics=pre['inconds'])
P = sim.get_ss(DS)

DS.set(pars=post['parameters'],
       ics=P,
       tdata=[0,30])
P = sim.simulateModel(DS)
ax[1].plot(P['t'], P['Sch9'],'k--',lw=2)

#####
DS_mut.set(pars=pre['parameters'],
       ics=pre['inconds'])
P = sim.get_ss(DS_mut)
DS_mut.set(pars=post['parameters'],
       ics=P,
       tdata=[0,30])
P = sim.simulateModel(DS_mut)
ax[1].plot(P['t'], P['Sch9'],'#777777',lw=2)
ax[1].plot([0.0,2.0,5.0,15.0,30.0],[0.20,0.63,0.92,0.91,1.0],'k.',markersize=8)

ax[0].set_ylabel('%Sch9-P')
ax[0].legend()
ax[1].set_ylabel('%Sch9-P')
#ax[0].xlabel('min')
ax[0].set_title('Acute Glucose Starvation')
ax[1].set_xlabel('min')
ax[1].set_title('Relief from Glucose Starvation')
plt.suptitle('Snf1 mediated TORC1 inhibition')
#plt.tight_layout()
#plt.show()
plt.savefig(HOME + '/group/amogh-jalihal/papers/Nutrient-Signaling-Model/figs/Snf1-TORC1.png', dpi=300)
#return './figs/Snf1-TORC1.png'
