from nutrient_signaling.simulators import get_simulator
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')
#modelpath = home + '/jalihal_projects/Research/nutrient-signaling/data/2019-02-26'
modelpath = home + '/jalihal_projects/Research/nutrient-signaling/data/2019-06-21'

f, ax = plt.subplots(1,2)
# wt -3AT
model = get_simulator(modelpath=modelpath,
                      simulator='py')


model.set_attr(pars={'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0},tdata=[0,90])
model.set_attr(pars={'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0}, ics = model.get_ss(),
               tdata=[0,90])
P = model.simulate_and_get_points()
eIF = [(1-e )/(e + 1e-3) for e in P['eIF']]
ax[0].plot(P['t'], eIF,'k--',lw=3,label='wt')
# wt -3AT
model = get_simulator(modelpath=modelpath,
                      simulator='py')


model.set_attr(pars={'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0},tdata=[0,90])
model.set_attr(pars={'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':0.0},
               ics = model.get_ss(),
               tdata=[0,90])
P = model.simulate_and_get_points()
eIF = [(1-e )/(e + 1e-3) for e in P['eIF']]
ax[1].plot(P['t'], eIF,'k--',lw=3,label='wt')

# S577A
model = get_simulator(modelpath=modelpath,
                      simulator='py')


model.set_attr(pars={'w_gcn_torc':0.0,'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0},tdata=[0,90])
model.set_attr(pars={'w_gcn_torc':0.0,'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0}, ics = model.get_ss(),
               tdata=[0,90])

P = model.simulate_and_get_points()
eIF = [(1-e )/(e + 1e-3) for e in P['eIF']]
ax[0].plot(P['t'], eIF,'r',label='S577A')

# S577A + 3AT
model = get_simulator(modelpath=modelpath,
                      simulator='py')


model.set_attr(pars={'w_gcn_torc':0.0,'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0},tdata=[0,90])
model.set_attr(pars={'w_gcn_torc':0.0,'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':0.0}, ics = model.get_ss(),
               tdata=[0,90])

P = model.simulate_and_get_points()
eIF = [(1-e )/(e + 1e-3) for e in P['eIF']]
ax[1].plot(P['t'], eIF,'r',label='S577A')

# snf1 delta - SC
model = get_simulator(modelpath=modelpath,
                      simulator='py')


model.set_attr(pars={'Snf1_T':0.0,
                     'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0},
               ics={'Snf1':0},
               tdata=[0,90])
model.set_attr(pars={'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0},
               ics = model.get_ss(),
               tdata=[0,90])

P = model.simulate_and_get_points()
eIF = [(1-e )/(e + 1e-3) for e in P['eIF']]
ax[0].plot(P['t'], eIF,'b',label='snf1')

# snf1 delta - 3AT
model = get_simulator(modelpath=modelpath,
                      simulator='py')

model.set_attr(pars={'Snf1_T':0.0,
                     'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':1.0},
               ics={'Snf1':0},
               tdata=[0,90])
model.set_attr(pars={'Carbon':1.0,'ATP':1.0,
                     'Glutamine_ext':0.0},
               ics = model.get_ss(),
               tdata=[0,90])

P = model.simulate_and_get_points()
eIF = [(1-e )/(e + 1e-3) for e in P['eIF']]
ax[1].plot(P['t'], eIF,'b',label='snf1')

ax[0].set_title('SC')
ax[1].set_title('SC+3AT')
ax[1].legend()
ax[0].set_ylabel('eIF-P/eIF')
ax[0].set_xlabel('time (minutes)')
ax[1].set_xlabel('time (minutes)')
plt.legend()
#plt.ylim([0,1])
plt.savefig('./img/claim1.png')
