from nutrient_signaling.simulators import get_simulator
import matplotlib.pyplot as plt
import os
home = os.path.expanduser('~')
modelpath = home + '/jalihal_projects/Research/nutrient-signaling/data/2019-02-26'

preshift = {'Glutamine_ext':1.0,'Carbon':1.0,'ATP':1.0}
treatment = {'Glutamine_ext':1.0,'Carbon':1.0,'ATP':1.0,'TORC1_T':0}
straindef = {'w_gln_sit':0.0}

plotspec = {'wt':{'c':'k','label':'wt'},
            'mutant':{'c':'r','label':'$sit4$'}}
timeAtTreatment = 60
timePostTreatment = 90

f, ax = plt.subplots(1,1,figsize=(3.5,3))
toplot = 'Gln3'

for strain in ['wt','mutant']:
    # Initialize
    model = get_simulator(modelpath=modelpath,
                          simulator='cpp',
                          execpath=modelpath + '/',
                          executable='nutsig.o')
    if strain == 'mutant':
        model.set_attr(pars = straindef)
        
    model.set_attr(pars=preshift,tdata=[0,90])
    # preshift
    model.set_attr(ics=model.get_ss(),tdata=[0,timeAtTreatment]) 
    Pre = model.simulate_and_get_points()
    # treatment, start from previous ss
    model.set_attr(ics=model.get_ss(),
               pars=treatment,tdata=[0,timePostTreatment])
    Post = model.simulate_and_get_points()
    
    # Plot
    plotcolor = plotspec[strain]['c']
    plotlabel = plotspec[strain]['label']
    posttime = [timeAtTreatment + t for t in Post['t']]
    ax.plot(Pre['t'], Pre[toplot],c=plotcolor,lw=3)
    ax.plot(posttime, Post[toplot],
            c=plotcolor,lw=3,
            label=plotlabel)
    ax.plot(Pre['t'][-1], Pre[toplot][-1],c=plotcolor,
            marker='o',ms=10,markerfacecolor="white")            
    ax.plot(posttime[-1],Post[toplot][-1],
            c=plotcolor,lw=3,
            marker='>',
            markerfacecolor="white",
            #markeredgecolor=plotcolor,
            ms=10)

drawuntil = 0.4
ax.set_yticks([0.0,0.5,1.0])
ax.set_xticks([0,30,60,90,120,150, 180])
ax.axvline(x=timeAtTreatment,ymin=drawuntil,ymax=0.9,linestyle='--',c='k')
ax.plot(60,drawuntil,'v',c='k')
ax.set_ylim([-0.1,1.1])
ax.set_ylabel('Gln3 activity')
ax.set_xlabel('Time (minutes)')
ax.legend(frameon=False)
ax.set_title('Rapamycin treatment')
plt.tight_layout()
plt.savefig('img/fig2b.pdf')
