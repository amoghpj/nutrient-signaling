import matplotlib
import matplotlib.pyplot as plt
from nutrient_signaling.timecourse import TimeCourse
font = {'size':24}
matplotlib.rc('font', **font)

datapath = "data/yaml/time-course-data.yaml"
alternatePars = "data/parameter-set-ensembles/expansion_iter4.txt"
modelpath = "data/2019-06-21"
comparisontype = "time"
simulator = "cpp"

# Initialize TimeCourse Comparison object
tccomp = TimeCourse()
# Set simulator type.
tccomp.setSimulator(modelpath, simulatortype='cpp')

print("Reading experimental data")
tccomp.readData(datapath)

print("Loading parameter sets")
tccomp.loadAlternateParameterSets(alternatePars)
tccomp.setNumberOfParameterSets(100)
print(tccomp.data[0])

#f = plt.figure(figsize=(,10))
order = [
    #1, 2, 3, 4, 5, 6, 7,
         0,
    #8, 9
]
titles = ["Glucose Addition",
          "Glucose Addition - $sch9\Delta$",
          "High Glutamine",
          "Low Glutamine",
          "High Glutamine - $gtr1\Delta$",
          "Glucose Starvation",
          "Glucose Addition",
          "Glucose Starvation",
          "Glucose Addition",
          "Rapamycin Treatment"]

ylabels = ["cAMP", "cAMP",
           "Sch9-P","Sch9-P",
           "Sch9-P","Sch9-P",
           "Sch9-P","Snf1-P",
           "log(nMig1/cMig1)", "Rel. RPL32 mRNA"]

for i, simid in enumerate(order):
    print(simid)
    f, ax = plt.subplots(1,1,figsize=(6,4))
    experiment = tccomp.data[simid]
    tccomp.compareToExperiment(experiment, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    #ax.set_ylabel(ylabels[i])
    ax.set_title("")
    #ax.set_title(titles[i])
    plt.tight_layout()
    name = experiment['description']
    plt.savefig(f'img/{name}.pdf')
    plt.close()    

