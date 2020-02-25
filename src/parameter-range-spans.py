"""
Generate Parameter Range plots
1. Sort representative set of values,
2. obtain min and max for each parameter 
3. Calculate fold change wrt min
4. Plot log 10 fold change
"""
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
font = {'size': 13}
matplotlib.rc('font', **font)

path = 'data/parameter-set-ensembles/expansion_iter4.txt'
df = pd.read_csv(path,index_col=False,sep='\t')

df['relcost'] = df['cost']/df.cost.min()
df = df[df['relcost'] < 2.0]

dfmin = df[df.cost == df.cost.min()]
dfmin.reset_index(inplace=True)
parameters = [par for par in df.columns if par not in ['c_t', 'c_p','relcost','cost']]

ranges = {par:{} for par in parameters}

for par in parameters:
    mi = float(df[par].min())/dfmin[par][0]
    mx = float(df[par].max())/dfmin[par][0]
    lmi = np.log10(mi)
    lmx = np.log10(mx)
    ranges[par]['min'] = lmi
    ranges[par]['max'] = lmx
    ranges[par]['range'] = lmx - lmi

rangedf = pd.DataFrame(ranges)
rangedf = rangedf.T
rangedf = rangedf.sort_values(by='range',
                              axis='index')

counter = 0
f,ax = plt.subplots(1,1, figsize=(7,16))

for par, row in rangedf.iterrows():
    r = Rectangle((row['min'],counter-0.5),
                  width=row['range'],facecolor='#000000',
                  height=1.0)
    ax.add_patch(r)
    counter +=1

print(rangedf['range'].min())
print(rangedf['range'].max())
ax.set_yticks([i for i in range(rangedf.shape[0])])        
ax.set_yticklabels(list(rangedf.index))    
plt.xlim([-1.3,1.3])
plt.ylim([-2, 82])
rangedf.to_csv('test.csv', sep='\t')
plt.xlabel('log$_{10}(p/p_{min})$')
plt.tight_layout()
plt.savefig('img/parameter-ranges-representative-set.png', dpi=300)
