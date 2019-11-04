import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
df = pd.read_csv('summary_1000_of_full.csv',index_col=0)
cols = [c for c in df.columns if c != 'cost']
df['sum'] = df[cols].sum(axis=1)

f, ax = plt.subplots(2,1,figsize=(8,8))
ax[0].axhline(100.*32./48.,c='r',zorder=1)
ax[0].scatter(df['cost']/df['cost'].min(), 100.*df['sum']/48., c='k',s=2.0,zorder=2)

ax[0].set_title('Percent Experiments Explained ({} parameter sets)'.format(df.shape[0]))
ax[0].set_xlabel('x $C_{min}$')
ax[0].set_ylabel('% experiments explained')

percentexplained = sorted(df['sum'].unique())
print(percentexplained)
data = []
for pe in percentexplained:
    data.append(df[df['sum'] == pe]['cost']/df['cost'].min())
ax[1].axvline(100.*32./48.,c='r',zorder=1)    
ax[1].boxplot( data,zorder=2)
ax[1].set_xticklabels([str(round(100.*pe/48.,1)) for pe in percentexplained])
ax[1].set_ylabel('x $C_{min}$')
ax[1].set_xlabel('% experiments explained')

plt.tight_layout()
plt.savefig('fits-vis-1000.png')

