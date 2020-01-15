import pandas as pd
import matplotlib
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import sys
font = {'size'   : 17}

matplotlib.rc('font', **font)

# path = 'output/summary_24066.csv'
#path = 'output/combined_summary_24066.csv'
path = 'output/combined_summary_24066_truncated.csv'
df = pd.read_csv(path,index_col=0)
df['relcost'] = df['cost']/df.cost.min()
df = df[df['relcost'] < 2.0]
print(df.shape)

f, ax = plt.subplots(1,1,figsize=(8,4))
cols = [c for c in df.columns if c not in  ['cost','relcost']]
print(len(cols))
df = df.assign(sum=pd.Series(df[cols].sum(axis=1)))
df['explained'] = df['sum'] #100.*df['sum']/(len(cols))
#print([d*float(len(cols))/100 for d in df['explained'].unique()])

df = df[(df['explained'] < 29) & (df['explained'] > 22)]
minc = df[df['cost'] == min(df['cost'])]
# If using output/summary_24066.csv, use the following to filter out outliers
#df = df[(df['explained'] < 69.9) & (df['explained'] > 57.0)]

sns.boxplot(x='explained',y='relcost',data=df,ax=ax,showfliers=False,palette='coolwarm')
sns.stripplot(data=df,x='explained',y='relcost',ax=ax,jitter=0.4,alpha=0.05,size=2,color='k')
ax.plot(4, minc['relcost'],'ro',ms=7)
#ax.set_xticks([round(ne,1) for ne in sorted(df['explained'].unique())])
ax.set_xticklabels([round(ne,1) for ne in sorted(df['explained'].unique())])

ax.set_ylabel('x $C_{min}$')
ax.set_xlabel('Number of Experiments Explained (total 41)')

plt.tight_layout()
plt.savefig('cost-vs-explanatory-capacity.pdf',dpi=400)


## Subplots, showing cost scatter and box plot distributions
# f, ax = plt.subplots(2,1,figsize=(8,8))
# ax[0].axhline(100.*32./48.,c='r',zorder=1)
# ax[0].scatter(df['cost']/df['cost'].min(), 100.*df['sum']/48., c='k',s=2.0,zorder=2)

# ax[0].set_title('Percent Experiments Explained ({} parameter sets)'.format(df.shape[0]))
# ax[0].set_xlabel('x $C_{min}$')
# ax[0].set_ylabel('% experiments explained')

# percentexplained = sorted(df['sum'].unique())

# data = []
# for pe in percentexplained:
#     data.append(df[df['sum'] == pe]['cost']/df['cost'].min())
# #print(df[df['sum'] == 22])
    
# #ax[1].axvline(100.*32./48.,c='r',zorder=1)    
# ax[1].boxplot( data,zorder=2)
# ax[1].set_xticklabels([str(round(100.*pe/48.,1)) for pe in percentexplained])
# ax[1].set_ylabel('x $C_{min}$')
# ax[1].set_xlabel('% experiments explained')

# plt.tight_layout()
# plt.savefig('fits-vis-full-cutoff.png')

# numberexplained = sorted(df['sum'].unique())

# data = []
# for ne in numberexplained:
#     if 100*ne/48. < 57. or 100*ne/48.> 69.:
#         continue
#     else:
#         data.append(df[df['sum'] == ne]['cost']/df['cost'].min())
#print(df[df['sum'] == 22])
    
#ax.boxplot( data,showfliers=False )
