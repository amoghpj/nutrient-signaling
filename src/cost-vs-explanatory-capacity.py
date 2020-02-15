import pandas as pd
import matplotlib
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import sys
font = {'size'   : 17}

matplotlib.rc('font', **font)

rapamycinExps = ['6-gln3 gat1',
                 '7-2\mu URE2',
                 '8-wt',
                 '9-gln3',
                 '10-gat1',
                 '14-bcy1',
                 '15-ira1',
                 '16-ira1 ira2',
                 '17-ras2',
                 '18-tpk1',
                 '19-RAS2v19 gln3 gat1',
                 '20-TPK1 gln3 gat1',
                 '21-bcy1',
                 '22-bcy1 gln3 gat1',
                 '27-gln3 gcn4',
                 '28-gcn4',
                 '34-snf1',
                 '35-reg1',
                 '36-ure2',
                 '37-tap42-11',
                 '38-SCH9^{DE}',
                 '39-SCH9^{DE} tap42',
                 ]


## Plot non rapamycin experiments
path = "output/summary_24066_gln3_tap42_redefined.csv"
#path = "output/summary_24066_rap-is-dot6.csv"

rapamycinExps.append('25-snf1')
df = pd.read_csv(path,index_col=0)
df['relcost'] = df['cost']/df.cost.min()
df = df[df['relcost'] < 2.0]
df = df.drop(columns=rapamycinExps)

f, ax = plt.subplots(1,1,figsize=(8,4))
cols = [c for c in df.columns if c not in  ['cost','relcost']]
print(len(cols))
df = df.assign(sum=pd.Series(df[cols].sum(axis=1)))
df['explained'] = df['sum'] #100.*df['sum']/(len(cols))

minc = df[df['cost'] == min(df['cost'])]
print(minc.explained)

# If using output/summary_24066.csv, use the following to filter out outliers
# df = df[(df['explained'] < 69.9) & (df['explained'] > 57.0)]

sns.boxplot(x='explained',y='relcost',data=df,ax=ax,showfliers=False,palette='coolwarm')
sns.stripplot(data=df,x='explained',y='relcost',ax=ax,jitter=0.4,alpha=0.05,size=2,color='k')
# plot the mincost point
#ax.plot(minc['explained'] - df.explained.min(), minc['relcost'],'ro',ms=7)
ax.set_xticklabels([round(ne,1) for ne in sorted(df['explained'].unique())])

ax.set_ylabel('x $C_{min}$')
ax.set_xlabel(f'Number of Experiments Explained (total {len(cols)})')

plt.tight_layout()
plt.savefig('cost-vs-explanatory-capacity-drop-rap-gln3-drop25snf1.pdf',dpi=400)
# #plt.savefig('cost-vs-explanatory-capacity-drop-rap-dot6.pdf',dpi=400)


## Plot rapamycin experiments for both definitions
rapamycinExps.append('cost')
# glngcnpath = "output/summary_24066_gln3_tap42_redefined.csv"
# dotpath = "output/summary_24066_rap-is-dot6.csv"

# glngcndf = pd.read_csv(glngcnpath,index_col=0)
# dropcols  = [c for c in glngcndf.columns if c not in rapamycinExps]
# glngcndf['relcost'] = glngcndf['cost']/glngcndf.cost.min()
# glngcndf = glngcndf[glngcndf['relcost'] < 2.0]
# glngcndf = glngcndf.drop(columns=dropcols)


# print(dotpath)
# dotdf = pd.read_csv(dotpath,index_col=0)
# dropcols =  [c for c in dotdf.columns if c not in rapamycinExps]
# dotdf['relcost'] = dotdf['cost']/dotdf.cost.min()
# dotdf = dotdf[dotdf['relcost'] < 2.0]
# dotdf = dotdf.drop(columns=dropcols)

# f, ax = plt.subplots(2,1,figsize=(8,8))
# cols = [c for c in dotdf.columns if c not in  ['cost','relcost']]

# dotdf = dotdf.assign(sum=pd.Series(dotdf[cols].sum(axis=1)))
# dotdf['explained'] = dotdf['sum']

# glngcndf = glngcndf.assign(sum=pd.Series(glngcndf[cols].sum(axis=1)))
# glngcndf['explained'] = glngcndf['sum']

# # minc = df[df['cost'] == min(df['cost'])]
# # print(minc.explained)

# # if using output/summary_24066.csv, use the following to filter out outliers
# # df = df[(df['explained'] < 69.9) & (df['explained'] > 57.0)]

# for df,axis,title in zip([dotdf, glngcndf],ax,['Dot6', 'Gln3 Gcn4']):
#     sns.boxplot(x='explained',y='relcost',data=df,ax=axis,showfliers=False,palette='coolwarm')
#     sns.stripplot(data=df,x='explained',y='relcost',ax=axis,jitter=0.4,alpha=0.05,size=2,color='k')
#     axis.set_xticklabels([round(ne,1) for ne in sorted(df['explained'].unique())])
#     axis.set_ylabel('x $C_{min}$')
#     axis.set_title(title)
#     axis.set_xlabel("")
# axis.set_xlabel(f'Number of Experiments Explained (total {len(cols)})')

# plt.tight_layout()
# plt.savefig('cost-vs-explanatory-capacity-rapamycin-exps.pdf',dpi=400)

