import pandas as pd

path = ['/data/amogh-jalihal/Nutrient-Signaling/automatic-iterations/full_paramlist_abs_eigs/iter4/expansion_iter4.txt']
outpath = '/home/jamogh/jalihal_projects/Research/nutrient-signaling/data/'
numSamples = 1000

for p in path:
    df = pd.read_csv(p, sep='\t')
    print(df.cost.min())
    outfname = p.split('/')[-1].split('.')[0]
    # dfsorted = df.sort_values(by=['cost'])[0:10000]
    # dfsample = dfsorted.sample(n=numSamples)
    dfsample = df.sample(n=numSamples)
    dfsample.to_csv(outpath + outfname + '_'+str(numSamples)+'_full.csv')
