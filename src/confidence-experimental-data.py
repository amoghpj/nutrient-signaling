import pandas as pd
import sys
#path = 'output/combined_summary_24066.csv'
# path = 'output/combined_summary_24066_truncated.csv'
# path = 'output/summary-24066-new-rapamycin.csv'
#suffix = "_rap-is-dot6"
#suffix = "_rap-is-gln3-gcn4"
suffix = "_gln3_tap42_redefined"

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
                 '25-snf1',
                 '31-mig1 snf1 pde2',
                 '32-mig1 snf1 pde2',
                 '33-mig1 snf1 pde2',
                 '35-reg1',
                 '36-ure2',
                 '37-tap42-11',
                 '38-SCH9^{DE}',
                 '39-SCH9^{DE} tap42','cost']

# path = 'output/summary_24066' + suffix + ".csv"
path = 'search-redone-sample/summary_43982_test.csv'

df = pd.read_csv(path, index_col=0)
df.drop(columns=rapamycinExps,inplace=True)
out = {}
for col in df.columns:
    #print(col, 100*df[col].sum()/df.shape[0])
    out[col] = [100*df[col].sum()/df.shape[0]]

outdf = pd.DataFrame(out)

outdf = outdf.T
outdf.columns = ['Confidence']
outdf = outdf.sort_values(by='Confidence', ascending=False)
# outdf.to_csv('confidence_scores.csv')
# outdf.to_csv('confidence_scores-with-rtg-truncated.csv')
# outdf.to_csv('confidence_scores'+suffix+'-drop.csv')
outname = 'confidence_scores-rebuttal.csv'
outdf.to_csv(outname)
with open(outname.replace('.csv', '.tex'), 'w') as outfile:
    # two columns
    full = outdf.shape[0]
    if full % 2 == 0:
        half = full/2
    else:
        half = int(full/2) +1
    for i in range(half):
        exp1 = outdf.index[i]
        val1 = outdf.iloc[i]['Confidence']
        if half + i < full:
            exp2 = outdf.index[half + i]
            val2 = outdf.iloc[half + i]['Confidence']
            outfile.write(f"%s & %.2f & %s & %.2f \\\\\n" % (exp1, val1, exp2, val2))
        else:
            outfile.write(f"%s & %.2f & & \\\\\n" % (exp1, val1))


