import pandas as pd
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
'35-reg1',
'36-ure2',
'37-tap42-11',
'38-SCH9^{DE}',
'39-SCH9^{DE} tap42']

path = 'output/summary_24066' + suffix + ".csv"

df = pd.read_csv(path)
df.drop(columns=rapamycinExps,inplace=True)
out = {}
for col in df.columns:
    print(col, 100*df[col].sum()/df.shape[0])
    out[col] = [100*df[col].sum()/df.shape[0]]

outdf = pd.DataFrame(out)

outdf = outdf.T
outdf.columns = ['Confidence']
outdf = outdf.sort_values(by='Confidence')
# outdf.to_csv('confidence_scores.csv')
# outdf.to_csv('confidence_scores-with-rtg-truncated.csv')
outdf.to_csv('confidence_scores'+suffix+'-drop.csv')
