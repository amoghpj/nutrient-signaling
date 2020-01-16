import pandas as pd
#path = 'output/combined_summary_24066.csv'
path = 'output/combined_summary_24066_truncated.csv'

df = pd.read_csv(path)
out = {}
for col in df.columns:
    print(col, 100*df[col].sum()/df.shape[0])
    out[col] = [100*df[col].sum()/df.shape[0]]

outdf = pd.DataFrame(out)

outdf = outdf.T
outdf.columns = ['Confidence']
outdf = outdf.sort_values(by='Confidence')
# outdf.to_csv('confidence_scores.csv')
outdf.to_csv('confidence_scores-with-rtg-truncated.csv')
