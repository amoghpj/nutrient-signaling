import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
path  = 'output/confidence_scores.csv'

df = pd.read_csv(path)
print(df.Confidence)
plt.hist(df.Confidence, bins=20)
plt.xlabel('% of parameter sets agreeing with experiment')
plt.ylabel('Number of Experiments')
plt.title('Summary of model agreement (48 experiments)')

plt.savefig('img/prediction-confidence-histogram.png')


