import pandas as pd
import matplotlib.pyplot as plt
df=pd.read_csv('/Users/vikash/git-hub/Fertility_screen/data/wt21.txt',sep='\t')
tmp=df[['GCKO2_RGR','GCKO2_diff_min','g145480_RGR','g145480_diff_min']]
tmp.hist(bins=5)
plt.show()
import pdb; pdb.set_trace()
