import pickle
import pandas as pd

prev_to_new=pickle.load(open('/Users/vikash/git_hub/Fertility_screen/data/prevTonew_PBANKA.pickle','rb'))
kk=prev_to_new
df=pd.DataFrame(prev_to_new.items(), columns=['old_id', 'new_id'])
df.to_csv('old_to_new_ids.txt','\t')
