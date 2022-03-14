import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 12})

def read_pickle_file(f):
    res=pickle.load(open(f,'rb'))

    return res

folder='/Users/vpandey/projects/githubs/Fertility_screen/ErrorPickle/'
pilot=read_pickle_file(folder+'pilot.p')
pool1=read_pickle_file(folder+'pool1.p')
pool2=read_pickle_file(folder+'pool2.p')
pool3=read_pickle_file(folder+'pool3.p')
pool4=read_pickle_file(folder+'pool4.p')
pool6=read_pickle_file(folder+'pool6.p')
pool5_1=read_pickle_file(folder+'pool5_1.p')
pool7_1=read_pickle_file(folder+'pool7_1.p')
pool5_2=read_pickle_file(folder+'pool5_2.p')
pool7_2=read_pickle_file(folder+'pool7_2.p')

all_data=[pilot,pool1,pool2,pool3,pool4,pool6,pool5_1,pool7_1,pool5_2,pool7_1]\


list_data=[[],[],[],[],[],[],[],[]]
list_labels=['mf1_day0_mean','mf1_day0_relErr','mf2_day0_mean','mf2_day0_relErr','mf1_day13_mean','mf1_day13_relErr','mf2_day13_mean','mf2_day13_relErr']


for item in all_data:
    data_day0=item[0]
    data_day13=item[1]

    for i in range(4):
        list_data[i].append(data_day0[i])
        list_data[4+i].append(data_day13[i])

## concat the data_frames
concat_data=[]
for item in list_data:
    result = pd.concat(item)
    concat_data.append(result)

## combine mosquito feed
sex=['GCKO2','145480']



fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(16, 16))
fig.tight_layout(pad=3.0)
feed_day=['mf1_day0_','mf2_day0_','mf1_day13_','mf2_day13_']
index=[0,2,4,6]

for counter,i in enumerate(index):
    ## scatter plot
    axes[counter, 0].plot(concat_data[i].iloc[:,0].values, concat_data[i+1].iloc[:,0].values, 'b.')
    axes[counter, 0].set_xlim([-20, 0])
    # axes[counter, 0].set_ylim([-20 ,0])
    axes[counter, 0].set_title(feed_day[counter]+sex[0])
    axes[counter, 0].set_xlabel('log2(rel. abundance)')
    axes[counter, 0].set_ylabel('Relative error')


    axes[counter, 1].hist(concat_data[i].iloc[:,0].values, 20, density=True,color= 'b')
    axes[counter, 1].set_xlim([-20, 0])
    axes[counter, 1].set_title(feed_day[counter]+sex[0])
    axes[counter, 1].set_xlabel('log2(rel. abundance)')
    axes[counter, 1].set_ylabel('Frequency')

    # axes[counter, 0].set_ylim([-20, 0])


    axes[counter, 2].plot(concat_data[i].iloc[:,1].values, concat_data[i+1].iloc[:,1].values, 'b.')
    axes[counter, 2].set_xlim([-20 ,0])
    axes[counter, 2].set_title(feed_day[counter]+sex[1])
    axes[counter, 2].set_xlabel('log2(rel. abundance)')
    axes[counter, 2].set_ylabel('Relative error')

    # axes[counter, 0].set_ylim([-20, 0])

    axes[counter, 3].hist(concat_data[i].iloc[:,1].values, 20, density=True,color= 'b')
    axes[counter, 3].set_xlim([-20, 0])
    axes[counter, 3].set_title(feed_day[counter]+sex[1])
    axes[counter, 3].set_xlabel('log2(rel. abundance)')
    axes[counter, 3].set_ylabel('Frequency')
    # axes[counter, 0].set_ylim([-20, 0])
plt.savefig('error_analysis.pdf')












## mf1
