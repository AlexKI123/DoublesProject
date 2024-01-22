import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title


data=pd.read_csv("dataEmpirical/data_P53.txt", sep=",")


data = data.dropna(subset=['mutation'])


#####################################################
#####################################################

singles_df = data[data['mutation'] == "single"]
singles = singles_df['fitness'].to_numpy()

singles_count = singles_df['GENIE_Mutation_Counts'].to_numpy()
print("singles_count", sum(singles_count))

#####################################################
#####################################################

doubles_df = data[data['mutation'] == "tandem"]
doubles = doubles_df['fitness'].to_numpy()

doubles_count = doubles_df['GENIE_Mutation_Counts'].to_numpy()
print("doubles_count", sum(doubles_count))

print(sum(doubles_count)/(sum(doubles_count)+sum(singles_count)))


#np.save('singles_P53.npy', singles)
#np.save('doubles_P53.npy', doubles)

print("len(singles)", len(singles))
print("len(doubles)", len(doubles))

print("np.mean(singles)", np.mean(singles))
print("np.mean(doubles)", np.mean(doubles))

#print("t-test", stats.ttest_ind(singles,doubles, equal_var=True, alternative="greater", trim=0.1))
#print("t-test", stats.ttest_ind(singles,doubles, equal_var=False, alternative="greater", trim=0.1))



data.rename(columns = {'mutation':'type'}, inplace = True)
data['type'] = data['type'].replace({'tandem': 'double'})


print(data)


# create a figure with two subplots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, figsize=(5, 3.7), gridspec_kw={'height_ratios': [3, 0.18,1,1], "hspace":0},)

# create a histogram using seaborn
sns.histplot(data=data, x="fitness", hue="type", bins=40, palette=["#62AAC5", "#DD903B"], log_scale=(False, True),hue_order=['single', 'double'],stat='probability', alpha=0.7, common_norm=False, legend=False, ax=ax1)
ax1.set_xlabel('')
ax1.set_ylabel('Frequency')
ax1.set_title('P53')

ax2.axis("off")

sns.boxplot(data=data, x="fitness", y="type", palette=["#62AAC5", "#DD903B"], order=['single', 'double'], showmeans=True, meanprops={"markerfacecolor":"black", "markeredgecolor":"black"}, ax=ax3,)

xmin, xmax = ax1.get_xlim()
ax3.set_xlim([xmin, xmax])
ax3.set_xlabel('')
ax3.set_ylabel('All')


data = data[data['fitness'] > 1]
sns.boxplot(data=data, x="fitness", y="type", palette=["#62AAC5", "#DD903B"], order=['single', 'double'], showmeans=True, meanprops={"markerfacecolor":"black", "markeredgecolor":"black"}, ax=ax4,)
xmin, xmax = ax1.get_xlim()
ax4.set_xlim([xmin, xmax])
ax4.set_xlabel('Fitness')
ax4.set_ylabel('Benef.')


plt.savefig("plotP53.pdf", bbox_inches = "tight")
plt.show()
