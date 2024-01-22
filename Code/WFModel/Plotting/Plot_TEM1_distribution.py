# Import necessary libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define font sizes for the plot
SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

# Configure the plot font sizes using matplotlib settings
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title


data=pd.read_csv("dataEmpirical/data_TEM1_Codon_Substitutions.csv")
dataCrossCodons=pd.read_csv("dataEmpirical/data_TEM1_CrossCodonsDoubles.csv")


count = ((data['mutation'] == 'single') ).sum()
print("AA accessible via SN", count)
count = ((data['mutation'] == 'single') & (data['fitness'] == 'unknown')).sum()
print("Missing fitness values", count)
count = ((data['mutation'] == 'double') ).sum()

print("One amino acid changes only accessible by doubles", count)
count = ((data['mutation'] == 'double') & (data['fitness'] == 'unknown')).sum()
print("Missing fitness values", count)

print("Two amino acid changes only accessible by doubles", len(dataCrossCodons))
count = ((dataCrossCodons['mutation'] == 'double') & (dataCrossCodons['fitness'] == 'unknown')).sum()
print("Missing fitness values", count)

# Remove rows with NA values in the mutation column
data = data.dropna(subset=['mutation'])

# Combine the two datasets into a single DataFrame
data = pd.concat([data, dataCrossCodons], ignore_index=True)

# Filter out rows with 'unknown' fitness values
data = data[data['fitness'] != 'unknown']

data['fitness'] = data['fitness'].astype(float)


#####################################################
#####################################################

# Rename the 'mutation' column to 'type'
data.rename(columns={'mutation':'type'}, inplace = True)

# Replace 'tandem' values with 'double' in the 'type' column
data['type'] = data['type'].replace({'tandem': 'double'})


# create a figure with two subplots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, figsize=(5, 3.7), gridspec_kw={'height_ratios': [3, 0.18,1,1], "hspace":0},)

# create a histogram using seaborn
sns.histplot(data=data, x="fitness", hue="type", bins=40, palette=["#62AAC5", "#DD903B"], log_scale=(False, True),hue_order=['single', 'double'],stat='probability', alpha=0.7, common_norm=False, legend=True, ax=ax1)
ax1.set_xlabel('')
ax1.set_ylabel('Frequency')
ax1.set_title('TEM1')

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

plt.savefig("TEM1_Dist.pdf", bbox_inches = "tight")

plt.show()


