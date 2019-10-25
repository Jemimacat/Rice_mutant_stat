
# coding: utf-8

# # Stat Plots for GATK
# 2019-10-25

# In[9]:


import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[11]:


chem = pd.read_csv('snp.chem.stat.num.tsv',sep='\t')
chem_total = np.sum(chem.VARIANCE_NUM)
chem_total


# In[6]:


nature = pd.read_csv('snp.nature.stat.num.tsv',sep='\t')
nature_total = np.sum(nature.VARIANCE_NUM)
nature_total


# ## Chemical Induction & Natural Status (SNPs)

# In[12]:


names = ['Chemical','Natural']
counts = [chem_total,nature_total]
y_pos = np.arange(len(names))

plt.bar(y_pos, counts, align='center', alpha=0.5)
plt.xticks(y_pos, names)
plt.ylabel('SNPs')
plt.title('Chemical Induction vs. Natural Status')

plt.show()


# In[24]:


## fetch mutation numbers
def mutation_type_num_dict(df):
    change_byte = list(df.CHANGE_BYTE)
    variance_num = list(df.VARIANCE_NUM)
    num_dict = dict(zip(change_byte,variance_num))
    return num_dict

chem_dict = mutation_type_num_dict(chem)
nature_dict = mutation_type_num_dict(nature)

mtype = sorted(chem_dict.keys()) ## Mutation type
chem_nums = [] ## chemical induced mutation numbers
nature_nums = [] ## natural mutation numbers
for mt in mtype:
    chem_nums.append(chem_dict[mt])
    if mt in nature_dict:
        nature_nums.append(nature_dict[mt])
    else:
        nature_nums.append(0)

# create plot
fig, ax = plt.subplots()
index = np.arange(len(mtype))
bar_width = 0.35
opacity = 0.8

rects1 = plt.bar(index, nature_nums, bar_width,
alpha=opacity,
color='b',
label='Natural')

rects2 = plt.bar(index + bar_width, chem_nums, bar_width,
alpha=opacity,
color='g',
label='Chemical')

plt.ylabel('SNPs')
plt.title('Chemical Induction vs. Natural Status')
plt.xticks(index + bar_width, mtype)
plt.legend()

plt.tight_layout()
plt.show()






