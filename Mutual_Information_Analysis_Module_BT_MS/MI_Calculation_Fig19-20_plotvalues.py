#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import csv
import math
import collections
import natsort
import itertools
import json
import operator
from itertools import islice
import matplotlib.pyplot as plt
from collections import OrderedDict
from itertools import groupby
#import plotly.plotly as py
#import plotly.graph_objs as go

# os.system('clear')

######### BIOMASS MUTUAL INFORMATION CALCULATION #############

# -----Run Python Scripts from a Command Prompt---------
# !/usr/bin/env python
# $ chmod a+x hello.py
# $ ./BT_1.py

# load the files MS_7bits.csv or BT_7bits.csv with AllFBAs_7.csv media files with correct path. BT = B.theta MS=M.smithii
# carefully load the corresponding file for each organism
#reaction_states = pd.read_csv('/Users/Zarzeeth/Box/ZS_Research/MetabolicMC/MetabolicMC_Git/Mutual_Information_Analysis_Module/BT_7bits.csv')
reaction_states = pd.read_csv('/Users/Zarzeeth/Box/ZS_Research/MetabolicMC/MetabolicMC_Git/Mutual_Information_Analysis_Module/MS_7bits.csv')
reaction_states.as_matrix()

media_combinations = pd.read_csv('/Users/Zarzeeth/Box/ZS_Research/MetabolicMC/MetabolicMC_Git/Mutual_Information_Analysis_Module/AllFBAs_7.csv')
media_combinations.as_matrix()

# ----Input validation of Media/FBAs with Binary Matrix FBAs------
# 1.0 Number of rows in Media.csv file =  (Number of columns -1)
#   1.0. If they are different: Through an ERROR saying missed match number of FBAs in media and binary matrix. 
# 1.1 Check whether the elements in Media.csv file contains only binary values (i.e. 0 and 1)
#   1.1. If the elements are different: Through an ERROR saying not approapriate input values
# 1.2 Check whether the compounds in Media.csv file match with number of FBAs
#   1.2. If the compounds are different from number of FBAs: Through an ERROR saying not appropriate input values

s_df1 = reaction_states.shape
s_df2 = media_combinations.shape

Temp_df2 = np.array(media_combinations.values)
# Create matrix with only the elements remove first column and all the rows
Temp_df2 = Temp_df2[0:, 1:]

Bin_val_check = np.array_equal(Temp_df2, Temp_df2.astype(bool))
num_compounds = (s_df2[1]) - 1

if ((s_df1[1] - 1) != s_df2[0]) or (Bin_val_check != True) or (int(math.log(s_df2[0], 2)) != num_compounds):
    print('invalid input values')

# -----All possible combination of the chemical compounds----------------------
# 2.0 Sperating m0 from rest of the lables

Temp1_df2 = media_combinations

cols = Temp1_df2.columns
for i in range(1, len(cols)):
    Temp1_df2.loc[Temp1_df2[cols[i]] == 1, cols[i]] = cols[i]

# print Temp1_df2

# 2.1 Creating a disctionary for all FBAs except m0
# print len(Temp1_df2)
mydict = {}
for x in range(0, len(Temp1_df2)):
    for i in range(1, s_df2[1]):
        currentvalue = Temp1_df2.iloc[x, i]
        currentid = Temp1_df2.iloc[x, 0]
        currentvalue = Temp1_df2.iloc[x, i]
        mydict.setdefault(currentid, [])
        if float(currentvalue) > 0:
            mydict[currentid].append(currentvalue)

# Add the first key as m0
media_0_name = 'm0'
mydict[media_0_name] = "['0']"
# Sort the keys
mydict = collections.OrderedDict(natsort.natsorted(mydict.items()))
# print mydict

for k, v in mydict.items():
    print(k, v)

# List of Compounds combination in the list
my_combi_list = []
Compounds_Combi = list(range(1, num_compounds + 1))
for L in range(0, len(Compounds_Combi) + 1):
    for subset in itertools.combinations(Compounds_Combi, L):
        my_combi_list.append(list(subset))
# print my_combi_list


# Created a dictionary where the keys: 
# list of compounds combination 
# values are corresponding FBAs list in df2
result_dict = {}
for element in my_combi_list[1:]:
    for k, v in mydict.items():
        if set(v).issubset(set(map(lambda x: str(x), element))):
            key = ','.join(map(lambda x: str(x), element))
            if result_dict.get(key):
                media_list = result_dict[key]
                media_list.append(k)
                media_list = list(set(media_list))
                result_dict.update({key: media_list})
            else:
                result_dict.update({key: [media_0_name, k]})
# print result_dict

# Created a dictionary where the keys are: 
# list of compounds combination 
# values are compounds combination FBAs with df1 vaules 
All_Comp_Combi_dic = {}
for column, value in result_dict.items():
    All_Comp_Combi_dic.update({column: reaction_states.get(value)})

# To print an item from the All_Comp_Combi_dic
df = (pd.DataFrame(list(All_Comp_Combi_dic.items())))

# Load the file contain the information of FBAs(media) along with corresponding Biomass (growth)
#df3 = pd.read_csv('/Users/Zarzeeth/Box/ZS_Research/MetabolicMC/MetabolicMC_Git/Mutual_Information_Analysis_Module/BT_B_U_S_All.csv')
df3 = pd.read_csv('/Users/Zarzeeth/Box/ZS_Research/MetabolicMC/MetabolicMC_Git/Mutual_Information_Analysis_Module/MS_B_U_S_All.csv')

df3.as_matrix()


MI_dict_b_u_s = {}
for r in range(0, len(df[0])):

    reaction_states = df[1][r]

    def get_groups(flux_df):
        groups = collections.defaultdict(list)
        unique = flux_df.aggregate(lambda x: hash(str(x.values)))
        for k, v in unique[0:].iteritems():
             groups[v].append(k)
        return dict([(i, g) for i, g in enumerate(groups.values())])


    n_group = collections.defaultdict(int)
    groups = get_groups(reaction_states)
    for group in groups.values():
        n_group[len(group)] += 1
    # print n_group
    # print groups

    groups_count = {}
    for key, values in groups.items():
        groups_count[key] = len(values)
    # print groups_count

    # Take first FBA label of every group
    group_id = {}
    for k, v in groups.items():
        group_id.update({k: groups.values()[k][0]})

    # Obtain the Biomass of each Group
    cols_df = group_id.values()
    cols_df3 = df3.columns
    #print cols_df

    # Dictionary of first FBA label of every group and its corresponding number of members
    groups_label_count = {}
    for k, v in groups_count.items():
        groups_label_count.update({cols_df[k]: v})
    print groups_label_count

    # Extract FBA Groups biomass inside df3
    Groups_Biomass = df3[df3['FBAs'].isin(cols_df)].reset_index(drop=True)
    #print cols_df
    #print Groups_Biomass

    # Regroup based on the biomass values for B.theta
    #re_group = Groups_Biomass.groupby(['cpd00254_e0','cpd00067_e0','cpd00205_e0','cpd00034_e0','cpd00009_e0','cpd00013_e0','cpd10515_e0','cpd00099_e0','cpd00030_e0','cpd00012_e0','cpd00048_e0','cpd00149_e0','cpd00073_e0','cpd10516_e0','cpd00001_e0','cpd00011_e0','cpd00007_e0','cpd11416_c0','cpd00058_e0','cpd00084_e0','cpd00063_e0','cpd00027_e0','cpd00029_e0','cpd00028_e0','Biomass'])

    # Regroup based on the biomass values for M.smithii
    re_group = Groups_Biomass.groupby(
        ['cpd00129_e0', 'cpd00001_e0', 'cpd00058_e0', 'cpd00013_e0', 'cpd10515_e0', 'cpd00009_e0', 'cpd00067_e0',
         'cpd00254_e0', 'cpd00205_e0', 'cpd00034_e0', 'cpd00047_e0', 'cpd10516_e0', 'cpd00226_e0', 'cpd00007_e0',
         'cpd00011_e0', 'cpd00133_e0', 'cpd00028_e0', 'cpd00030_e0', 'cpd00048_e0', 'cpd00149_e0', 'cpd00099_e0',
         'cpd00063_e0', 'cpd11416_c0', 'cpd00084_e0', 'cpd00268_e0', 'cpd00239_e0', 'cpd00027_e0', 'cpd00075_e0',
         'cpd00029_e0', 'cpd01414_e0', 'cpd00119_e0', 'biomass'])

    my_list=[]
    for index,values in re_group:
        my_list.append(values['FBAs'].values)

    B_U_S_dict = {}
    for media in my_list:
        if len(media) > 1:
            media_cond = 0
            for i in (0,len(media)-1):
                media_cond += groups_label_count[media[i]]
            B_U_S_dict.update({str(media)[1:-1]: media_cond})
            #final_my_dict.update({tuple(media.tolist()):media_cond})
        else:
            B_U_S_dict.update({str(media)[1:-1]:groups_label_count[str(tuple(media.tolist()))[1:-1][:-1][1:-1]]})
    B_U_S_dict = {eval(k):v for k, v in B_U_S_dict.iteritems()}
    print B_U_S_dict
    Summery = pd.DataFrame(B_U_S_dict.items(), columns=['FBAs_x', 'FBAs_y'])


    Data_4_CondMI = Summery.groupby('FBAs_y').count()
    Data_4_CondMI = Data_4_CondMI.to_dict(orient='dict')

    print Data_4_CondMI

    for k, v in Data_4_CondMI.items():
        Data_4_CondMI = v

    Num_of_FBAs = Data_4_CondMI.keys()
    Count_Num_of_FBAs = Data_4_CondMI.values()
    print (Num_of_FBAs)
    print (Count_Num_of_FBAs)

    # -------Mutual Inforamtion (MI) calculation Stage II (input compounds respect to Biomass, Uptake and Secretion-------------
    # Pure Entropy
    FBAs = len(df[1][r].columns)
    pure_entropy = math.log(FBAs, 2)

    conditional_entropy = 0.0
    for l in range(0, len(Count_Num_of_FBAs)):
        temp = -1 * Count_Num_of_FBAs[l] * (Num_of_FBAs[l] / float(FBAs)) * Num_of_FBAs[l] * (
            1.0 / float(Num_of_FBAs[l]) * (math.log(1.0 / float(Num_of_FBAs[l]), 2)))
        conditional_entropy += temp

    Mutual_Info_B_U_S = pure_entropy - conditional_entropy
    # print('Mutaul Info:', Mutual_Info_B_U_S)

    MI_dict_b_u_s.update({df[0][r]: Mutual_Info_B_U_S})

# Sorted MI_dict_biomass
MI_dict_b_u_s = sorted(MI_dict_b_u_s.items(), key=lambda x: (-len(x[0]), x[0]))
MI_dict_b_u_s = OrderedDict(MI_dict_b_u_s)

print( MI_dict_b_u_s)

x_groups = [[] for x in range(num_compounds)]
y_groups = [[] for x in range(num_compounds)]
names = [[] for x in range(num_compounds)]
Comp_Mapping = [[] for x in range(num_compounds)]

for key, val in MI_dict_b_u_s.iteritems():
    del_count = key.count(',')
    x_groups[del_count].append(key)
    y_groups[del_count].append(val)

#for x, y in zip(x_groups, y_groups):
    #data.append(go.Bar(x=x, y=y, name='test'))

compound_IDs = ['H2', 'Vitamin K', 'Hematin', 'Glucose', 'Acetate', 'Formate', 'B12']

pdata = []
for i in range(0, len(x_groups)):
    names[i] = str(i + 1) + ' Compound Combination'
    Comp_Mapping = str(i + 1) + '-' + compound_IDs[i]

    record = {}
    record["x"] = []
    for e in x_groups[i]:
        record["x"].append("c" + e)
    record["y"] = y_groups[i]
    record["names"] = names[i]
    record["Comp_Mapping"] = Comp_Mapping
    pdata.append(record)

print pdata

#with open('pdata_BT_B_U_S.json', 'w') as outfile:
  #  json.dump(pdata, outfile)
#print json.dumps(pdata)

