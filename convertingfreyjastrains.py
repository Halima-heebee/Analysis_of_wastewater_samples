
import pandas as pd
import ast

import pandas as pd
import ast
import glob


df = pd.read_csv('/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/ontdemixed-ww-demixed.tsv', sep='\t')
df.rename(columns={'Unnamed: 0': 'SAMPLE'}, inplace=True)
df["facet"] = df["SAMPLE"].str.extract("^([^-]+)")[0]
df['group'] = df.SAMPLE.str[:7]
df["lineages"] = df["lineages"].apply(lambda x:str(x))
df["abundances"] = df["abundances"].apply(lambda x:str(x))
# Create a list of the lineage names we want to extract
#lineages = []

#for index, row in df.iterrows():
#	lineage_labels = row['lineages'].split(' ')
#	lineages.extend(lineage_labels)

df['plate_well'] = df['SAMPLE'].str[4:7]
df['group'] = df['SAMPLE'].str[0:3]
df['strain'] = df['lineages'].str.split()
df['abundance'] = df['abundances'].str.split()
new_df = df.explode(['strain', 'abundance'])
new_df = new_df.reset_index(drop=True)
new_df = new_df[["SAMPLE", "plate_well", "group", "strain", "abundance"]]
#new_df.columns[["SAMPLE", "plate_well", "group", "strain", "abundance"]]
print(new_df)
new_df.to_csv('updated_strains.csv', index=False)

