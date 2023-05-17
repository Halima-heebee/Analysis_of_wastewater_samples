
import pandas as pd
import ast

df = pd.read_csv('/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/ontdemixed-ww-demixed.tsv', sep='\t')
df.rename(columns={'Unnamed: 0': 'SAMPLE'}, inplace=True)
df["facet"] = df["SAMPLE"].str.extract("^([^-]+)")[0]
df['group'] = df.SAMPLE.str[:7]

# Create a list of the lineage names we want to extract
lineages = ['Alpha', 'Beta', 'Delta', 'Omicron', 'Eta', 'Gamma', 'Lambda', 'Other']

# Loop over each row in the dataframe and extract the lineage abundances
for index, row in df.iterrows():
    # Initialize an empty list to hold the lineage abundances for this row
    lineage_abundances = []
    # Convert the summarized column from a string to a list of tuples
    summarized_list = ast.literal_eval(row['summarized'])
    # Loop over each lineage we want to extract
    for lineage in lineages:
        # Search for the lineage in the summarized list
        found = False
        for tup in summarized_list:
            if tup[0] == lineage:
                # Add the lineage abundance to the list for this row
                lineage_abundances.append(tup[1])
                found = True
                break
        # If the lineage wasn't found, add a zero to the list for this row
        if not found:
            lineage_abundances.append(0)
    # Add the lineage abundances to the dataframe as new columns
    for i, lineage in enumerate(lineages):
        df.at[index, lineage] = lineage_abundances[i]

lineages = ['Alpha', 'Beta', 'Delta', 'Omicron', 'Eta', 'Lambda', 'Gamma', 'Other']
new_rows = []

for index, row in df.iterrows():
    sample = row['SAMPLE']
    facet = row['facet']
    group = row['group']
    for lineage in lineages:
        new_row = {
            'SAMPLE': sample,
            'facet': facet,
            'group': group,
            'lineage': lineage,
            'abundance': row[lineage]
        }
        new_rows.append(new_row)

new_df = pd.DataFrame(new_rows)
new_df['group'] = new_df['group'].apply(lambda x: x if x.startswith('ill') else 'ONT')
new_df["group"] = new_df["group"].apply(lambda x: "ILL" if x.startswith("ill") else x)
new_df['facet'] = new_df['facet'].str.replace('ill.', '')

# Write the updated dataframe to a new file
new_df.to_csv('updated_file_all.csv', index=False)

