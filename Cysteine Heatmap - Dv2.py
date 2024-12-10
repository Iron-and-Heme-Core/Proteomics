import os
import glob
import re
import pandas as pd
import plotly.graph_objects as go
import plotly.figure_factory as ff
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist



# Step 0: Define variables
filter_percent = 0              #### UNUSED
filter_min_oxidation = 0       # Minimum oxidation percentage to keep (Removes anything with average average <= 0)
filter_nondetect_max = 75      # Maximum number of non-detects allowed

tissue = 'PBMC'
# tissue = 'SKM'

cmap = 'blues'  # Colormap for the heatmap


#desired_order = [f'{tissue}_Average_OLD-24hr POST', f'{tissue}_Average_OLD-POST', f'{tissue}_Average_OLD-PRE',
#                 f'{tissue}_Average_YOUNG-24hr POST', f'{tissue}_Average_YOUNG-POST', f'{tissue}_Average_YOUNG-PRE']  # Replace with actual names
desired_order = [f'OLD-POST_Average_{tissue}', f'OLD-PRE_Average_{tissue}',
                 f'YOUNG-POST_Average_{tissue}', f'YOUNG-PRE_Average_{tissue}']  # Replace with actual names


# Step 1: Define the directory where the CSV files are stored and create a list for all DataFrames
directory_path = r'C:\Users\\'
output_directory_path = r'C:\Users\\'

# Define group using metadata group column name
group_column = 'Age-Exercise'
metadata_path = r'C:\Users\sample_metadata_concat.csv'

file_pattern = os.path.join(directory_path, f'PROCESSED*{tissue}*.csv')  # Match all CSV files in the directory
file_paths = glob.glob(file_pattern)


# Step 2: Read the metadata file and build a dictionary of groups of DataFrames
metadata_df = pd.read_csv(metadata_path)
metadata_df['Number'] = metadata_df['Number'].astype(str)  # Ensure 'Number' is a string for matching

# Initialize a dictionary to hold groups of DataFrames
grouped_dfs = {group: [] for group in metadata_df[group_column].unique()}

# Read all the files into DataFrames and group them based on metadata
for file_path in file_paths:
    file_name = os.path.basename(file_path).replace('.csv', '')
    file_number = re.search(r'\d+', file_name).group()  # Extract the number from the file name
    group = metadata_df.loc[metadata_df['Number'] == file_number, group_column].values[0]
    df = pd.read_csv(file_path)
    grouped_dfs[group].append(df)

print(grouped_dfs.keys())  # Should print the unique groups

#######################################################################################################################

merged_df = pd.DataFrame()
# Step 3: Iterate over each file and merge the data
for file_path in file_paths:
    # Extract the DataFrame name from the file name
    df_name = os.path.basename(file_path).replace('.csv', '')
    df_name = re.search(r'\d+', df_name).group()
    df_name = df_name + " " + metadata_df.loc[metadata_df['Number'] == df_name, group_column].values[0]

    # Load the DataFrame
    df = pd.read_csv(file_path)
    temp_df = df[['Protein_Peptide_Concat', 'Percent Oxidation']].copy()     # Retain only the necessary columns

    # Rename the 'Percent Oxidation' column to include the DataFrame name
    temp_df.rename(columns={'Percent Oxidation': f'{df_name}'}, inplace=True)

    # Merge the temporary DataFrame into the master DataFrame
    if merged_df.empty:
        merged_df = temp_df
    else:
        merged_df = pd.merge(
            merged_df, temp_df,
            on='Protein_Peptide_Concat',
            how='outer'
        )

# Step 4: Save the final merged DataFrame
output_file = os.path.join(output_directory_path, 'final_merged_dataframe.csv')
merged_df.to_csv(output_file, index=False)

print(f"Final merged DataFrame saved to {output_file}")

# Step 5: Construct a group merged heatmap for each group from the full merged DataFrame
for group in metadata_df[group_column].unique():
    group_columns = [col for col in merged_df.columns if group in col]
    group_df = merged_df[['Protein_Peptide_Concat'] + group_columns].copy()

    average_column = f'{group}_Average'
    group_df[average_column] = group_df[group_columns].mean(axis=1)

    non_detect_column = f'{group}_NDs'
    group_df[non_detect_column] = group_df[group_columns].isna().sum(axis=1)
    group_df[non_detect_column] = group_df[non_detect_column] / len(group_columns) * 100
    group_df.to_csv(output_directory_path + f'{group}_merged_heatmap.csv', index=False)
    group_df = group_df[group_df[non_detect_column] < filter_nondetect_max]
    group_df.to_csv(output_directory_path + f'{group}_merged_heatmap_filtered.csv', index=False)

# Step 6: Construct a dataframe for the heatmap. This will have the average values for each group
heatmap_df2 = pd.DataFrame()
for group in metadata_df[group_column].unique():
    group_df = pd.read_csv(output_directory_path + f'{group}_merged_heatmap_filtered.csv')
    group_df = group_df[['Protein_Peptide_Concat', f'{group}_Average']]
    group_df = group_df.rename(columns={f'{group}_Average': f'{group}_Average_{tissue}'})
    if heatmap_df2.empty:
        heatmap_df2 = group_df
    else:
        heatmap_df2 = pd.merge(heatmap_df2, group_df, on='Protein_Peptide_Concat', how='outer')


# Step 7: Force removal of 24hr columns from the heatmap dataframe
# THIS IS A TEMPORARY FIX, REMOVE THIS ONCE THE PROCESSED DATA DIRECTORY IS CLEANED
heatmap_df2 = heatmap_df2[heatmap_df2.columns.drop(list(heatmap_df2.filter(regex='24hr')))]


# Step 8: Remove rows where where values are below the filter_min_oxidation threshold
heatmap_df2 = heatmap_df2[heatmap_df2.iloc[:, 1:].min(axis=1) > filter_min_oxidation]



# Step 9: Initialize the heatmap data DataFrame
heatmap_df2 = heatmap_df2.fillna(0)
heatmap_peptides = heatmap_df2['Protein_Peptide_Concat']
heatmap_peptides = heatmap_peptides.reset_index(drop=False)
heatmap_data = heatmap_df2.set_index('Protein_Peptide_Concat').T
heatmap_data.to_csv(output_directory_path + f'{tissue}_heatmap_data.csv')

# Step 10: Perform hierarchical clustering on the heatmap data
row_linkage = linkage(pdist(heatmap_data.values, metric='euclidean'), method='ward')
row_dendrogram = ff.create_dendrogram(heatmap_data.values, orientation='right', linkagefun=lambda x: row_linkage)
row_dendro_indices = list(map(int, row_dendrogram['layout']['yaxis']['ticktext']))  # Numeric indices for rows

col_linkage = linkage(pdist(heatmap_data.T.values, metric='euclidean'), method='ward')
col_dendrogram = ff.create_dendrogram(heatmap_data.T.values, orientation='bottom', linkagefun=lambda x: col_linkage)
col_dendro_indices = list(map(int, col_dendrogram['layout']['xaxis']['ticktext']))  # Numeric indices for columns

# Reorder heatmap data based on clustering results
row_peptide_order = heatmap_data.index[row_dendro_indices]  # Reordered row peptide names
col_peptide_order = heatmap_data.columns[col_dendro_indices]  # Reordered column peptide names
reordered_heatmap_data = heatmap_data.loc[desired_order, col_peptide_order] # Reordered heatmap data


# Step 11: Create and save the heatmap
fig = go.Figure(data=go.Heatmap(
    z=reordered_heatmap_data.values,
    x=reordered_heatmap_data.columns,
    y=reordered_heatmap_data.index,
    colorscale=cmap,
    zauto=False
))
fig.update_layout(
    title=dict(text=f"Cysteine Oxidation - {tissue}", font=dict(size=24)),
    xaxis_title="Peptides",
    xaxis=dict(showticklabels=False),
    font=dict(family="Arial", size=11)
)
fig.update_traces(connectgaps=False)
fig.update_traces(hovertemplate='Peptide: %{x}<br>Sample: %{y}<br>Oxidation: %{z:.2f}%')
fig.update_traces(xgap=.1, ygap=.1)
fig.write_html(output_directory_path + f'{tissue}_Cysteine_Ox_Heatmap.html', auto_open=True)
fig.write_image(output_directory_path + f"{tissue}_heatmap.png", width=2500, height=1500, scale=3)

print(f"Heatmap saved to {output_directory_path + f'{tissue}_Cysteine_Ox_Heatmap.html'}")
print(f"Image saved to {output_directory_path + f'{tissue}_heatmap.png'}")

