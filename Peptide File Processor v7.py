from itertools import product
import pandas as pd
import re
import os
import numpy as np
from scipy.stats import ttest_ind
import warnings


# Step 1: Define directories and parameters
tissue_name = 'PBMC'
range_slot = ' 29-50'  # range_slot = ' 1-28'   range_slot = ' 51-78_99-100'    range_slot = ' 79-98'   range_slot = ' 29-50'
age_1, age_2 = 'YOUNG', 'OLD'
exercise_1, exercise_2, exercise_3 = 'PRE', 'POST', '24hr POST'
# ages, exercises = [age_1, age_2], [exercise_1, exercise_2]     # exercises = [exercise_1, exercise_2, exercise_3]
ages, exercises = [age_1, age_2], [exercise_1, exercise_2, exercise_3]

# root_dir = r'C:\Users\u0716965.AD\PycharmProjects\python-lectures-master\Protein Scripts\\'
# dataset_directory = root_dir + 'Raw Data\\' + tissue_name + range_slot
metadata_directory = r'C:\Users\u0716965.AD\PycharmProjects\python-lectures-master\Protein Scripts\Volcano Plot - FINAL\sample_metadata_concat.csv'
dataset_directory = r'C:\Users\u0716965.AD\PycharmProjects\python-lectures-master\Protein Scripts\Raw Data\\' + tissue_name + range_slot
output_directory = r'C:\Users\u0716965.AD\PycharmProjects\python-lectures-master\Protein Scripts\Volcano Plot - FINAL\v7 test\process test\\' + tissue_name + range_slot
merged_output_directory = r'C:\Users\u0716965.AD\PycharmProjects\python-lectures-master\Protein Scripts\Volcano Plot - FINAL\v7 test\process test\\'


# Step 2: Load metadata
metadata_df = pd.read_csv(metadata_directory)
group_column = 'Age-Exercise'


#Step 3: Define helper functions
#
# Function to split peptides with cysteine positions
def split_peptides_with_cysteine(row):
    """
    Splits a peptide row into multiple rows based on cysteine positions.

    For each reduced or oxidized cysteine site, a new row is created. The
    cysteine type (Reduced or Oxidized) is recorded, and the cysteine position
    is appended to the peptide sequence.

    Args:
        row (pd.Series): A row containing 'Reduced Cysteine' and 'Oxidized Cysteine' columns.

    Returns:
        list: A list of new rows with updated 'Peptide Sequence' and 'Cysteine Type' information.
    """
    new_rows = []
    if row['Reduced Cysteine']:
        for site in row['Reduced Cysteine'].split(', '):
            new_row = row.copy()
            new_row['Peptide Sequence'] = f"{row['Peptide Sequence']} - {site}C"
            new_row['Cysteine Type'] = 'Reduced'
            new_rows.append(new_row)
    if row['Oxidized Cysteine']:
        for site in row['Oxidized Cysteine'].split(', '):
            new_row = row.copy()
            new_row['Peptide Sequence'] = f"{row['Peptide Sequence']} - {site}C"
            new_row['Cysteine Type'] = 'Oxidized'
            new_rows.append(new_row)
    return new_rows


def clean_modifications(seq):
    """
    Removes specific modification annotations from a peptide sequence.

    Args:
        seq (str): A peptide sequence with modification annotations (e.g., n[42.x]).

    Returns:
        str: A cleaned peptide sequence without modifications.
    """
    seq = re.sub(r'n\[42\.\d+\]', '', seq)
    seq = re.sub(r'M\[15\.\d+\]', 'M', seq)
    return seq


def normalize_leading_trailing(seq):
    """
    Normalizes a peptide sequence by removing excessive leading or trailing amino acids.

    Leading K/R residues are stripped, and trailing patterns (KK, RR, KR, RK)
    are normalized to a single amino acid.

    Args:
        seq (str): A peptide sequence.

    Returns:
        str: The normalized sequence.
    """
    seq = re.sub(r'^(K|R)+', '', seq)
    seq = re.sub(r'(KK|RR|KR|RK)$', lambda m: m.group(0)[0], seq)
    return seq


def calculate_percent_oxidation(group):
    """
    Calculates the percent oxidation of cysteine residues.

    Args:
        group (pd.DataFrame): A DataFrame group with 'Intensity' and 'Cysteine Type' columns.

    Returns:
        pd.Series: A series with total intensity, oxidized intensity, and percent oxidation.
    """
    total_intensity = group['Intensity'].sum()
    oxidized_intensity = group[group['Cysteine Type'] == 'Oxidized']['Intensity'].sum()
    percent_oxidation = (oxidized_intensity / total_intensity) * 100 if total_intensity > 0 else 0
    return pd.Series({
        'Total Intensity': total_intensity,
        'Oxidized Intensity': oxidized_intensity,
        'Percent Oxidation': round(percent_oxidation, 2)
    })


def clean_sample_name(file_name):
    """
    Extracts and formats a sample name using metadata.

    Args:
        file_name (str): The filename containing a slot number.

    Returns:
        str: A formatted sample name or 'Unknown' if not found.
    """
    match = re.search(r'Slot (\d+)', file_name)
    if match:
        number = match.group(1)
        group = metadata_dict.get(number, 'Unknown')
        return f'{number}_{group}'
    return 'Unknown'


def replace_cysteines_beta(row):
    """
    Replaces cysteine modifications in a peptide sequence with placeholders.

    Args:
        row (pd.Series): A row with 'Modified Sequence'.

    Returns:
        pd.Series: The updated row with 'Modified Sequence Placeholder'.
    """
    modified_seq = re.sub(r'C\[57\.\d+\]', 'J', row['Modified Sequence'])
    modified_seq = re.sub(r'C\[61\.\d+\]', 'X', modified_seq)
    row['Modified Sequence Placeholder'] = modified_seq
    return row


def count_cysteines_beta(row):
    """
    Counts the positions of reduced (J) and oxidized (X) cysteines in a peptide sequence.

    Args:
        row (pd.Series): A row with 'Modified Sequence Placeholder'.

    Returns:
        pd.Series: The updated row with 'Reduced Cysteine' and 'Oxidized Cysteine' positions.
    """
    row['Reduced Cysteine'] = ', '.join(
        [str(i + 1) for i, letter in enumerate(row['Modified Sequence Placeholder']) if letter == 'J'])
    row['Oxidized Cysteine'] = ', '.join(
        [str(i + 1) for i, letter in enumerate(row['Modified Sequence Placeholder']) if letter == 'X'])
    return row


def split_peptide(peptide):
    """
    Splits a peptide into fragments based on miscleavage sites (K or R not followed by P).

    Args:
        peptide (str): A peptide sequence.

    Returns:
        list: A list of peptide fragments.
    """
    pattern = r"(?<=[KR])(?!P)(?=[A-Z])"
    split_points = [0] + [m.end() for m in re.finditer(pattern, peptide)] + [len(peptide)]
    return [peptide[split_points[i]:split_points[i + 1]] for i in range(len(split_points) - 1)]


def process_daughter_peptides(df, peptide_column):
    """
    Identifies miscleavages and splits rows into daughter peptides.

    Args:
        df (pd.DataFrame): The input DataFrame.
        peptide_column (str): The column containing peptide sequences.

    Returns:
        pd.DataFrame: A DataFrame with additional rows for each daughter peptide.
    """
    split_rows = []
    for _, row in df.iterrows():
        daughters = split_peptide(row[peptide_column])

        if len(daughters) > 1:  # If miscleavages detected, process daughter peptides
            daughter_peptides = ', '.join(daughters)  # Join daughter peptides with a comma
            row['Daughter Peptides'] = daughter_peptides  # Record daughter peptides in the row
            for daughter in daughters:
                new_row = row.copy()
                new_row[peptide_column] = daughter
                split_rows.append(new_row)

        else:
            # If no daughter peptides, keep the row as is
            row['Daughter Peptides'] = ''
            split_rows.append(row)

    return pd.DataFrame(split_rows)


def process_dataframe_beta(df, filename):
    """
        Processes a DataFrame of peptide data to clean, modify, split, and calculate oxidation.

        Args:
            df (pd.DataFrame): Input DataFrame containing peptide data with columns such as
                'Peptide Sequence', 'Modified Sequence', 'Assigned Modifications', and 'Intensity'.
            filename (str): The filename of the input data (for logging/debugging).

        Returns:
            pd.DataFrame: A cleaned and processed DataFrame with key columns such as
                'Peptide Sequence', 'Protein ID', 'Percent Oxidation', and concatenated IDs.
    """
    # Step 1: Drop rows with all NaN values and select key columns
    df = df.dropna(how='all')
    intensity_col = [col for col in df.columns if "Intensity" in col][0]

    # Step 2: Select key columns and filter out rows with no cysteines or zero intensity
    df = df[['Peptide Sequence', 'Modified Sequence', 'Assigned Modifications', 'Protein ID', intensity_col]]
    df = df.rename(columns={intensity_col: 'Intensity'})
    df = df[df['Peptide Sequence'].str.contains('C')]
    df = df[df['Assigned Modifications'].str.contains('C', na=False)]
    df = df[df['Intensity'] != 0]

    # Step 3: Clean and normalize peptide sequences
    df['Modified Sequence'] = df['Modified Sequence'].apply(clean_modifications)  # Clean modifications
    df['Peptide Sequence'] = df['Peptide Sequence'].apply(normalize_leading_trailing)  # Normalize leading/trailing
    df['Modified Sequence'] = df['Modified Sequence'].apply(
        normalize_leading_trailing)  # Normalize leading/trailing
    df['pre_concatenated_peptide_sequence'] = df[
        'Peptide Sequence']  # Save original (cleaned and normalized) Peptide Sequence before appending cysteine information

    # Step 4: Replace cysteines with placeholders and count cysteine positions
    df = df.apply(replace_cysteines_beta, axis=1)

    # Step 5: Drop unnecessary columns
    df = df.drop(columns=['Modified Sequence', 'Assigned Modifications'])

    # Step 6: Process daughter peptides
    df = process_daughter_peptides(df, 'Modified Sequence Placeholder')

    # Step 7: Drop columns whose Modified Sequence Placeholder is less than 6 characters
    df = df.drop(columns=['Daughter Peptides'])

    # Step 8: Drop rows with no J or X in the Modified Sequence Placeholder
    df = df[df['Modified Sequence Placeholder'].str.len() >= 6]
    df = df.apply(count_cysteines_beta, axis=1)

    # Step 9: Replace Peptide Sequence with Modified Sequence Placeholder and replace J and X with C in the Peptide Sequence
    df = df[(df['Modified Sequence Placeholder'].str.contains('J')) |
            (df['Modified Sequence Placeholder'].str.contains('X'))]

    # Step 10: Replace Peptide Sequence with Modified Sequence Placeholder and replace J and X with C in the Peptide Sequence
    df['Peptide Sequence'] = df['Modified Sequence Placeholder'].str.replace('J', 'C').str.replace('X', 'C')
    df['pre_concatenated_peptide_sequence'] = df['Modified Sequence Placeholder'].str.replace('J', 'C').str.replace(
        'X', 'C')

    # Step 11: Split peptides with cysteine positions
    df_split = pd.DataFrame(
        [new_row for index, row in df.iterrows() for new_row in split_peptides_with_cysteine(row)])

    # Step 12: Group by peptide sequence, protein ID, and cysteine type and sum intensities
    df_split = df_split.groupby(['Peptide Sequence', 'Protein ID', 'Cysteine Type'], as_index=False).agg(
        {'Intensity': 'sum'})

    # Step 13: Calculate percent oxidation for each peptide sequence
    df_with_percent_oxidation = df_split.groupby('Peptide Sequence').apply(
        calculate_percent_oxidation).reset_index()
    df_final = pd.merge(df_split, df_with_percent_oxidation[['Peptide Sequence', 'Percent Oxidation']],
                        on='Peptide Sequence', how='left')

    # Step 14: Concatenate protein ID and peptide sequence
    df_final['pre_concatenated_peptide_sequence'] = df_final['Peptide Sequence'].str.split(' - ').str[0]
    df_final['Protein ID - Peptide'] = df_final['Protein ID'] + ' - ' + df_final[
        'pre_concatenated_peptide_sequence']
    df_final['Protein_Peptide_Concat'] = df_final['Protein ID'] + ' - ' + df_final['Peptide Sequence']

    # Step 15: Select key columns and drop duplicates
    df_reduced = df_final[['pre_concatenated_peptide_sequence', 'Peptide Sequence', 'Protein ID - Peptide',
                           'Protein ID', 'Protein_Peptide_Concat', 'Percent Oxidation']]
    df_reduced = df_reduced.sort_values(by='Protein_Peptide_Concat').drop_duplicates(
        subset=['Peptide Sequence', 'Protein ID'])

    return df_reduced


# Step 3: Check if all files are processed
# Flag to track whether the initial processed check has been done
# Filter phrase
start_phrase = 'Slot'
csv_dict = {}
summed_csv_dict = {}

initial_check_done = False
all_files_processed = False


# Check if all files are already processed
if not initial_check_done:
    all_files_processed = all(
        os.path.exists(os.path.join(output_directory, f"PROCESSED_{file}"))
        for file in os.listdir(dataset_directory)
        if file.endswith('.csv') and file.startswith(start_phrase)
    )
    initial_check_done = True

if all_files_processed:
    print("All files are already processed. Skipping further checks.")
else:
    # Process unprocessed files
    for file in os.listdir(dataset_directory):
        if file.endswith('.csv') and file.startswith(start_phrase):
            processed_file = f"PROCESSED_{file}"
            processed_output_path = os.path.join(output_directory, processed_file)

            if os.path.exists(processed_output_path):
                print(f"Skipping {file} as processed file already exists.")
                continue

            file_path = os.path.join(dataset_directory, file)
            df = pd.read_csv(file_path)

            processed_df = process_dataframe_beta(df, file)
            processed_df.to_csv(processed_output_path, index=False)
            print(f"Processed and saved {processed_file}.")


# Step 5: Generate group comparisons
groupcomps = []
for (age_1, exercise_1), (age_2, exercise_2) in product(product(ages, exercises), repeat=2):
    group_1 = f"{age_1}-{exercise_1}"
    group_2 = f"{age_2}-{exercise_2}"
    if group_1 != group_2:
        groupcomps.append(f"{group_1}_{group_2}")

# Print all valid group combinations
for groupcomp in groupcomps:
    print(groupcomp)


# Begin loop
for groupcomp in groupcomps:
    print(f'Processing group comparison: {groupcomp}')             # Print the group comparison
    groupcomp1, groupcomp2 = groupcomp.split('_')                # Split the group comparison into two groups
    file_tag = ' ' + tissue_name + range_slot + '_' + groupcomp1 + '_' + groupcomp2     # Create a file tag for the output files

    # Define output file names
    output_file_merge = 'All Oxidation Ratio' + file_tag + '.csv'
    output_file_with_averages = 'All Oxidation Ratio with Averages' + file_tag + '.csv'
    output_file_ox_pvals = 'All Oxidation Ratio with p-vals' + file_tag + '.csv'
    output_file_volcano_ox_pvals = 'Volcano Data Oxidation p-vals' + file_tag + '.csv'


    # Update metadata filtering to dynamically use groupcomp1 and groupcomp2
    metadata_df_filter = metadata_df[
        metadata_df[group_column].str.contains(groupcomp1) | metadata_df[group_column].str.contains(groupcomp2)
        ]

    # Update metadata dict and group names dynamically
    metadata_dict = dict(zip(metadata_df_filter['Number'].astype(str), metadata_df_filter[group_column]))
    group_names = [groupcomp1, groupcomp2]
    metadata_df_filter = metadata_df_filter[metadata_df_filter[group_column].isin(group_names)]

    csv_dict = {}  # Initialize a dictionary to hold all processed DataFrames

    # Load all processed files and merge them
    for file in os.listdir(output_directory):
        if file.endswith('.csv'):
            processed_file_path = os.path.join(output_directory, file)  # Get the full file path
            df = pd.read_csv(processed_file_path)                       # Read the file into a DataFrame
            if file.startswith("PROCESSED_") and 'Percent Oxidation' in df.columns:
                csv_dict[file] = df                                     # Add the DataFrame to the dictionary

    # Merge all processed files for percent oxidation
    merged_df = None
    for file_name, df in csv_dict.items():
        # Check if the DataFrame contains the necessary columns
        if 'Peptide Sequence' in df.columns and 'Protein ID' in df.columns:
            oxidation_column = [col for col in df.columns if 'Percent Oxidation' in col]
            #
            if oxidation_column:
                oxidation_column = oxidation_column[0]
                sample_name = clean_sample_name(file_name)
                # Select only the necessary columns
                df = df[['pre_concatenated_peptide_sequence', 'Peptide Sequence', 'Protein ID',
                         'Protein ID - Peptide', 'Protein_Peptide_Concat', oxidation_column]]
                # Rename the oxidation column to the sample name
                df = df.rename(columns={oxidation_column: f'{sample_name}'})

                # Merge the DataFrame into the master DataFrame
                if merged_df is None:
                    merged_df = df
                else:
                    merged_df = pd.merge(merged_df, df, on=['pre_concatenated_peptide_sequence', 'Peptide Sequence',
                                                            'Protein ID', 'Protein ID - Peptide',
                                                            'Protein_Peptide_Concat'], how='outer')
                # remove columns with the word "Unknown" in them
                merged_df = merged_df.loc[:, ~merged_df.columns.str.contains('Unknown')]

    # Final adjustments and save merged oxidation DataFrame
    test_df = merged_df
    if merged_df is not None:
        # Iterate over groups dynamically
        group_columns = {}
        for group in group_names:
            group_columns[group] = [col for col in merged_df.columns if f'_{group}' in col]

        # Sort the DataFrame by 'Peptide Sequence'
        merged_df.sort_values(by='Peptide Sequence', inplace=True)
        merged_output_file = os.path.join(merged_output_directory, output_file_merge)
        merged_df.to_csv(merged_output_file, index=False)

        # Dynamically calculate group averages
        for group, columns in group_columns.items():
            if columns:
                merged_df[f'{group}_Average'] = merged_df[columns].mean(axis=1)    # Calculate the average for each group
                blank_count_average = merged_df[f'{group}_Average'].isna().sum()   # Count the number of blanks
                print(f'Group {group} Blank Averages: {blank_count_average}')      # Print the number of blanks

        # Save the final DataFrame with averages
        merged_output_file_with_averages = os.path.join(merged_output_directory, output_file_with_averages)
        merged_df.to_csv(merged_output_file_with_averages, index=False)
        print(f'Merged oxidation DataFrame with averages saved as {output_file_with_averages}.')

        # Create filter_diff_merged_df
        if len(group_names) == 2:  # Ensure only two groups exist
            group1, group2 = group_names    # Dynamically assign the two groups
            filter_diff_merged_df = merged_df[[
                'Protein ID - Peptide',
                'Protein_Peptide_Concat',
                f'{group1}_Average',
                f'{group2}_Average'
            ]].copy()

            # Calculate differences
            filter_diff_merged_df[f'{group1}_minus_{group2}'] = (
                filter_diff_merged_df[f'{group1}_Average'] - filter_diff_merged_df[f'{group2}_Average']
            )
            # Calculate the negative difference, negative difference is the same as Group 1 vs Group 2
            filter_diff_merged_df[f'{group2}_minus_{group1}'] = -filter_diff_merged_df[f'{group1}_minus_{group2}']


            # Save the filtered DataFrame
            filter_diff_output_file = os.path.join(merged_output_directory, f'Filtered Difference {output_file_with_averages}')
            filter_diff_merged_df.to_csv(filter_diff_output_file, index=False)
            print(f'Filtered differences DataFrame saved as Filtered Difference {output_file_with_averages}')

    print(f"Group Names: {group_names}, Length: {len(group_names)}")
    print(f"Number of groups with data: {len([cols for cols in group_columns.values() if cols])}")


    merged_df = test_df
    # Perform t-tests and save results
    if merged_df is not None:
        # Iterate over groups dynamically
        group_columns = {}
        for group in group_names:
            group_columns[group] = [col for col in merged_df.columns if f'_{group}' in col]
        # Ensure exactly two groups exist for the t-test
        if len(group_columns) == 2:
            group1, group2 = group_names  # Dynamically assign the two groups
            print(f"Performing t-tests between {group1} and {group2}.")

            # Add a column for p-values
            merged_df['p_value'] = np.nan

            # Iterate over each row in the DataFrame to perform t-tests
            for index, row in merged_df.iterrows():
                # Extract values for the first group, replacing missing data with NaN and converting to float
                group1_values = row[group_columns[group1]].replace('-', np.nan).dropna().astype(float)
                # Extract values for the second group in the same manner
                group2_values = row[group_columns[group2]].replace('-', np.nan).dropna().astype(float)

                # Check if there are sufficient data points in both groups
                if len(group1_values) > 1 and len(group2_values) > 1:
                    try:
                        # Suppress warnings related to t-test calculations (e.g., divide by zero warnings)
                        with warnings.catch_warnings():
                            warnings.filterwarnings('ignore', category=RuntimeWarning)
                            # Perform the t-test assuming equal variance (change to `equal_var=False` if unequal variances expected)
                            _, p_value = ttest_ind(group1_values, group2_values, equal_var=True)

                            # Store the resulting p-value in the DataFrame
                            merged_df.at[index, 'p_value'] = p_value
                    except Exception as e:  # Handle any errors during t-test calculations
                        print(f"Error performing t-test for row {index}: {e}")
                        merged_df.at[index, 'p_value'] = np.nan  # Assign NaN if an error occurs
                else:
                    # Assign NaN if there are insufficient data points for the t-test
                    merged_df.at[index, 'p_value'] = np.nan

            # Save the DataFrame with p-values
            oxidation_with_pvals_file = os.path.join(merged_output_directory, output_file_ox_pvals)
            merged_df.to_csv(oxidation_with_pvals_file, index=False)
            print(f'Merged oxidation DataFrame with p-values saved as {oxidation_with_pvals_file}.')
        else:
            print(f"Skipping t-tests: Exactly two groups are required, but {len(group_columns)} were found.")


    print(f"Group Names: {group_names}, Length: {len(group_names)}")
    print(f"Number of groups with data: {len([cols for cols in group_columns.values() if cols])}")

    ox_change_df = pd.read_csv(filter_diff_output_file)
    columns_to_keep_volc = ['Protein_Peptide_Concat', f'{group1}_minus_{group2}', f'{group2}_minus_{group1}', 'p_value']
    volc_df_ox_pvalues = pd.read_csv(oxidation_with_pvals_file)

    # Calculate differences
    volc_df_ox_pvalues[f'{group1}_minus_{group2}'] = (volc_df_ox_pvalues[f'{group1}_Average'] - volc_df_ox_pvalues[f'{group2}_Average'])
    volc_df_ox_pvalues[f'{group2}_minus_{group1}'] = -volc_df_ox_pvalues[f'{group1}_minus_{group2}']
    volc_df_ox_pvalues = volc_df_ox_pvalues[columns_to_keep_volc]

    # Calculate count of peptides remaining after each step
    volc_rows_1 = len(volc_df_ox_pvalues)
    volc_df_ox_pvalues = volc_df_ox_pvalues.dropna(subset=['p_value'])
    volc_rows_2 = len(volc_df_ox_pvalues)
    volc_df_ox_pvalues = volc_df_ox_pvalues.dropna(subset=[f'{group1}_minus_{group2}'])
    volc_rows_3 = len(volc_df_ox_pvalues)

    # Print the remaining peptide counts.
    print(f'Peptides in Raw Volcano Data:  {volc_rows_1}')
    print(f'Peptides in p-value Filtered Volcano Data:  {volc_rows_2}')
    print(f'Peptides in p-value and Oxidation Change Filtered Volcano Data:  {volc_rows_3}')

    # Save the final DataFrame
    volc_file_ox_pvals = os.path.join(merged_output_directory, output_file_volcano_ox_pvals)
    volc_df_ox_pvalues.to_csv(volc_file_ox_pvals, index=False)

# End loop that began
print('Data Processing Complete')
