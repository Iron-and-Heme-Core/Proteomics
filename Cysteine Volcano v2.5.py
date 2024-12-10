import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
import plotly.express as px
import matplotlib as mpl
import matplotlib.pyplot as plt
from adjustText import adjust_text
from itertools import permutations

# Step 0: Define the groups and tissue
tissue = 'SKM'                                                # Define tissue name
# tissue = 'PBMC'

groups = ['YOUNG-PRE', 'YOUNG-POST', 'OLD-PRE', 'OLD-POST']   # groups = ['YOUNG-PRE', 'YOUNG-POST', 'OLD-PRE', 'OLD-POST', 'YOUNG-24hr POST', 'OLD-24hr POST']  # All possible groups
range_slot = ''                                               # range_slot = ' 79-98'  range_slot = ' 51-78_99-100' # range_slot = ' 51-78_99-100' # range_slot = ' 29-50' # range_slot = ' 1-28'

directory = r'C:\Users\\'
directory_figures = r'C:\Users\\' + tissue + range_slot

### Begin MAIN LOOP for all pairwise comparisons between groups and generate volcano plots
for group_1, group_2 in permutations(groups, 2):              # Generate all pairwise permutations
    try:
        # Define file names and tags
        png_name = f'{tissue}{range_slot} - {group_2} vs {group_1} Volcano.png'
        data_tag = f'Volcano Data Oxidation p-vals {tissue}{range_slot}_{group_1}_{group_2}.csv'
        file_path = directory + data_tag

        # Step 1: Load the data
        try:
            data = pd.read_csv(file_path)
        except FileNotFoundError:
            print(f"File not found for {group_1} vs {group_2}: {data_tag}")
            continue

        #print(f"Processing: {group_1} vs {group_2}")

        # Step 2: Define columns and drop missing values
        percent_change = fr'{group_2}_minus_{group_1}'
        p_val = 'p_value'
        data[percent_change] = data[percent_change].replace(0, np.nan)
        data[p_val] = data[p_val].replace(0, np.nan)
        data = data.dropna()

        # Step 3: Define thresholds and calculate -log10(p)
        p_val_threshold = 0.05
        threshold_percent_change = 10
        x_axis_range = 100
        threshold_log_pval = -np.log10(p_val_threshold)
        data['-log10(P_Value)'] = -np.log10(data['p_value'])

        # Step 4: Perform FDR correction
        fdr_corrected = multipletests(data['p_value'], method='fdr_bh')
        data['FDR_P_Value'] = fdr_corrected[1]
        data['-log10(FDR_P_Value)'] = -np.log10(data['FDR_P_Value'])

        # Step 5: Define classification function
        def classify_point_p(row):
            if row[percent_change] < -threshold_percent_change and row['-log10(P_Value)'] >= threshold_log_pval:
                return 'Decreased Oxidation'
            elif row[percent_change] > threshold_percent_change and row['-log10(P_Value)'] >= threshold_log_pval:
                return 'Increased Oxidation'
            else:
                return 'Not Significant'

        # Step 6: Apply classification function and extract protein names
        volcano_data = data
        volcano_data['Color_matplotlib'] = volcano_data.apply(classify_point_p, axis=1)
        volcano_data['Protein_Name'] = volcano_data['Protein_Peptide_Concat'].str.split(' -').str[0]  # Extract protein name

        # Step 7: Generate the volcano plot
        fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
        colors = {'Increased Oxidation': 'red', 'Decreased Oxidation': 'blue', 'Not Significant': 'gray'}
        mpl.rcParams['font.family'] = 'Arial'

        for color in colors:
            subset = volcano_data[volcano_data['Color_matplotlib'] == color]
            ax.scatter(subset[percent_change], subset['-log10(P_Value)'],
                       c=colors[color], label=color, alpha=0.8, s=30, edgecolor='white', linewidth=0.5)

        # Add text labels to significant points
        texts = []
        for _, row in volcano_data.iterrows():
            if abs(row[percent_change]) >= threshold_percent_change and row['-log10(P_Value)'] >= threshold_log_pval:
                texts.append(ax.text(row[percent_change], row['-log10(P_Value)'], row['Protein_Name'], fontsize=6))

        # Dynamically adjust text positions to avoid overlap
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='grey', alpha=0.5, lw=0.5))

        # Finalize plot details
        ax.set_xlabel('Change in Oxidation (%)')
        ax.set_ylabel('-log (p)')
        ax.set_title(f'{tissue} - % Change in Cysteine Oxidation - {group_2} vs. {group_1}')
        ax.axhline(y=threshold_log_pval, color='gray', linestyle='--', linewidth=0.6)
        ax.axvline(x=-threshold_percent_change, color='gray', linestyle='--', linewidth=0.6)
        ax.axvline(x=threshold_percent_change, color='gray', linestyle='--', linewidth=0.6)
        ax.set_xlim([-x_axis_range, x_axis_range])
        plt.legend(loc='lower right', fontsize='small', fancybox=True)
        fig.tight_layout()

        # Save and show the plot
        fig.savefig(directory_figures + png_name, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved: {png_name}")

    except Exception as e:
        print(f"An error occurred for {group_1} vs {group_2}: {e}")


print("All pairwise comparisons completed")

