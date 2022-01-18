import pandas as pd
import statistics as stats
import scikit_posthocs as sc_post
import re


class ImportedData:

    def __init__(self, df):
        self.main_dataframe = df
        self.DKOAng = df.loc[df['label'] == 'DKO_Ang'].iloc[:, 2:] # Removes first two cols
        self.DKOSham = df.loc[df['label'] == 'DKO_Sham'].iloc[:, 2:]
        self.TKOSham = df.loc[df['label'] == 'TKO_Sham'].iloc[:, 2:]
        self.WTAng = df.loc[df['label'] == 'WT_Ang'].iloc[:, 2:]
        self.WTSham = df.loc[df['label'] == 'WT_Sham'].iloc[:, 2:]
        self.results = pd.DataFrame()

    def calculate_values(self):
        # Calculate FC and q_values for every metabolite
        # Comparisons are:
        #     DKO_Sham vs WT_Sham
        #     DKO_Sham vs TKO_Sham
        #     TKO_Sham vs WT_Sham

        comparisons = [[self.DKOSham, self.WTSham], [self.DKOSham, self.TKOSham], [self.TKOSham, self.WTSham]]
        comp_strings = ['DKO_Sham vs WT_Sham', 'DKO_Sham vs TKO_Sham', 'TKO_Sham vs WT_Sham']

        metabolite_list = []
        comparison_list = []
        fc_list = []
        p_value_list = []
        q_value_list = []

        i = 0
        for comparison in comparisons:
            group1 = comparison[0]
            group2 = comparison[1]
            comparison_string = comp_strings[i]
            i += 1

            for colname in group1:
                # Calculate everything first
                group1_avg = stats.mean(group1[colname])
                group2_avg = stats.mean(group2[colname])
                if not group2_avg == 0:
                    fold_change = group1_avg / group2_avg
                else:
                    fold_change = 'Undefined'
                p_val = sc_post.posthoc_dunn([group1[colname], group2[colname]]).iloc[1][1]
                q_val = sc_post.posthoc_dunn([group1[colname], group2[colname]], p_adjust='fdr_bh').iloc[1][1]

                # Append values to lists
                metabolite_list.append(colname)
                comparison_list.append(comparison_string)
                fc_list.append(fold_change)
                p_value_list.append(p_val)
                q_value_list.append(q_val)

        # Check for duplicate values and correct the names
        i = 0
        j = 2
        first_match = True
        for colname in metabolite_list:
            if colname[-1].isdigit():
                base_name = re.match(r'(.*)(?=\.\d+)', colname)
                if base_name is not None:
                    base_name = base_name.group()
                    if first_match:
                        metabolite_list[i-1] = base_name + '_1'
                        first_match = False
                    metabolite_list[i] = base_name + '_' + str(j)
                    j += 1
                else:
                    first_match = True
                    j = 2
            else:
                first_match = True
                j = 2
            i += 1

        # Modify self.results with the final dataframe
        self.results = pd.DataFrame(list(zip(metabolite_list, comparison_list, fc_list, p_value_list, q_value_list)),
                                    columns=['Metabolite', 'Comparison', 'FC', 'p-value', 'q-value-BH'])
        return 0

    def save_results(self, filepath):
        self.results.to_csv(filepath, encoding='utf-8', index=False)
        return 0