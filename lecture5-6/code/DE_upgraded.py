# import modules
import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind, mannwhitneyu

print("The tool works only with two gene expression datasets! The tables should contain ONLY numeric gene expression values.\n Each column stands for a certain gene. The column order should be the SAME in the two tables. \n")
# read paths to the files
first_cell_type_expressions_path = str(input('Enter the path to the first dataset, including the file name (i.e. a table with gene expression values for cell type 1): '))
second_cell_type_expressions_path = str(input('Enter the path to the second dataset, including the file name (i.e. a table with gene expression values for cell type 2): '))
save_results_table = str(input('Enter the path to the output (without csv extension) '))

print("Select on of the following tests to perform: ztest, ci_test, ttest_ind, mannwhitneyu ")
test = input("Select the test (mandatory!): ")

print("Choose one of the methods for testing and adjustment of pvalues: bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky. \nPress Enter if no correction is needed.")
correction = input('Choose the correction method: ')

# read the data
first_table = pd.read_csv(first_cell_type_expressions_path, index_col=0)
second_table = pd.read_csv(second_cell_type_expressions_path, index_col=0)

if test == 'ztest': 
  method = ztest
elif test == 'ci_test':   
  method = 'ci_test'
elif test == 'ttest_ind':   
  method = ttest_ind
elif test == 'mannwhitneyu':   
  method = mannwhitneyu
# function that checks intervals intersection
def check_intervals_intersect(first_ci, second_ci):   
  are_intersect = first_ci[0] < second_ci[0] < first_ci[1] or second_ci[0] < first_ci[0] < second_ci[1]
  return are_intersect

# function that counts test results
def DE_tools(first_table, second_table, save_results_table, method=method, correction=None):

# lists for test results
    test_p_values, H0_rejected_list = [], []
    mean_diff = []


    for i in first_table.columns: 
    # ci test
        if method == 'ci_test':
        # for table 1
            ci_1 = st.t.interval(alpha=0.95, # 95% CI
                    df=len(first_table[i]) - 1, # df - 1
                    loc=np.mean(first_table[i]), # mean
                    scale=st.sem(first_table[i])) # se
            # for table 2
            ci_2 = st.t.interval(alpha=0.95, # 95% CI
                    df=len(second_table[i]) - 1, # df - 1
                    loc=np.mean(second_table[i]), # mean
                    scale=st.sem(second_table[i])) # se
            result = check_intervals_intersect(ci_1, ci_2)
            H0_rejected_list.append(result)

        else:
            #test
            result = method(first_table[i], second_table[i])
            H0_rejected = (result[1] < 0.05)
            H0_rejected_list.append(H0_rejected)
            test_p_values.append(result[1])

        # mean_diff
        mean_diff_each = np.mean(first_table[i]) - np.mean(second_table[i])
        mean_diff.append(mean_diff_each)

    # correction for multiple tests
    names = list(first_table.columns)
    if method != 'ci_test':
    #dictionary for df columns
        results = {
            "gene_name": names,
             test + "_p_values": test_p_values,
            "adjusted_p_values": [],
            "mean_diff": mean_diff
        }
    else:
        results = {
            "gene_name": names,
            "mean_diff": mean_diff,
            "H0_rejected": H0_rejected_list
        }
    if correction:
        adjusted_p_value = multipletests(test_p_values, method=correction)
        results["adjusted_p_values"] = adjusted_p_value[1]
        results["H0_rejected"] = adjusted_p_value[0]
    if correction == None:
        results["H0_rejected"] = H0_rejected_list



    # dataframe from dictionary
    results = pd.DataFrame(results)
    # save results to .csv
    results.to_csv(f"{save_results_table}.csv")


DE_tools(first_table, second_table, save_results_table, method=method, correction=correction)
print("Done!")
