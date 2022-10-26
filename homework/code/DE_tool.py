# import modules
import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest

# read paths to the files
first_cell_type_expressions_path = str(input('Введите полный путь до таблицы с экспрессиями генов для одного клеточного типа'))
second_cell_type_expressions_path = str(input('Введите полный путь до таблицы с экспрессиями генов для второго клеточного типа. Названия и нормера колонок должны быть как и в первой таблице!'))
save_results_table = str(input('Введите название таблицы с результатами'))

# read the data
first_table = pd.read_csv(first_cell_type_expressions_path, index_col=0)
second_table = pd.read_csv(second_cell_type_expressions_path, index_col=0)

# lists for test results
ci_test_results = []
z_test_results = []
z_test_p_values = []
mean_diff = []

# function that checks intervals intersection
def check_intervals_intersect(first_ci, second_ci):   
  are_intersect = first_ci[0] < second_ci[0] < first_ci[1] or second_ci[0] < first_ci[0] < second_ci[1]
  return are_intersect

# function that counts test results
def DE_tools(first_table, second_table, ci_test_results = ci_test_results, z_test_results = z_test_results,
             z_test_p_values = z_test_p_values, mean_diff = mean_diff):

    for i in first_table.columns: 
    # ci test
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
        ci_test_result = check_intervals_intersect(ci_1, ci_2)
        ci_test_results.append(ci_test_result)

    # ztest
        z_test_result = ztest(first_table[i], second_table[i])[1] < 0.05
        z_test_results.append(z_test_result)
        z_test_p_value = ztest(first_table[i], second_table[i])[1]
        z_test_p_values.append(z_test_p_value)
  
   # mean_diff
        mean_diff_each = np.mean(first_table[i]) - np.mean(second_table[i])
        mean_diff.append(mean_diff_each)

DE_tools(first_table, second_table)
names = list(first_table.columns)
# Созданим словарь {'название колонки': список_значений}
results = {
    "gene_name": names,
    "ci_test_results": ci_test_results,
    "z_test_results": z_test_results,
    "z_test_p_values": z_test_p_values,
    "mean_diff": mean_diff
}

# Из словаря делаем датафрейм
results = pd.DataFrame(results)
# Сохраним таблицу в .csv файл
results.to_csv(f"{save_results_table}.csv")
