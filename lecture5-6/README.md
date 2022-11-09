# DE_tool_upgraded
Compares two datasets and finds differentially expressed genes
## Usage
The tool works with two gene expression datasets. The tables should contain ONLY numeric gene expression values.
Each column stands for a certain gene. The column order should be the SAME in the two tables.

## Parameters:
The tool takes the paths to the first and second datasets and to the output table with test results
The following tests are available: z-test, confidence interval test, t-test for independent variables and Mann-Whitney U test.
The user can also choose one of the methods to correct for multiple comparisons: 
bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky.
No correction is performed by default.
The output table provides gene names, test_p_values (except for CI test), adjusted_p_values (if any), mean difference and if H0 is rejected.
