R and C++ codes for all the analyses from the manuscript
==========================================================

##### C++ files######

1) joint_test.cpp

This is the c++ file that describes our main algorithm. Also included are the functions to run the variance component score test (that tests Genotype x Methylation x Tissue interaction effect)

2) jaguar.cpp

This is the c++ file that describes JAGUAR, which is available at https://cran.r-project.org/web/packages/JAGUAR/index.html

##### R files #######

1) simulations_joint_test.R

This R script is necessary to run Monte Carlo simulations on our joint model. 

2) simulations_variance_comp_test.R

This R script runs Monte Carlo simulations comparing statistical power obtained from the variance component score test (for the Genotype x Methylation x Tissue interaction effect) with our joint model and for good measure, likelihoojd ratio test (RLRsim).

3) simulations_JAGUAR_vs_JT.R

This R scripts runs Monte Carlo simulations comparing statistical power obtained from JAGUAR and our joint model both in the absence and presence of any effects due to methylation data.

4) preprocess_gibbs.R

This R scripts preprocesses Gibbs et al adult human brain data from four brain regions, originally published at http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000952
This script essentially creates multiple sub-directories and deposits even partitions of gene expression and methylation data.

5) run_joint_test.R

This R script runs our model on the Gibbs et al data fragments.

6) run_separate_tbt.R

This R script runs a separate region-by-region analysis of gene expression and methylation data.

7) run_TBT.R

This R script runs a region-by-region analysis using a linear model where methylation data is a covariate.

8) run_vc_test.R

This R script runs a variance component score test to investigate any Genotype x Methylation x Tissue interaction effect in the Gibbs et al data. 

9) final_results.R

This script aggregates all the text files produced in each sub-directory after running all the above scripts and identifies statistically significant assocations. 
