-------
predicting_cognition_cortical_morphometry_functional_networks
-------

__Note__

- This repository only contains the `scripts/` directory from my local git repo but I have included descriptions of all files and subdirectories for better understanding of my workflow.

- The input is a set of subjects run through the Freesurfer reconstruction stream (recon-all) and resampled to fsaverage using the recon-all -qcache flag in addition to subject data in csv file with cognitive scores (ie.., memory and executive function ability scores and covariates of age and gender).

__________

## Directories (with file contents)
__________

`manuscript_inprep.(Rmd|pdf)` : the working version of markdown manuscript displaying results. Much of the code for producing R figures and tables is taken from individual R scripts. All figures created in python are loaded as PNG files.
To view document, see the __pdf__ version.


### scripts/

__The following scripts are in order of workflow__
1. `make_npdata_from_pca.R` : runs principal component analysis from individual tasks and makes composite scores.

2. `make_corticalmorph_dfs.(ipynb|py)`: this notebook uploads each subject's morphometry surface files for surface area and thickness in fsaverage space and (1) concatenates all subjects, (2) joins by an annotation (e.g., Yeo et al., 2011 networks), (3) saves a file for each annotation label (e.g., each network in Yeo et al., 2011) containing a dataframe for each measure/network combination. In other
words, it saves a concatenated pandas dataframes (ie., subjects (rows) by all cortical morphometry
measures (surface area and thickness) for vertices across both hemispheres (columns) in pickle files. At the end of the script, it also computes whole brain thickness and surface area and adds it to
npdata dataframe so whole brain predictions and use as covariate is easier than going through entire
pipeline.

3. `cv_results.(ipynb|py)` :  uses outputs of `make_corticalmorph_dfs` and runs cross validation
procedure
using sci-kit learn tools---this script uses a set of train/test set and iterates through these
iterations/splits and, in parallel, outputs prediction values for each iteration/split and results
(e.g., r squared and mean squared error metrics) in a pandas dataframe for different feature
selection
thresholds/methods and different cortical morphometry measures/networks of interest (i.e., percentile
and false positive rate sci kit learn options). For the percentile feature selection method, only the
"max percentile" was used from a feature selection search (file with these values was created from
cv_results_percentile_prelim.R script). 

4. `make_cv_dfs.R` : functions that make nested dataframe used for analyses in cv_predict consisting of necessary
prediction values obtained from cv_results_*_.py and npdata.

5. `cv_predict.R` : runs bootstrap statistics on predictions obtained from "cv_results_*_.ipynb" and
produces figures. Currently contains other tables included in manuscript as well. The code for figures and tables is also in `manuscript_inprep` to increase ease of reproducibility.

6. `regression_mediation_pervertex_wholebrain.ipynb` : Inputs the whole brain
concatenated pandas dataframe from `make_corticalmorph_dfs.(ipynb|py)` and computes either a vertex-wise linear regression
or mediation for each vertex after resampling (with replacement; bootstrap samples) and keeps a
running count of the mean of the either the beta coefficient (if per vertex regression) or indirect
effect (if per vertex mediation) and standard error. Then bootstrap ratio scores are calculated.
Freesurfer overlay files are created for visualization and further analyses (i.e., cluster correction
using mri_cluster). Z-scores are computed by shuffling annotation labels to create a null
distribution. That is, given the number of significant vertices obtained, how many, by chance, should
each network obtain.

7. `visualize_numvertex_stats_wholebrain.R` : Plots the bootstrap statistics (with
shuffled stats) obtained from `regression_mediation_pervertex_wholebrain.ipynb`

### curvoverlays/ : contains the outputted freesurfer curv files for visualization of each hemisphere statistics including
  1. bootstrap ratio scores (which were inputted into mri_surfcluster at specified cut offs)
  2. p value average and percent selected across splits/iterations
### cv_splits/ : contains indices for splits across iterations (i.e., pickle file created from balanced_cv fxn from cv_results*.ipynb)

### data/ :
1. output from the python scripts for cv results and prediction values which were inputted into R
scripts for stats/figures
2. wholebrain_bootstrap/
    - pickle files are the latest output of wholebrain bootstrap script (regression and per vertex mediation) in  dir of scripts/Cortical_Thickness_ML_1_19_18_BOOTSTRAP_MEDIATION_pervertex.ipynb
    - wholebrain_bootstrap_bootscore_permtests_frompearsoncoeffs.csv are the z-score tests from using pearson correlation coefficient as bootstrapped variable in prelim script version (use for checking newer scripts)
  
### figures/ 
  1. contains output of cv_predict.R (tables, stats, and figures)

### old/
  1. old files from early development (just in case I need them)
  
### reconall/
  1. initial files used to run Freesurfer processing (i.e., recon-all)

__________
## Branches
_________

***Note see README in prelimV2 and prelimV3_leaveoneout for reasoning of major switch in methods/approach***

- dev : (working branch) used as the step before merging into the master branch 
- prelim : files used for my Chapter 1 prelim as presented at preliminary examination and document (see prelim_files folder on server just in case as well for same thing).
- prelimV2 : script changes after prelim exam but before major changes to cv_results approach (documented below)
- prelimV3_leaveoneout : script fixes and change to leave one out predictions (figures revealed same results as KFold previously used and less computationally intensive given that I did repeated KFold and now will do permutation testing in upcoming version)
- covars_leaveoneoutV4 : uses leaveoneout cross validation, has the no intercept option for correlations. However, found that there was high variance in MSE for leave one out so may not be best for comparing models based on MSE.

##Decisions made
-  scale y and x to fix intercept to zero (Poldrack blog)
-  run permutation tests to determine null distribution, get pvals, get measure of bias (see all papers that do this)
    - pvals of r stat are the proportion of permutations greater than or equal to empirical stat
- mediation analyses: using prediction values from cv pipeline, switching sign of Age (as Age is negatively correlated with cognition to make consistent with predictive ability-- i.e., negative predictive ability means worse than chance whereas positive is better than chance). Doing this approach rather than more complex approach of train-test split.

**taken out:**

-  tentative: not assessing mse of nested models as unclear how to run inferential tests (keeping it simpler and just doing rs with permutation testing to keep more similar to prev papers)
    - This two pronged approach was done in a 2011 JoCN paper by Tor Wager (see B--the "multiple regression analysis"-- and C-- the correlations---below from this paper). 
    - Still have this in code though
  





