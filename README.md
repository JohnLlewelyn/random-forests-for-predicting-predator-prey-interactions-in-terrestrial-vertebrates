## Predicting predator-prey interactions in terrestrial endotherms using random forest
<img align="right" src="network figure.jpg" alt="contraception" width="400" style="margin-top: 20px">

Accompanies paper:

<a href="https://globalecologyflinders.com/people/#JL">Llewelyn, J</a>, G Strona, CR Dickman, A Greenville, G Wardle, MSY Lee, S Doherty, F Shabani, F Saltré, CJA Bradshaw. 2022. Predicting predator-prey interactions in terrestrial endotherms using random forest. <em>bioRχiv</em> doi:<a href="http://doi.org/10.1101/2022.09.02.506446">10.1101/2022.09.02.506446</a>

The files include a folder with data files (.rds format) and a folder with R scripts. File paths will need to be adjusted once downloaded, and can be found in the scripts by searching for “###”.

### DATA FOLDER
The data folder contains 4 files.:
-	GloBIplus_Int20EVs.RDS contains the enhanced global interaction records,
-	allNon_sameCont_1.RDS and allNon_sameCont_2.RDS contain the enhanced global non-interaction records – these two datasets need to be combined (rbind) and renamed allNon_sameCont.RDS or Non for running scripts,
-	allperms_cut2_20EVs.RDS contains the interaction and non-interaction records for the seven focal predators from the Simpson Desert.
Each of these files include all ecomorphological traits and phylogenetic eigenvectors used in analyses. 

### FUNCTIONS FOLDER
The functions folder contains a file (all_functions_ranger.R) with functions required by the script files.

### SCRIPT FOLDER
The script folder contains three more folders.
1)	OPTIMIZATION OF RANDOM FOREST MODELS
contains scripts that identify optimal parameter values for 6 different random forest models (that differ in terms of the trait data used).
2)	APPLIED TO GLOBAL AND SIMPSON DESERT DATASETS
contains scripts that apply the 6 optimised model (as determined in step 1) to the enhanced global interaction/non-interaction data and the Simpson Desert data (for the seven focal predators).
3)	DATA QUALITY MANIPULATION AND MODEL PERFORMANCE
contains two more folders of scripts that test the effect of modifying training data quality on model performance (when training on the enhanced global data and applied to the Simpson Desert data).
<em>i</em>. RECORDREMOVAL&REPLACE_MODELPERFORMANCE
contains scripts that test the effect of removing records or switching interaction records to non-interactions (false negatives) on model performance. These modification to training data quality were made to different subsets of the data including: the whole dataset, focal prey species only, focal predator species only, and non-focal species (non-Simpson Desert) only.
<em>ii</em>.	CORRELATION&CHANGE_PROBABILITY
contains scripts testing the effects of modifying the focal-predator component of training data (removing records or switching interactions to noninteractions) on (a) relative suitability of different prey for each predator and (b) the mean probability assigned to potential prey for each predator.

### SUPPORTING INFORMATION DATA FOLDER
The supporting information data folder contains 3 files.:
- interactions_between_Simpson_Desert_species.xlsx contains observed predator-prey interactions between the 7 focal predators and prey species in the Simpson Desert species assemblage (birds and mammals only).
- Simpson_Desert_predators_with_nonSD_prey.xlsx contains observed predator-prey interactions between the 7 focal predators and non-Simpson Desert species (birds and mammals only).
- Simpson_Desert_sp_traits.xlsx contains trait data for birds and mammals from the Simpson Desert species assemblage.
