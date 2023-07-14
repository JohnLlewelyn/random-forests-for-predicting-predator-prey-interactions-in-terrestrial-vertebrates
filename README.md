## Predicting predator-prey interactions in terrestrial endotherms using random forest
<img align="right" src="./network-figure.png" alt="contraception" width="400" style="margin-top: 20px">

Accompanies paper:

<a href="https://globalecologyflinders.com/people/#JL">Llewelyn, J</a>, G Strona, CR Dickman, A Greenville, G Wardle, MSY Lee, S Doherty, F Shabani, F Saltré, CJA Bradshaw. 2022. Predicting predator-prey interactions in terrestrial endotherms using random forest. <em>bioRχiv</em> doi:<a href="http://doi.org/10.1101/2022.09.02.506446">10.1101/2022.09.02.506446</a>

The files include a folder with data files (.rds format) and a folder with R scripts. File paths will need to be adjusted once downloaded, and can be found in the scripts by searching for “###”.

### DATA FOLDER
The data folder contains 9 files.:
-	<code>GloBIplus_Int20EVs.RDS</code> contains the global interaction records,
-	<code>allNon_sameCont.RDS</code> contains the global non-interaction records,
-	<code>allperms_cut2_20EVs.RDS</code> contains the interaction and non-interaction records for the seven focal predators from the <a href="https://en.wikipedia.org/wiki/Simpson_Desert">Simpson Desert</a>,
-	and 6x 'KeepVar' files that contain the variables used in the few-variable models (i.e., the most important variables as identified by the variable importance script).

Each of the interaction and non-interaction files include all ecomorphological traits and phylogenetic eigenvectors used in analyses. 

### FUNCTIONS FOLDER
The functions folder contains three files with functions required by the script files:
- (<code>all_functions_ranger.R</code>),
- (<code>opt.functions.R</code>),
- (<code>variable_importance_functions.R</code>).

### SCRIPT FOLDER
The script folder contains three more folders.
1)	OPTIMISE ON GLOBAL_APPLY TO GLOBAL AND SIMPSON DESERT
contains scripts that (i) identify optimal parameter values for 12 different random forest models (that differ in terms of the trait data used and the global training dataset) and (ii) applies the 12 optimised model (as determined in step i) to the global interaction/non-interaction data and the Simpson Desert data (for the seven focal predators).
2)	DATA QUALITY MANIPULATION AND MODEL PERFORMANCE
contains two more folders of scripts that test the effect of modifying training data quality on model performance (when training on the enhanced global data and applied to the Simpson Desert data).
<em>i</em>. RECORDREMOVAL&REPLACE_MODELPERFORMANCE
contains scripts that test the effect of removing records or switching interaction records to non-interactions (false negatives) on model performance. These modification to training data quality were made to different subsets of the data including: the whole dataset, focal prey species only, focal predator species only, and non-focal species (non-Simpson Desert) only.
<em>ii</em>.	CORRELATION&CHANGE_PROBABILITY
contains scripts testing the effects of modifying the focal-predator component of training data (removing records or switching interactions to noninteractions) on (a) relative suitability of different prey for each predator and (b) the mean probability assigned to potential prey for each predator.
3) VARIABLE_IMPORTANCE
contains 6 scripts for identifying the most important variables to retain in the few-variable models.

### SUPPORTING INFORMATION DATA FOLDER
The supporting information data folder contains 3 files.:
- <strong>[interactions_between_Simpson_Desert_species.xlsx](supporting%20%information%20data%20files/interactions_between_Simpson_Desert_species.xlsx)</strong> [.tsv](supporting%20information%20data%20files/interactions_between_Simpson_Desert_species.tsv),[.csv](supporting%20information%20data%20files/interactions_between_Simpson_Desert_species.csv) contains observed predator-prey interactions between the 7 focal predators and prey species in the Simpson Desert species assemblage (birds and mammals only).
- <strong>[Simpson_Desert_predators_with_nonSD_prey.xlsx](supporting%20information%20data%20files/Simpson_Desert_predators_with_nonSD_prey.xlsx)</strong> [.tsv](supporting%20information%20data%20files/Simpson_Desert_predators_with_nonSD_prey.tsv),[.csv](supporting%20information%20data%20files/Simpson_Desert_predators_with_nonSD_prey.csv) contains observed predator-prey interactions between the 7 focal predators and non-Simpson Desert species (birds and mammals only).
- <strong>[Simpson_Desert_sp_traits.xlsx](supporting%20information%20data%20files/Simpson_Desert_sp_traits.xlsx)</strong> [.tsv](supporting%20information%20data%20files/Simpson_Desert_sp_traits.tsv),[.csv](supporting%20information%20data%20files/Simpson_Desert_sp_traits.csv) contains trait data for birds and mammals from the Simpson Desert species assemblage.

### Indexing

[![GloBI Review by Elton](../../actions/workflows/review.yml/badge.svg)](../../actions/workflows/review.yml) [![GloBI](https://api.globalbioticinteractions.org/interaction.svg?accordingTo=globi:JohnLlewelyn/random-forests-for-predicting-predator-prey-interactions-in-terrestrial-vertebrates&refutes=true&refutes=false)](https://globalbioticinteractions.org/?accordingTo=globi:JohnLlewelyn/random-forests-for-predicting-predator-prey-interactions-in-terrestrial-vertebrates)

This repository contains configuration to help Global Biotic Interactions (GloBI, https://globalbioticinteractions.org) data associated with:

Llewelyn, J., Strona, G., Dickman, C.R., Greenville, A.C., Wardle, G.M., Lee, M.S.Y., Doherty, S., Shabani, F., Saltré, F. and Bradshaw, C.J.A. (2023), Predicting predator–prey interactions in terrestrial endotherms using random forest. Ecography e06619. https://doi.org/10.1111/ecog.06619
