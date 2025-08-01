<h1>This directory contains the R scripts used for Mavrommatis et al. </h1>
 
 URL: 

**Script1_data_prep_7July2025.r** <br>
Includes data cleaning steps to calculate the time to censoring, impute missing values for the covariates using KNN, and calculates age acceleration residuals that additionally adjust for estimated white blood cell proportions and relatedness (kinship matrix).

**Script2_Cox_disease_7July2025.r** <br>
Cox regerssion analyses for each clock/disease pairing. We consider diagnoses made after the blood draw and over the first 10 years of follow up to electronic health records. Death is treated as a censoring event. The code cycles through each disease outcome, runs the regression for each clock (adjusting for age, sex, deprivation, BMI, smoking, education and alcohol) and tests the regression assumptions (Schoenfeld residuals). 

**Script3_logistic_AUC_7July2025.r** <br>
As per Script 2 but this time running a glm (logistic regression) for a model with covariates only and then adding an epigenetic clock. Disease status is considered over the first 10 years of follow up.  

**Script4_mortality_Cox_AUC_7July2025.r** <br>
As per scripts 2 and 3 but, instead of disease outcomes, we consider all-cause mortality.

**Script5_average_effect_size_7July2025.r** <br>
Code to generate the output in Supplementary Tables 3 and 4 (average logHR for each clock across all 174 diseases and a comparison of the average effects with GrimAge v1 as the reference category).
