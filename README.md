<h1>This directory contains the R scripts used for Mavrommatis et al. </h1>
 
 Pre-print URL: https://www.medrxiv.org/content/10.1101/2025.07.14.25331494v1 <br>
 Published paper URL: 

**Script1_data_prep_25Aug2025.r** <br>
Includes data cleaning steps to calculate the time to censoring, impute missing values for the covariates using KNN, and calculates age acceleration residuals that additionally adjust for estimated white blood cell proportions and relatedness (kinship matrix).

**Script1a_disease_meta_data_25Aug2025.r** <br>
Calculates the descriptive statistics for each disease 

**Script2_Cox_disease_15Aug2025.r** <br>
Cox regerssion analyses for each clock/disease pairing. We consider diagnoses made after the blood draw and over the first 10 years of follow up to electronic health records. Death is treated as a censoring event. The code cycles through each disease outcome, runs the regression for each clock (adjusting for age, sex, deprivation, BMI, smoking, education and alcohol) and tests the regression assumptions (Schoenfeld residuals). Basic models adjust for age and sex, fully adjusted models include all covariates. Sex and smoking-stratified and interaction models are also considered.

**Script3_logistic_AUC_13Aug2025.r** <br>
Glm (logistic regression) for a model with covariates only and then adding an epigenetic clock. Disease status is considered over the first 10 years of follow up.  

**Script4_mortality_Cox_AUC_25Aug2025.r** <br>
As per scripts 2 and 3 but, instead of disease outcomes, we consider all-cause mortality.

**Script5_average_effect_size_7July2025.r** <br>
Code to generate average logHR for each clock across all 174 diseases and a comparison of the average effects with GrimAge v1 as the reference category.

**Tables_and_Figures_code_25Aug2025.r** <br>
Code to generate Supplementary Tables and Figures 1 and 2.
