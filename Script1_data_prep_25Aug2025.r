### Load Required Packages ###
packages <- c("dplyr", "coxme", "stringr")
lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
  }
  library(x, character.only = TRUE)
})

## if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## BiocManager::install("impute")
library(impute)


### load in Epigenetic Clocks ###
clock_data <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/gs_clocks.csv")

### load in disease, Covariate and Mortality data  ###
disease_data <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-03-29_diseases.csv")
cov_data <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-04-19_covariates.csv")
mortality_data <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-04-19_deaths.csv")

### Load in file linking ids to meth dis ###
GSK_MAP <- readRDS("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/GS20k_Targets_18869.rds")
names(GSK_MAP)[1] <- "DNAm_ID"

### load in appt records ###
appt_rds <- readRDS("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/202-04-26_gs_appt.rds")


## Filter Clocks to relevant predictors ##
names(clock_data)[1] <- "DNAm_ID"
filt_clock_data<- select(clock_data, c(DNAm_ID, Horvathv1, Hannum, Lin, PhenoAge, YingCausAge, YingDamAge,
									   YingAdaptAge, Horvathv2, Zhang_10, DunedinPoAm38, DunedinPACE, 
									   DNAmGrimAge, DNAmGrimAge.1, DNAmTL))


### merge appt data to GSK ###
appt_GSK <- merge(appt_rds, GSK_MAP[,c("DNAm_ID","Sample_Name")], by.x= "id", by.y= "Sample_Name")

merged_clock_data_1 <- merge(filt_clock_data,appt_GSK, by = "DNAm_ID")

### merge merged_clock_data with mortality data ### # if columns are same, no need for by.x or by.y just by = id ##
merged_clock_data_2<- merge(merged_clock_data_1, mortality_data[,2:3], by.x= "id", by.y= "id", all.x= TRUE )

#Check NA= alive (not in mortality data)# 
merged_clock_data_2_alive <- merged_clock_data_2 %>% 
  filter(!is.na(dod_ym))

### merge in covariates cov ###
merged_clock_data_3 <- merge(merged_clock_data_2, cov_data[,2:11], by.x= "id", by.y= "id")

### Calculate Time to censor: ###

### Substr APPT dates ###
merged_clock_data_3$GS_yoa <- substr(merged_clock_data_3$ym, 1,4) %>% as.numeric()
merged_clock_data_3$GS_moa <- substr(merged_clock_data_3$ym, 5,6) %>% as.numeric()

### substr YOD dates ### 
merged_clock_data_3$GS_yod <- substr(merged_clock_data_3$dod_ym, 1,4) %>% as.numeric()
merged_clock_data_3$GS_mod <- substr(merged_clock_data_3$dod_ym, 5,6) %>% as.numeric()

# recheck Check NA= as people have not died# 
merged_clock_data_3_alive <- merged_clock_data_3 %>% 
  filter(!is.na(GS_yod))

### new column : 202204- GS APPT ###
merged_clock_data_3$cutoff_minus_GSAPPT <- (2022 - merged_clock_data_3$GS_yoa) + ((04- merged_clock_data_3$GS_moa)/12) %>% as.numeric()

### new column : dod- GSappt ### 
merged_clock_data_3$DOD_minus_GSAPPT <- (merged_clock_data_3$GS_yod - merged_clock_data_3$GS_yoa) + ((merged_clock_data_3$GS_mod - merged_clock_data_3$GS_moa)/12) %>% as.numeric()

## New column T_Censor ##
merged_clock_data_3$t_censor <- ifelse(!is.na(merged_clock_data_3$dod_ym), merged_clock_data_3$DOD_minus_GSAPPT,  merged_clock_data_3$cutoff_minus_GSAPPT)

##  Create new DF trimmed  ###
merged_clock_data_4 <- subset(merged_clock_data_3, select = - c(ym, dod_ym, GS_yoa, GS_moa, GS_yod, GS_mod, cutoff_minus_GSAPPT, DOD_minus_GSAPPT))

### read in and merge WBC proportions ###
wbc <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/gs_cell_prop.csv")
wbc1 <- wbc[,2:ncol(wbc)]
for(i in 1:ncol(wbc1)){
wbc1[,i] <- as.numeric(wbc1[,i])
}

wbc2 <- t(wbc1)
colnames(wbc2) <- wbc$X
wbc3 <- as.data.frame(wbc2)

wbc3$DNAm_ID <- as.character(str_sub(rownames(wbc3), start=2))

### read in the kinship matrix ###
kinship <- readRDS("/Cluster_Filespace/Marioni_Group/Josie/Proteins/ewas/input_files/kinship_matrix_using_fixed_2022-01-28-pedigree.rds")

### Merge WBCs and calculate AgeAccel residuals for each clock ###
merged_clock_data_5 <- merge(merged_clock_data_4, wbc3, by="DNAm_ID", all.x=T)

### regress out WBCs and kinship from each clock ###
for (i in 3:16){
merged_clock_data_5[,i] <- resid(lmekin(merged_clock_data_5[,i] ~ age + Bmem + 
								Bnv + CD4mem + CD4nv + CD8mem + CD8nv + Eos + 
								Mono + NK + Neu + Treg + (1|merged_clock_data_5$id), 
								varlist = kinship*2, data=merged_clock_data_5, na.action="na.exclude"))
}

merged_clock_data_6 <- merged_clock_data_5[,-c(22,25,27:38)]

##########################################
### Sensitivity Analysis - WBC EpiDISH ###
##########################################

### Calculate WBC proportions via EpiDish ###
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EpiDISH")

library(EpiDISH)

GS_beta <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/Biolearn_GS_18860_17April2025_REM.RDS")
row.names(GS_beta) <- GS_beta$V1 
GS_beta <- GS_beta[,-1]
WBC <- epidish(beta.m = GS_beta, ref.m = cent12CT.m, method = "RPC")$estF
saveRDS(WBC, file="/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/EpiDish_WBCs_25Aug2025_REM.rds")

### compare clock residuals with those after regressing out WBCs estimated via EpiDish ###

WBC_epidish <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/Standardised_Data/EpiDish_WBCs_25Aug2025_REM.rds")
WBC_epidish <- as.data.frame(WBC_epidish)
WBC_epidish$DNAm_ID <- row.names(WBC_epidish)

merged_clock_data_5a <- merge(merged_clock_data_4, WBC_epidish, by="DNAm_ID", all.x=T)

# regress out EpiDish WBCs and kinship from each clock ###
for (i in 3:16){
merged_clock_data_5a[,i] <- resid(lmekin(merged_clock_data_5a[,i] ~ age + 	
                                 CD4Tnv + CD4Tmem + Bmem + Bnv + Treg + CD8Tmem + CD8Tnv + Eos + NK + Neu + Mono
                                 + (1|merged_clock_data_5a$id), 
								varlist = kinship*2, data=merged_clock_data_5a, na.action="na.exclude"))
}

for (i in 3:16) {
  correlation <- cor.test(merged_clock_data_5[, i], merged_clock_data_5a[, i])$estimate
  var_name <- names(merged_clock_data_5)[i]
  print(paste0("Correlation for ", var_name, ": ", round(correlation, 4)))
}

# "Correlation for Horvathv1: 0.9837"
# "Correlation for Hannum: 0.9641"
# "Correlation for Lin: 0.9866"
# "Correlation for PhenoAge: 0.9838"
# "Correlation for YingCausAge: 0.99"
# "Correlation for YingDamAge: 0.9933"
# "Correlation for YingAdaptAge: 0.9099"
# "Correlation for Horvathv2: 0.989"
# "Correlation for Zhang_10: 0.9092"
# "Correlation for DunedinPoAm38: 0.9729"
# "Correlation for DunedinPACE: 0.975"
# "Correlation for DNAmGrimAge: 0.9785"
# "Correlation for DNAmGrimAge.1: 0.9766"
# "Correlation for DNAmTL: 0.9556"

###################################
### End of sensitivity analysis ###
###################################

### KNN imputation for missing covariates ###

covs <- as.data.frame(merged_clock_data_6[,c("age","sex","rank","bmi","years","pack_years","units")])
covs$sex <- ifelse(covs$sex=="F", 1, 0)

covs$sex <- scale(covs$sex)
covs$age <- scale(covs$age)
covs$rank <- scale(covs$rank)
covs$bmi <- scale(log(covs$bmi))
covs$years <- scale(covs$years)
covs$pack_years <- scale(log(1 + covs$pack_years))
covs$units <- scale(log(1 + covs$units))

set.seed(1.234)
covs_imp <- as.data.frame(impute.knn(as.matrix(covs), rowmax=0.6)$data)

table(rownames(covs_imp) == rownames(merged_clock_data_6))

### merge imputed covariates with clocks ###
merged_clock_data_7 <- cbind(merged_clock_data_6[,c(1:18,24)], covs_imp[,3:7])

### binary code for sex ###
merged_clock_data_7$sex <- ifelse(merged_clock_data_7$sex=="F", 1, 0)

saveRDS(merged_clock_data_7, file="/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/Merged_clock_data_17Oct2024.RDS")