# Load Required Packages
packages <- c("dplyr", "survival", "pROC")
lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
  }
  library(x, character.only = TRUE)
})

# Read in clock data
clock_data <- readRDS("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/Merged_clock_data_17Oct2024.RDS")

# List of clocks to loop through
clocks <- c("Horvathv1", "Hannum", "Lin", "PhenoAge", "YingCausAge", 
            "YingDamAge", "YingAdaptAge", "Horvathv2", "Zhang_10", 
            "DunedinPoAm38", "DunedinPACE", "DNAmGrimAge", 
            "DNAmGrimAge.1", "DNAmTL")

# Read in diseases
diseases <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-03-29_diseases.csv", stringsAsFactors = FALSE)
diseases <- diseases[-which(diseases$Source == "Primary_Care" & diseases$GP_Consent==0),]

# Get unique diseases
dis_list <- unique(diseases$Disease)

# Create storage lists for results
auc_est <- setNames(vector("list", length(clocks)), clocks)
auc_null <- setNames(vector("list", length(clocks)), clocks)
auc_p_comp <- setNames(vector("list", length(clocks)), clocks)

# Disease summary table
disease_summary <- list()

# Loop through each disease
for (disease in dis_list) {
  disease_data <- diseases %>% filter(Disease == disease)
  cat(sprintf("Processing disease %s\n", disease))
  
  # Calculate time to incident
  disease_data <- disease_data %>%
    mutate(yod = as.numeric(substr(dt1_ym, 1, 4)),
           mod = as.numeric(substr(dt1_ym, 5, 6)),
           yoa = as.numeric(substr(gs_appt, 1, 4)),
           moa = as.numeric(substr(gs_appt, 5, 6)),
           t_disease = (yod - yoa) + ((mod - moa) / 12)) %>%
    select(-mod, -yoa, -moa)  # Remove excess columns

  # Merge with clock data
  merged_data <- merge(disease_data, clock_data, by = "id", all.y = TRUE)
  
  # Filter to incident cases
  merged_data_2 <- merged_data %>% filter(incident != 0 | is.na(incident))
  merged_data_2 <- merged_data_2 %>% filter(t_disease > 0 | is.na(t_disease))
  merged_data_2 <- merged_data_2 %>%
    mutate(DP = ifelse(t_disease > 0, 1, 0)) %>%
    filter(is.na(yod) | yod > 1980) 
 
  # Order data and remove duplicates
  merged_data_ord <- merged_data_2[order(merged_data_2$dt1_ym), ]
  merged_data_ord <- merged_data_ord[!duplicated(na.omit(merged_data_ord$id)), ]
  
  # Create event column and time to event
  merged_data_ord <- merged_data_ord %>%
    mutate(event = ifelse(is.na(incident), 0, 1),
           t_event = ifelse(event == 1, t_disease, t_censor))

  # Trim data
  final_data <- merged_data_ord %>% select(-yod, -dt1_ym, -gs_appt)

  # Skip disease if no events
  if (sum(final_data$event == 1) == 0) next

  # Gender filtering
  prop_F <-  length(which(final_data$sex == 1 & final_data$event ==1))/ length(which(final_data$event ==1 ))
  if (prop_F > 0.9) {
    final_data <- final_data %>% filter(sex == 1)
  } else if (prop_F < 0.1) {
    final_data <- final_data %>% filter(sex == 0)
  }

  # Filter to 10 year analyses
  final_data$event[final_data$t_event > 10 & final_data$event == 1] <- 0
  final_data$t_event[final_data$t_event > 10] <- 10 

  # Count events
  n_event <- sum(final_data$event == 1)
  final_data$n_event <- n_event

  # Skip if n_event < 30
  if (n_event < 30) next



  ## Run logistic regression model
  if (prop_F > 0.1 & prop_F < 0.9) {
    # Case 1: prop_F between 0.1 and 0.9
    for (clock in clocks) {
      model_formula <- as.formula(paste0("event ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex"))
      m1 <- glm(model_formula, data = final_data)
      probs <- predict(m1, type="response")
      roc_obj <- roc(final_data$event, probs)
      auc_est[[clock]][[disease]] <- auc(roc_obj) 
      
      model_formula2 <- as.formula(paste0("event ~ age + bmi + years + pack_years + units + rank + sex"))
      m2 <- glm(model_formula2, data = final_data)
      probs2 <- predict(m2, type="response")
      roc_obj2 <- roc(final_data$event, probs2)
      auc_null[[clock]][[disease]] <- auc(roc_obj2)  
	  
	  auc_p_comp[[clock]][[disease]] <- roc.test(roc_obj, roc_obj2)$p.value
	  
      }
  } else if (prop_F <= 0.1) {
    # Case 2: prop_F less than 0.1
    for (clock in clocks) {
      model_formula <- as.formula(paste0("event ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank"))
      m1 <- glm(model_formula, data = final_data)
      probs <- predict(m1, type="response")
      roc_obj <- roc(final_data$event, probs)
      auc_est[[clock]][[disease]] <- auc(roc_obj) 
      
      model_formula2 <- as.formula(paste0("event ~ age + bmi + years + pack_years + units + rank"))
      m2 <- glm(model_formula2, data = final_data)
      probs2 <- predict(m2, type="response")
      roc_obj2 <- roc(final_data$event, probs2)
      auc_null[[clock]][[disease]] <- auc(roc_obj2) 
	  
	  auc_p_comp[[clock]][[disease]] <- roc.test(roc_obj, roc_obj2)$p.value
    }
  } else {
    # Case 3: prop_F greater than or equal to 0.9
    for (clock in clocks) {
      model_formula <- as.formula(paste0("event ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank"))
      m1 <- glm(model_formula, data = final_data)
      probs <- predict(m1, type="response")
      roc_obj <- roc(final_data$event, probs)
      auc_est[[clock]][[disease]] <- auc(roc_obj) 
      
      model_formula2 <- as.formula(paste0("event ~ age + bmi + years + pack_years + units + rank"))
      m2 <- glm(model_formula2, data = final_data)
      probs2 <- predict(m2, type="response")
      roc_obj2 <- roc(final_data$event, probs2)
      auc_null[[clock]][[disease]] <- auc(roc_obj2) 
	  
	  auc_p_comp[[clock]][[disease]] <- roc.test(roc_obj, roc_obj2)$p.value
    }
  }

  # Cleanup
  rm(final_data)
  gc() 
}


### compile results ###
auc_full <- as.data.frame(do.call('rbind', auc_est))
auc_null <- as.data.frame(do.call('rbind', auc_null))
auc_p_comp <- as.data.frame(do.call('rbind', auc_p_comp))

df1_num <- as.data.frame(lapply(auc_full, unlist))
df2_num <- as.data.frame(lapply(auc_null, unlist))
result <- df1_num - df2_num
p_comp <- as.data.frame(lapply(auc_p_comp, unlist))

df1 <- as.data.frame(t(df1_num))
df2 <- as.data.frame(t(df2_num))
auc_diff <- as.data.frame(t(result))
p_comp1 <- as.data.frame(t(p_comp))

df1$disease <- rownames(df1)
df2$disease <- rownames(df2)
auc_diff$disease <- rownames(auc_diff)
p_comp1$disease <- rownames(p_comp1)

df1 <- df1[,c("disease", setdiff(names(df1), "disease"))]
df2 <- df2[,c("disease", setdiff(names(df1), "disease"))]
auc_diff <- auc_diff[,c("disease", setdiff(names(auc_diff), "disease"))]
p_comp1 <- p_comp1[,c("disease", setdiff(names(p_comp1), "disease"))]

write.csv(df1, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_full_16NMay2025.csv", row.names=F)
write.csv(df2, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_null_16NMay2025.csv", row.names=F)
write.csv(auc_diff, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_diff_16NMay2025.csv", row.names=F)
write.csv(p_comp1, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_p_comp1_13Aug2025.csv", row.names=F)


d1 <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_full_16NMay2025.csv")
d2 <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_null_16NMay2025.csv")
d3 <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_diff_16NMay2025.csv")
d4 <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_p_comp1_13Aug2025.csv")


# Function to transform each data frame into the desired format
transform_df <- function(df) {
  # Ensure the data frame is a matrix for consistent handling
  df <- as.matrix(df)
  
  # Create the repeated first column
  first_col_repeated <- rep(df[, 1], each = 14)
  
  # Create the second column (flatten the matrix and convert to numeric)
  second_col <- as.numeric(as.vector(t(df[, -1])))
  
  # Create the repeated column names vector (repeating 175 times for each column)
  col_names_repeated <- rep(colnames(df)[2:ncol(df)], times = 175)
  
  # Create the final data frame
  result_df <- data.frame(first_col_repeated, second_col, col_names_repeated)
  
  return(result_df)
}

# Apply the function to each data frame (d1, d2, d3)
d11 <- transform_df(d1)
d22 <- transform_df(d2)
d33 <- transform_df(d3)
d44 <- transform_df(d4)

d4 <- cbind(d22[,1:2], d11[,2], d33[,2], d44[,2:3])
names(d4) <- c("disease","auc_null","auc_full","auc_diff","P_auc_diff","clock")
d4 <- d4[,c(1,6,2:5)]

write.csv(d4, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_long_13Aug2025.csv", row.names=F)


auc_results <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_long_13Aug2025.csv")
hrs <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/All_Results_22Aug2025.csv")
hrs <- hrs[hrs$model=="full",]
names(auc_results)[2] <- "Clock"

final_out <- merge(hrs, auc_results, by=c("disease","Clock"))

write.csv(final_out, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_HR_AUC_merged_25Aug2025.csv", row.names=F)

