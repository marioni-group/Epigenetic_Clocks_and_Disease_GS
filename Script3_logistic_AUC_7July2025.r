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
    }
  }

  # Cleanup
  rm(final_data)
  gc() 
}


### compile results ###
auc_full <- as.data.frame(do.call('rbind', auc_est))
auc_null <- as.data.frame(do.call('rbind', auc_null))

df1_num <- as.data.frame(lapply(auc_full, unlist))
df2_num <- as.data.frame(lapply(auc_null, unlist))
result <- df1_num - df2_num

df1 <- as.data.frame(t(df1_num))
df2 <- as.data.frame(t(df2_num))
auc_diff <- as.data.frame(t(result))

df1$disease <- rownames(df1)
df2$disease <- rownames(df2)
auc_diff$disease <- rownames(auc_diff)

df1 <- df1[,c("disease", setdiff(names(df1), "disease"))]
df2 <- df2[,c("disease", setdiff(names(df1), "disease"))]
auc_diff <- auc_diff[,c("disease", setdiff(names(auc_diff), "disease"))]

write.csv(df1, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_full_16NMay2025.csv", row.names=F)
write.csv(df2, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_null_16NMay2025.csv", row.names=F)
write.csv(auc_diff, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_diff_16NMay2025.csv", row.names=F)


## create single data frame

names(df1) <- paste0(names(df1), "_auc2_full")
names(df2) <- paste0(names(df2), "_auc1_null")
names(auc_diff) <- paste0(names(auc_diff), "_auc3_diff")

df3 <- cbind(df1, df2, auc_diff)
df3 <- df3[, order(names(df3))]
df3 <- df3[,-c(2:3)]
names(df3)[1] <- "disease"

write.csv(df3, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Disease_AUC_all_16NMay2025.csv", row.names=F)