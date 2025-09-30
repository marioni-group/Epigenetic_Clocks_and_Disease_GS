# Load Required Packages
packages <- c("dplyr")
lapply(packages, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) {
    install.packages(x)
  }
  library(x, character.only = TRUE)
})

# Read in clock data
clock_data <- readRDS("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/Merged_clock_data_17Oct2024.RDS")

# Read in diseases
diseases <- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-03-29_diseases.csv", stringsAsFactors = FALSE)
diseases <- diseases[-which(diseases$Source == "Primary_Care" & diseases$GP_Consent==0),]

# Get unique diseases
dis_list <- unique(diseases$Disease)

# Check number of unique diseases
cat("Number of unique diseases:", length(dis_list), "\n")


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
  
  # Filter for incident cases
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

  # Summary statistics
  mean_age <- mean(final_data$age[final_data$event==1], na.rm = TRUE)
  sd_age <- sd(final_data$age[final_data$event==1], na.rm = TRUE)
  ncontrol <- nrow(final_data) - n_event
  percent_female <- length(which(final_data$sex == 1 & final_data$event ==1))/ length(which(final_data$event ==1 )) * 100  

  mean_age_ctrl <- mean(final_data$age[final_data$event==0], na.rm = TRUE)
  sd_age_ctrl <- sd(final_data$age[final_data$event==0], na.rm = TRUE)
  percent_female_ctrl <- length(which(final_data$sex == 1 & final_data$event ==0))/ length(which(final_data$event ==0 )) * 100  

  mean_age_event <- mean(final_data$age[final_data$event==1] + final_data$t_event[final_data$event==1], na.rm = TRUE)
  sd_age_event <- sd(final_data$age[final_data$event==1] + final_data$t_event[final_data$event==1], na.rm = TRUE)
  mean_tte <- mean(final_data$t_event[final_data$event==1], na.rm = TRUE)
  sd_tte <- sd(final_data$t_event[final_data$event==1], na.rm = TRUE)

  # Append to summary list
  disease_summary[[disease]] <- data.frame(
    disease = disease,
    ncase = n_event,
    ncontrol = ncontrol,
    mean_age_case = mean_age,
    sd_age_case = sd_age,
    mean_age_disease = mean_age_event,
    sd_age_disease = sd_age_event,
    mean_time_to_disease = mean_tte,
    sd_time_to_disease = sd_tte,
    mean_age_control = mean_age_ctrl,
    sd_age_control = sd_age_ctrl,
    percent_female_case = percent_female,    
    percent_female_control = percent_female_ctrl
  )


  # Cleanup
  rm(final_data)
  gc() 
}


# Write out the mean ages, % female for case/control of each disease
disease_summary_df <- do.call("rbind", disease_summary)

write.csv(disease_summary_df, file = "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/disease_meta_data_21Aug2025.csv", row.names = FALSE)

