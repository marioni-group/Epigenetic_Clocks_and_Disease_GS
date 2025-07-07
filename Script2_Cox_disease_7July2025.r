# Load Required Packages
packages <- c("dplyr", "ggplot2", "survival")
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

# Check number of unique diseases
cat("Number of unique diseases:", length(dis_list), "\n")

# Create storage lists for results
res.cox <- setNames(vector("list", length(clocks)), clocks)
female_proportions <- list()
age_percentiles <- list()
n_event <- list()
median_age_per_disease <- list()
res.cox.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.zph.global <- setNames(vector("list", length(clocks)), clocks)

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
  female_proportions[[disease]] <- prop_F

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
  percent_female <- prop_F * 100

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


  ## Run Cox regression model
  if (prop_F > 0.1 & prop_F < 0.9) {
    # Case 1: prop_F between 0.1 and 0.9
    for (clock in clocks) {
      model_formula <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex"))
      res.cox[[clock]][[disease]] <- coxph(model_formula, data = final_data)
    }
  } else if (prop_F <= 0.1) {
    # Case 2: prop_F less than 0.1
    for (clock in clocks) {
      model_formula <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank"))
      res.cox[[clock]][[disease]] <- coxph(model_formula, data = final_data[final_data$sex == 0, ])
    }
  } else {
    # Case 3: prop_F greater than or equal to 0.9
    for (clock in clocks) {
      model_formula <- as.formula(paste0("Surv(t_event, event) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank"))
      res.cox[[clock]][[disease]] <- coxph(model_formula, data = final_data[final_data$sex == 1, ])
    }
  }

  # Proportional Hazards Assumption test
  cat(sprintf("Running zph for %s\n", disease))
  for (clock in clocks) {
    res.cox.zph.local[[clock]][[disease]] <- cox.zph(res.cox[[clock]][[disease]])$table[1,] 
    res.cox.zph.global[[clock]][[disease]] <- cox.zph(res.cox[[clock]][[disease]])$table["GLOBAL",]
  }

  # Cleanup
  rm(final_data)
  gc() 
}


## post processing of results

local_p = list()
for(clock in names(res.cox.zph.local)){
tmp = res.cox.zph.local[[clock]]
local_p[[clock]] = data.frame(do.call('rbind', lapply(tmp, function(x){x[3]})))
}

for(clock in names(res.cox.zph.local)){
local_p[[clock]]$clock = clock
local_p[[clock]]$disease = rownames(local_p[[clock]])
}

local_p = do.call('rbind', local_p)
names(local_p)[c(1,2)] <- c("local_p", "Clock")

global_p = list()
for(clock in names(res.cox.zph.global)){
tmp = res.cox.zph.global[[clock]]
global_p[[clock]] = data.frame(do.call('rbind', lapply(tmp, function(x){x[3]})))
}

for(clock in names(res.cox.zph.global)){
global_p[[clock]]$clock = clock
global_p[[clock]]$disease = rownames(global_p[[clock]])
}

global_p = do.call('rbind', global_p)
names(global_p)[c(1,2)] <- c("global_p", "Clock")

# Write out the mean ages, % female for case/control of each disease
disease_summary_df <- do.call("rbind", disease_summary)
write.csv(disease_summary_df, file = "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/disease_meta_data_18Mar2025.csv", row.names = FALSE)


## Extract summary output from Cox models for plotting 
extract_coefficients <- function(models_list) {
  # Initialize an empty list to store results
  results <- list()
  
  # Iterate over each disease
  for (disease in names(models_list)) {
    # Initialize a list to store coefficients for the current disease
    disease_coefficients <- list()
    
    # Iterate over each epigenetic clock model
    for (clock in names(models_list[[disease]])) {
      model <- models_list[[disease]][[clock]]
      
        # Extract the first row of coefficients
        coeff <- summary(model)$coefficients[1, ]
	n_event <- summary(model)$nevent
	n_total <- summary(model)$n
        disease_coefficients[[clock]] <- c(coeff, n_event, n_total)
    }
    
    # Combine the coefficients into a data frame for the current disease
    results[[disease]] <- do.call(rbind, disease_coefficients)
  }
  
  # Convert the results list to a data frame
  final_results <- do.call(rbind, lapply(names(results), function(disease) {
    cbind(Disease = disease, results[[disease]])
  }))
  
  return(as.data.frame(final_results))
}

out <- extract_coefficients(res.cox)
names(out)[c(1,3,4,6,7,8)] <- c("Clock","HR", "SE","P","N_cases","N_total")
out$disease <- sub("\\..*", "", rownames(out))

out <- out %>%
  mutate(across(c("HR", "coef", "SE", "P", "N_cases", "N_total"), as.numeric))

out$LCI <- exp(out$coef - 1.96*out$SE)
out$UCI <- exp(out$coef + 1.96*out$SE)

out1 <- out[,c("Clock","disease","N_cases","N_total","HR","LCI","UCI","P")]

out2 <- merge(out1, local_p, by=c("Clock", "disease")) 
out3 <- merge(out2, global_p, by=c("Clock", "disease")) 

out4 <- out3[out3$N_cases > 29, ]

write.csv(out4, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/All_Results_12Nov2024.csv", row.names=F)