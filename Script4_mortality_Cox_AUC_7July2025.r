### read in death data and QCd clock data ###
mortality_data<- read.csv("/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/2024-04-19_deaths.csv")
d <- readRDS(file="/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Data/Merged_clock_data_17Oct2024.RDS")

### create 10-year death phenotype and tte ###
temp <- which(d$id %in% mortality_data$id)
d$dead <- 0
d$dead[temp] <- 1
temp2 <- which(d$dead==1 & d$t_censor>10)  
d$dead[temp2] <- 0
d$t_censor[d$t_censor>10] <- 10

# List of clocks to loop through
clocks <- c("Horvathv1", "Hannum", "Lin", "PhenoAge", "YingCausAge", 
            "YingDamAge", "YingAdaptAge", "Horvathv2", "Zhang_10", 
            "DunedinPoAm38", "DunedinPACE", "DNAmGrimAge", 
            "DNAmGrimAge.1", "DNAmTL")

# objects to store results
res.cox <- setNames(vector("list", length(clocks)), clocks)
res.cox.zph.local <- setNames(vector("list", length(clocks)), clocks)
res.cox.zph.global <- setNames(vector("list", length(clocks)), clocks)

library(survival)

# Cox Regression model
for(clock in clocks){
      model_formula <- as.formula(paste0("Surv(t_censor, dead) ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex"))
      res.cox[[clock]] <- coxph(model_formula, data = d)
      res.cox.zph.local[[clock]] <- cox.zph(res.cox[[clock]])$table[1,] 
      res.cox.zph.global[[clock]] <- cox.zph(res.cox[[clock]])$table["GLOBAL",]
}

## prop hazards (local) ##
local_p = list()

for(clock in names(res.cox.zph.local)){
local_p[[clock]] <- res.cox.zph.local[[clock]][3]
}

local_p = as.data.frame(do.call('rbind', local_p))
local_p$Clock <- as.character(rownames(local_p))
names(local_p)[1] <- "local_p"

## prop hazards (global) ##
global_p = list()

for(clock in names(res.cox.zph.global)){
global_p[[clock]] <- res.cox.zph.global[[clock]][3]
}

global_p = as.data.frame(do.call('rbind', global_p))
global_p$Clock <- as.character(rownames(global_p))
names(global_p)[1] <- "global_p"


## Extract summary output from Cox models for plotting 
extract_coefficients <- function(models_list) {
  # Initialize an empty list to store results
  results <- list()
  mort_coefficients <- list()
   
    # Iterate over each epigenetic clock model
    for (clock in names(models_list)) {
      model <- models_list[[clock]]
      
        # Extract the first row of coefficients
        coeff <- summary(model)$coefficients[1, ]
	      n_event <- summary(model)$nevent
	      n_total <- summary(model)$n
        mort_coefficients[[clock]] <- c(coeff, n_event, n_total)
    }
    
    # Combine the coefficients into a data frame for the current disease
    results <- do.call(rbind, mort_coefficients)  

  return(as.data.frame(results))
}

out <- extract_coefficients(res.cox)
names(out)<- c("logHR","HR", "SE","Z","P","N_cases","N_total")
out$Clock <- rownames(out)

out$LCI <- exp(out$logHR - 1.96*out$SE)
out$UCI <- exp(out$logHR + 1.96*out$SE)

out1 <- out[,c("Clock","N_cases","N_total","HR","LCI","UCI","P")]

out2 <- merge(out1, local_p, by=c("Clock")) 
out3 <- merge(out2, global_p, by=c("Clock")) 

write.csv(out3, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/MortalityHR_16May2025.csv", row.names=F)


### Logistic regression/AUC analysis ###

library(pROC)


auc_est <- list()

# Regression model
for(clock in clocks){
      model_formula <- as.formula(paste0("dead ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex"))
      m1 <- glm(model_formula, data = d)
      probs <- predict(m1, type="response")
      roc_obj <- roc(d$dead, probs)
      auc_est[[clock]] <- auc(roc_obj) 
}

auc_all = as.data.frame(do.call('rbind', auc_est))

m1 <- glm(dead ~ age + bmi + years + pack_years + units + rank + sex, data = d)
probs <- predict(m1, type="response")
roc_obj <- roc(d$dead, probs)
auc_null <- auc(roc_obj) 

auc_all$null <- auc_null

auc_all$Clock <- rownames(auc_all)
names(auc_all)[1] <- "full"

auc_all$delta <- auc_all$full - auc_all$null

auc_all <- auc_all[,c(3,1,2,4)]

write.csv(auc_all, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Mortality_AUC_16May2025.csv", row.names=F)