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
local_p[[clock]] <- res.cox.zph.local[[clock]]
}

local_p = as.data.frame(do.call('rbind', local_p))
local_p$Clock <- as.character(rownames(local_p))
names(local_p)[1:3] <- c("chisq_local", "df_local", "local_p")


## prop hazards (global) ##
global_p = list()

for(clock in names(res.cox.zph.global)){
global_p[[clock]] <- res.cox.zph.global[[clock]]
}

global_p = as.data.frame(do.call('rbind', global_p))
global_p$Clock <- as.character(rownames(global_p))
names(global_p)[1:3] <- c("chisq_global", "df_global", "global_p")


## Extract summary output from Cox models for plotting 
extract_coefficients <- function(models_list) {
  # Initialize an empty list to store results
  results <- list()
  mort_coefficients <- list()
   
    # Iterate over each epigenetic clock model
    for (clock in names(models_list)) {
      model <- models_list[[clock]]
      
        # Extract the first row of coefficients
        coeff <- summary(model)$coefficients[1,]
	      n_event <- summary(model)$nevent
	      n_total <- summary(model)$n
        mort_coefficients[[clock]] <- c(coeff, n_event, n_total)
    }
    
    # Combine the coefficients into a data frame for the current disease
    results <- do.call(rbind, mort_coefficients)  

  return(as.data.frame(results))
}

out <- extract_coefficients(res.cox)
names(out) <- c("logHR","HR", "SE","Z","P","N_cases","N_total")
out$Clock <- rownames(out)

out$LCI <- exp(out$logHR - 1.96*out$SE)
out$UCI <- exp(out$logHR + 1.96*out$SE)

out1 <- out[,c("Clock","N_cases","N_total","HR","LCI","UCI","Z","P")]

out2 <- merge(out1, local_p, by=c("Clock")) 
out3 <- merge(out2, global_p, by=c("Clock")) 



### Logistic regression/AUC analysis ###

library(pROC)


auc_est <- list()
auc_null <- list()
auc_p_comp <- list()
auc_z <- list()

# Regression model
for(clock in clocks){
      model_formula <- as.formula(paste0("dead ~ scale(", clock, ") + age + bmi + years + pack_years + units + rank + sex"))
      m1 <- glm(model_formula, data = d)
      probs <- predict(m1, type="response")
      roc_obj <- roc(d$dead, probs)
      auc_est[[clock]] <- auc(roc_obj) 


	  model_formula2 <- as.formula(paste0("dead ~ age + bmi + years + pack_years + units + rank + sex"))
      m2 <- glm(model_formula2, data = d)
      probs2 <- predict(m2, type="response")
      roc_obj2 <- roc(d$dead, probs2)
      auc_null[[clock]] <- auc(roc_obj2)  
	  
	  auc_p_comp[[clock]] <- roc.test(roc_obj, roc_obj2)$p.value
	  auc_z[[clock]] <- roc.test(roc_obj, roc_obj2)$statistic
}

auc_all = as.data.frame(do.call('rbind', auc_est))
auc_null <- as.data.frame(do.call('rbind', auc_null))
auc_p_comp <- as.data.frame(do.call('rbind', auc_p_comp))
auc_z_comp <- as.data.frame(do.call('rbind', auc_z))

out4 <- cbind(auc_null, auc_all, auc_z_comp, auc_p_comp)
names(out4) <- c("Null_AUC","Full_AUC","AUC_Z_comp","AUC_P_comp")
out4$Clock <- row.names(out4)

out4$AUC_diff <- out4$Full_AUC - out4$Null_AUC
out4 <- out4[,c("Clock","Null_AUC","Full_AUC","AUC_diff","AUC_Z_comp","AUC_P_comp")]

out5 <- merge(out3, out4, by="Clock")

write.csv(out5, "/Cluster_Filespace/Marioni_Group/Riccardo/Christos_17Oct2024/Results/Mortality_HR_AUC_24Oct2025.csv", row.names=F)