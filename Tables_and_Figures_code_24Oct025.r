################
### Table S1 ###
################

data <- read.csv("disease_meta_data_21Aug2025.csv", header=T)

dis <- read.delim("combined_disease_names_and_groups.txt")
data <- merge(data, dis, by.x="disease", by.y="name")
names(data)[c(1,14)] <- c("Dis_short","Disease")

s1 <- data[,c(15,14,2:13)]


################
### Table S2 ###
################

### merge results tables for HR and AUC into single file ###
data <- read.csv("Disease_HR_AUC_merged_24Oct2025.csv", header=T)

dis <- read.delim("combined_disease_names_and_groups.txt")
data <- merge(data, dis, by.x="disease", by.y="name")
names(data)[c(1,23)] <- c("Dis_short","Disease")

data$Clock[data$Clock == "DNAmGrimAge"] <- "GrimAgeV1"
data$Clock[data$Clock == "DNAmGrimAge.1"] <- "GrimAgeV2"
data$Clock[data$Clock == "Horvathv1"] <- "Horvath"
data$Clock[data$Clock == "Horvathv2"] <- "Horvath_SkinBlood"

s2 <- data[,c(24, 23, 2:16, 17:22)]


################
### Table S3 ###
################

data <- read.csv("Mortality_HR_AUC_24Oct2025.csv", header=T)

data$Clock[data$Clock == "DNAmGrimAge"] <- "GrimAgeV1"
data$Clock[data$Clock == "DNAmGrimAge.1"] <- "GrimAgeV2"
data$Clock[data$Clock == "Horvathv1"] <- "Horvath"
data$Clock[data$Clock == "Horvathv2"] <- "Horvath_SkinBlood"

s3 <- data


#####################
### Table S4 & S5 ###
#####################

tmp <- s2

tmp$logHR <- log(tmp$HR)

library(dplyr)

result <- tmp %>%
  group_by(Clock) %>%
  summarise(
    mean_value = mean(logHR, na.rm = TRUE),
    sd_value = sd(logHR, na.rm = TRUE)
  )

s4 <- as.data.frame(result)


library(lmerTest)
tmp$Clock <- as.factor(tmp$Clock)
tmp$Clock <- relevel(tmp$Clock, ref="GrimAgeV1")
tmp$logHR[tmp$Clock=="DNAmTL"] <- tmp$logHR[tmp$Clock=="DNAmTL"] * -1
model <- lmer(logHR ~ Clock + (1|Disease), data=tmp)

s5_tmp <- as.data.frame(summary(model)$coefficients)

s5_tmp$Clock <- rownames(s5_tmp)

names(s5_tmp)[1:5] <- c("Estimate","SE","df","t","P") 

s5 <- s5_tmp[,c(6,1:5)]


################
### Table S6 ###
################

mort <- s3
mort$Group <- mort$Disease <- "Mortality"
mort <- mort[,c(21,20,1:8)]
mort$abslogHR <- abs(log(mort$HR))
mort <- mort[mort$P<0.05/174,]

d2 <- s2[s2$P<0.05/174,c(1:8,10,9)]
d2$abslogHR <- abs(log(d2$HR))

d3 <- rbind(d2, mort)

clocks <- unique(d3$Clock)
tmp <- data.frame(d3[0,])

for(i in 1:length(clocks)){
clock <- clocks[i]
d4 <- d3[d3$Clock == clock,]
d5 <- d4[which(d4$abslogHR >= d4$abslogHR[d4$Disease=="Mortality"]),]
ifelse(nrow(d5) >0, tmp <- rbind(tmp, d5), tmp <- tmp)
}

s6 <- tmp[order(tmp$Clock, -tmp$abslogHR),]

length(unique(s6$Disease)) - 1


################
### Table S7 ###
################

s7 <- s2[(s2$P<0.05/174 & s2$auc_diff>0.01 & s2$P_auc_diff<0.05),]


################
### Table S8 ###
################

data <- read.csv("All_Results_23Oct2025.csv", header=T)

dis <- read.delim("combined_disease_names_and_groups.txt")
data <- merge(data, dis, by.x="disease", by.y="name")
names(data)[c(1,18)] <- c("Dis_short","Disease")

data$Clock[data$Clock == "DNAmGrimAge"] <- "GrimAgeV1"
data$Clock[data$Clock == "DNAmGrimAge.1"] <- "GrimAgeV2"
data$Clock[data$Clock == "Horvathv1"] <- "Horvath"
data$Clock[data$Clock == "Horvathv2"] <- "Horvath_SkinBlood"

data <- data[,c(19,18,2:17)]

data <- data[order(data$Group, data$Disease, data$Clock),]


int <- read.csv("Interaction_Results_23Oct2025.csv", header=T)

int2 <- merge(int, dis, by.x="disease", by.y="name")
names(int2)[c(1,12)] <- c("Dis_short","Disease")
int2$Clock[int2$Clock == "DNAmGrimAge"] <- "GrimAgeV1"
int2$Clock[int2$Clock == "DNAmGrimAge.1"] <- "GrimAgeV2"
int2$Clock[int2$Clock == "Horvathv1"] <- "Horvath"
int2$Clock[int2$Clock == "Horvathv2"] <- "Horvath_SkinBlood"

int2$global_p <- int2$df_global <- int2$chisq_global <- int2$local_p <- int2$df_local <-  int2$chisq_local <- NA

int3 <- int2[c(13,12,2:10,14:19,11)]

data1 <- rbind(data, int3)


data1$Covariates <- ifelse(
  data1$model %in% c("full_smoke", "full_nosmoke", "full", "full_female", "full_male", "smoke_int", "sex_int"),
  "Full adjustment",
  "Basic adjustment"
)


data1$Strata <- ifelse(
  data1$model %in% c("full_male", "basic_male"),
  "Male", ifelse(
  data1$model %in% c("full_female", "basic_female"),
  "Female", ifelse(
  data1$model %in% c("full_smoke", "basic_smoke"),
  "Ever Smoker", ifelse(
  data1$model %in% c("full_nosmoke", "basic_nosmoke"),
  "Never Smoker", ifelse(
  data1$model %in% c("sex_int", "basic_sex_int"),
  "Interaction (Sex - Male reference)", ifelse(
  data1$model %in% c("smoke_int", "basic_smoke_int"),
  "Interaction (Smoking - Never smoker reference)",
  NA
))))))

       
s8 <- data1[,-18]

s8 <- s8[order(s8$Group, s8$Disease, s8$Clock),]

### analyses for Table 1 ###

bonf <- data1[data1$P<0.05/174 & (data1$model=="basic" | data1$model=="full"),]
nom <- data1[data1$P<0.05 & (data1$model=="basic" | data1$model=="full"),]

t_bonf <- table(bonf$model, bonf$Clock)
t_nom <- table(nom$model, nom$Clock)

out <- as.data.frame(rbind(t_nom, t_bonf))

out <- as.data.frame(t(out[,c("Hannum","Horvath","Lin","Zhang_10","PhenoAge","Horvath_SkinBlood","GrimAgeV1","DNAmTL","DunedinPoAm38","DunedinPACE","GrimAgeV2","YingCausAge","YingDamAge", "YingAdaptAge")]))

write.csv(out, file="T1_summary.csv", row.names=F)


sexint <- data1[data1$P<0.05/174 & data1$model=="sex_int",]
smokint <- data1[data1$P<0.05/174 & data1$model=="smoke_int",]






######################################
### write out to an excel workbook ###
######################################

install.packages("openxlsx")
library(openxlsx)

# Create workbook
wb <- createWorkbook()

# List of data frames
data_list <- list(s1 = s1, s2 = s2, s3 = s3, s4 = s4, s5 = s5, s6 = s6, s7 = s7, s8 = s8)

# Write each data frame to its own sheet
for (name in names(data_list)) {
  addWorksheet(wb, name)
  writeData(wb, name, data_list[[name]])
}

# Save the workbook
saveWorkbook(wb, "multiple_sheets.xlsx", overwrite = TRUE)












################
### Figure 1 ###
################


library(lme4)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)


tmp <- s2
tmp$logHR <- log(tmp$HR)

## positively code the DNAmTL effect sizes ##
tmp$logHR[tmp$Clock=="DNAmTL"] <- -1* tmp$logHR[tmp$Clock=="DNAmTL"]


# 1. Fit your mixed model
model <- lmer(logHR ~ Clock + (1|Disease), data = tmp)

# 2. Get estimated marginal means for Clock
emm <- emmeans(model, ~ Clock)

# 3. Get all pairwise comparisons
pairs_res <- pairs(emm)

# 4. Convert to data frame and separate contrasts
pairs_df <- as.data.frame(pairs_res) %>%
  separate(contrast, into = c("Clock1", "Clock2"), sep = " - ")

# 5. Create full grid of all clock pairs (for plotting)
clocks <- unique(c(pairs_df$Clock1, pairs_df$Clock2))
heatmap_df <- expand.grid(Clock1 = clocks, Clock2 = clocks, stringsAsFactors = FALSE)

# 6. Join pairwise results with the full grid
heatmap_df <- heatmap_df %>%
  left_join(pairs_df %>% select(Clock1, Clock2, estimate, p.value), by = c("Clock1", "Clock2"))

# 7. Fill diagonal with 0 differences and NA p-values
heatmap_df <- heatmap_df %>%
  mutate(
    estimate = ifelse(Clock1 == Clock2, 0, estimate),
    p.value = ifelse(Clock1 == Clock2, NA, p.value)
  )

# 8. Keep only the upper triangle (alphabetical order)
heatmap_df_upper <- heatmap_df %>%
  filter(Clock1 < Clock2)


heatmap_df_upper <- heatmap_df_upper %>%
  mutate(signif = ifelse(!is.na(p.value) & p.value < 0.05,
                         formatC(p.value, format = "e", digits = 1),  # 1 digit after decimal + exponent counts as 2 digits
                         ""))

write.csv(heatmap_df_upper, file="Fig1B_data.csv", row.names=F)

p2 <- ggplot(heatmap_df_upper, aes(x = Clock1, y = Clock2, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), na.rm = TRUE, size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Mean difference") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(title = "Heatmap of pairwise clock differences in average log hazards")


library(dplyr)
library(ggplot2)


d2 <- tmp

d2$gen <- ifelse(d2$Clock %in% c("DunedinPACE", "DunedinPoAm38"), "Third",
           ifelse(d2$Clock %in% c("GrimAgeV1", "GrimAgeV2", "PhenoAge", "Zhang_10"), "Second",
           ifelse(d2$Clock == "DNAmTL", "Telomere length",
           "First")))


# Reorder Clock by median logHR (descending)
d2 <- d2 %>%
  mutate(Clock = reorder(Clock, logHR, FUN = median, decreasing = TRUE))

write.csv(d2, file="Fig1A_data.csv", row.names=F)

p1 <- ggplot(d2, aes(x = Clock, y = logHR, fill = gen)) +
  geom_violin(trim = FALSE, color = "black") +
  stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +  # Median
  stat_summary(fun.data = function(x) {
    data.frame(y = median(x),
               ymin = quantile(x, 0.25),
               ymax = quantile(x, 0.75))
  }, geom = "errorbar", width = 0.2, color = "black") +  # IQR
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Null effect line
  labs(title = "",
       x = "Epigenetic Clock",
       y = "log(Hazard Ratio)",
       fill = "Generation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14)  # Only applies if using facets
  )


combined_plot <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 1)


ggsave("Figure1_24Oct2025.pdf", combined_plot, width = 20, height = 16)




################
### Figure 2 ###
################


df <- s2[s2$P<0.05/174 & s2$P_auc_diff<0.05 & s2$auc_diff >0.01,] 

# Create a combined label for disease-clock pairs
df <- df %>%
  mutate(disease_clock = paste(Disease, Clock, sep = " - "),
         disease_clock = factor(disease_clock, levels = unique(disease_clock)))


library(tidyr)
df_long <- df %>%
  select(disease_clock, auc_null, auc_full) %>%
  pivot_longer(cols = c(auc_null, auc_full),
               names_to = "model",
               values_to = "AUC") %>%
  mutate(
    model = recode(model,
                   "auc_null" = "Null Model",
                   "auc_full" = "Full Model")
  )


df_long <- df_long %>%
  mutate(model = factor(model, levels = c("Null Model", "Full Model")))

write.csv(df_long, file="Fig2_data.csv", row.names=F)

# Plot
fig <- ggplot(df_long, aes(x = AUC, y = disease_clock, color = model, fill = model)) +
  geom_segment(data = df,
               aes(x = auc_null, xend = auc_full, y = disease_clock, yend = disease_clock),
               inherit.aes = FALSE, color = "grey80", size = 2) +
  geom_point(shape = 21, size = 4, stroke = 1) +
  scale_color_manual(values = c("Null Model" = "darkblue", "Full Model" = "lightblue")) +
  scale_fill_manual(values = c("Null Model" = "darkblue", "Full Model" = "lightblue")) +
  labs(x = "AUC", y = NULL, color = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = c(0.97, 0.97),  # top-right inside plot
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.box.background = element_rect(color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 12)
  )

ggsave("Figure2_24Oct2025.pdf", fig, , width = 12, height = 10)