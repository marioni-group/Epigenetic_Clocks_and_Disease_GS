### aveage logHR by clock analysis ### 
library("readxl")
library(dplyr)
library(lmerTest)

d <- read_excel("Supplementary_Tables_7July2025.xlsx", sheet = 2)
d1 <- as.data.frame(d)  
d1$logHR <- log(d1$HR)

d2 <- d1[-(which(d1$Disease == "Mortality")),]

result <- d2 %>%
  group_by(Clock) %>%
  summarise(
    mean_value = mean(logHR, na.rm = TRUE),
    sd_value = sd(logHR, na.rm = TRUE)
  )

out <- as.data.frame(result)

d2$Clock <- as.factor(d2$Clock)
d2$Clock <- relevel(d2$Clock, ref="GrimAgeV1")
summary(lmer(logHR ~ Clock + (1|Disease), data=d2))