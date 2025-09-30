### aveage logHR by clock analysis ### 

library("readxl")
d <- read_excel("Supplementary_Tables_25June2025.xlsx", sheet = 3)
d1 <- as.data.frame(d)  
d1$logHR <- log(d1$HR)

d2 <- d1[-(which(d1$Disease == "Mortality")),]

library(dplyr)

result <- d2 %>%
  group_by(Clock) %>%
  summarise(
    mean_value = mean(logHR, na.rm = TRUE),
    sd_value = sd(logHR, na.rm = TRUE)
  )

out <- as.data.frame(result)


library(lmerTest)
d2$Clock <- as.factor(d2$Clock)
d2$Clock <- relevel(d2$Clock, ref="GrimAgeV1")
model <- lmer(logHR ~ Clock + (1|Disease), data=d2)
summary(model)

# Get estimated marginal means for Clock
emm <- emmeans(model, ~ Clock)

# Pairwise comparisons between all Clock levels (paired within Disease)
pairs(emm, adjust = "tukey")




library(lme4)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Fit your mixed model
model <- lmer(logHR ~ Clock + (1|Disease), data = d2)

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
  mutate(signif = case_when(
    !is.na(p.value) & p.value < 0.001 ~ "***",
    !is.na(p.value) & p.value < 0.01  ~ "**",
    !is.na(p.value) & p.value < 0.05  ~ "*",
    TRUE                             ~ ""
  ))

ggplot(heatmap_df_upper, aes(x = Clock1, y = Clock2, fill = estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), na.rm = TRUE, size = 6) +  # na.rm=TRUE avoids warning
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Mean difference") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  labs(title = "Upper Triangular Heatmap of Pairwise Clock Differences (logHR)")

