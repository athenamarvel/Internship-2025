#Analyse des vitesses de nage pour 3 traitements distincts entre 2 périodes de temps

#### 1. File directory ----------------------

setwd("/Users/athena.marvel/Documents/Stage M1 2025/R data")

#### 2. Load Packages ####

library(readxl)       # Load Excel tables
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)      # for reordering factors
library(lme4)         # linear mixed effects model
library(lmerTest)     # adds p-value to fixed effects
library(influence.ME) # check for influential data points in models
library(MuMIn)        # calculate R² (marginal and conditional)
library(car)          # model multicollinearity (VIF)
library(lmtest)       # lrtest() to compare models and the significance of interaction terms
library(performance)  # model diagnostics
library(see)          # model diagnostics
library(emmeans)      # explore effects / simple slopes / contrasts
library(ggeffects)    # for confidence intervals of model predictions
library(lspline)      # smoothing using splines
library(DHARMa)       # for simulated residuals function and for testing dispersion of the model
library(nlme)         # for non linear models


#### 3. Functions ####

# Function to compare models based on AIC
compare_aic <- function(models) {
  aic_results <- sapply(models, AIC)
  best_model <- which.min(aic_results)
  return(list(best_model = models[[best_model]], AIC = aic_results[best_model]))
}


#### 4. Data import ####

# 4.1. Data loading

dfC8 <- read.csv("data_speed_AHL_DMSO.csv", header = TRUE, na.strings = "NA")

dfmenthe <- read.csv("data_speed_menthe.csv", header = TRUE, na.strings = "NA")

df <- rbind(dfC8, dfmenthe)

# Pour distinguer les fish par traitement 
df$fish_ID <- paste0(df$fish, "_", df$treatment)
df$fish_ID <- as.factor(df$fish_ID)
df$treatment <- as.factor(df$treatment)
head(df)

# 4.3. Variables numeric/factor
df$part <- factor(df$part, levels = c("Before", "After"))

# 4.4. Centrer la vriable minute (pour garder minute et part dans le modèle)
df$minute_centered <- ave(df$minute, df$part, FUN = function(x) x - mean(x))


#### 5. Modeling ####
model_list <- list(
  model1 <- lmer(m1 ~ part * treatment * minute_centered + (1 | fish_ID), data = df),
  model2 <- lmer(m1 ~ part * treatment + minute_centered + (1 | fish_ID), data = df),
  model3 <- lmer(m1 ~ part + treatment * minute_centered + (1 | fish_ID), data = df),
  model4 <- lmer(m1 ~ part * minute_centered + treatment + (1 | fish_ID), data = df),
  model5 <- lmer(m1 ~ part + minute_centered + treatment + (1 | fish_ID), data = df)
)

# Run the comparison
best_model_aic <- compare_aic(model_list)
print(best_model_aic)


model <- model1
isSingular(model1, tol = 1e-4)

# Model checkings
# 1. Model Summary
summary(model) # extraire les p-values de là
anova(model)


# 2. Check residuals for normality & homoscedasticity
plot(model)  # Base method
qqnorm(resid(model))
qqline(resid(model))
hist(resid(model), breaks = 30)

# 4. Calculate R² (marginal and conditional)
r.squaredGLMM(model)

# 5. Check multicollinearity (VIF)
vif(model)

# 6. Check assumptions
check_model(model)

# 7. Check fixed effects
emm <- emmeans(model, ~ part * treatment * minute_centered)
pairs(emm, simple = "treatment")  # compare DMSO vs AHL pendant période
pairs(emm, simple = "part")      # compare before vs after pour chaque côté
pairs(emm, simple = "minute_centered")      # compare before vs after pour chaque côté
pairs(emm, simple = "m1") 



#### 6. Old Plots ####
emm_df <- as.data.frame(emm)
head(emm_df)

# Barplot : axe x = part, couleur = treatment, y = temps
ggplot(emm_df, aes(x = part, y = emmean,
                   fill = treatment, group = treatment)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                position = position_dodge(0.8), width = 0.2) +
  theme_minimal()

#labs(title = "Temps estimé par période et par côté",
#x = "Période (part)", y = "Temps estimé (s)",
#fill = "Côté de l’arène")

# Plot individuels par poisson
ggplot(df, aes(x = part, y = m1,
               color = treatment, group = interaction(fishID, treatment))) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ fishID, nrow = 3) +
  labs(title = "Évolution de la vitesse (par poisson)",
       x = "Période (part)", y = "m1",
       color = "Traitement") +
  theme_minimal()


# Summary plot
ggplot(emm_df, aes(x = part, y = emmean,
                   color = treatment, group = treatment)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.1) +
  labs(title = "Effets modélisés : vitesse par période et côté",
       x = "Période (part)", y = "vitesse",
       color = "Côté de l’arène") +
  theme_minimal() +
  theme(text = element_text(size = 13))



#### 7. New Plot ####

# graphe avec tous les effets
# Génère les prédictions à partir du modèle
pred <- ggpredict(model, terms = c("minute_centered", "treatment", "part"))

# Transforme en data.frame
pred_df <- as.data.frame(pred)

# Renomme les colonnes pour clarté
pred_df <- pred_df %>%
  rename(
    minute_centered = x,
    predicted = predicted,
    treatment = group,
    part = facet
  )
head(pred_df)

# Reorder des facteurs
pred_df$part <- factor(pred_df$part, levels = c("Before", "After"))
pred_df$treatment <- factor(pred_df$treatment, levels = c("C8AHL", "C8DMS0","menthe"))
head(pred_df)

# Ajoute la moyenne des minutes pour repasser sur l'échelle "minute" réelle
mean_minutes_df <- data.frame(
  part = c("Before", "After"),
  mean_minute = c(5.5, 15.5)
)

# Fusion + recalcul de la minute réelle
pred_df <- pred_df %>%
  left_join(mean_minutes_df, by = "part") %>%
  mutate(minute = minute_centered + mean_minute)
head(pred_df)

pred_df <- pred_df %>%
  mutate(minute = ceiling(minute))
head(pred_df)

# Couleurs manuelles
color_palette <- c(
  "C8AHL" = "#879BFF",   # bleu-violet moyen (ligne/point)
  "menthe" = "#27ae60",
  "C8DMS0" = "black"
)

fill_palette <- c(
  "C8AHL" = "#c4c9ff",   # bleu-violet très clair (fill/ribbon)
  "menthe" = "#a9dfbf",
  "C8DMS0" = "#DCDCDC"
)

border_palette <- c(
  "C8AHL" = "#5f6fdb",   # bleu-violet foncé (bordure)
  "menthe" = "#229954",
  "C8DMS0" = "#a9a9a9"
)

# Graphique final : effet modélisé du temps selon la minute, par part et côté
ggplot() +
  geom_jitter(data = df,
              aes(x = minute, y = m1, fill = treatment, color = treatment),
              shape = 21, alpha = 0.6, stroke = 1, size = 2,
              width = 0.05, height = 0.01) +
  
  geom_line(data = pred_df,
            aes(x = minute, y = predicted,
                color = treatment, group = interaction(treatment, part)),
            linewidth = 1.2) +
  
  geom_ribbon(data = pred_df,
              aes(x = minute, ymin = conf.low, ymax = conf.high,
                  fill = treatment, group = interaction(treatment, part)),
              alpha = 0.2, color = NA) +
  
  geom_line(data = pred_df, aes(x = minute, y = conf.high, color = treatment, group = interaction(treatment, part)),
            linewidth = 0.4, linetype = "solid", alpha = 0.8) +
  
  geom_line(data = pred_df, aes(x = minute, y = conf.low, color = treatment, group = interaction(treatment, part)),
            linewidth = 0.4, linetype = "solid", alpha = 0.8) +
  
  labs(title = "Effets modélisés de la vitesse (par minute, période, condition)",
       x = "Temps (min)", y = "vitesse de nage des poissons (mm/s)", color = "Condition", fill = "Condition") +
  
  scale_color_manual(values = color_palette) +  # uniquement color_palette
  scale_fill_manual(values = fill_palette) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
