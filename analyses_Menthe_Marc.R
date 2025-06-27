#### 1. File directory ----------------------

setwd("/Users/athena.marvel/Documents/Stage M1 2025/R data")


#### 2. Load Packages ####

library(readxl)       # Load Excel tables
library(ggplot2)
library(dplyr)
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
df <- read_xlsx(path="data_Menthe_Marc.xlsx", sheet="data_Menthe", col_names=TRUE, na="NA")

# 4.2. Variables numeric/factor
df$part <- factor(df$part, levels = c("Before", "After"))
df$fishSide <- factor(df$fishSide, levels = c("nocue", "cue"))


#### 5. Modeling ####
model_list <- list(
  model1 <- lmer(time ~ part * fishSide + (1 | fishID), data = df),
  model2 <- lmer(time ~ part + fishSide + (1 | fishID), data = df)
)

model3 <- glmer(time ~ part + fishSide + (1 | fishID), data = df, family=Gamma(link='log'))

# Run the comparison
best_model_aic <- compare_aic(model_list)
print(best_model_aic)


model <- model1
isSingular(model1, tol = 1e-4) #vérifie si le modèle est singulier, si les effets aléatoires posent problème

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
emm <- emmeans(model, ~ part * fishSide)
pairs(emm, simple = "fishSide")  # compare cue vs nocue dans chaque période
pairs(emm, simple = "part")      # compare before vs after pour chaque côté

#### 6. Plots ####
emm_df <- as.data.frame(emm)

# Barplot : axe x = part, couleur = fishSide, y = temps
ggplot(emm_df, aes(x = part, y = emmean,
                   fill = fishSide, group = fishSide)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                position = position_dodge(0.8), width = 0.2) +
  labs(title = "Temps estimé par période et par côté",
       x = "Période (part)", y = "Temps estimé (s)",
       fill = "Côté de l’arène") +
  theme_minimal()


# Plot individuels par poisson
ggplot(df, aes(x = part, y = time,
               color = fishSide, group = interaction(fishID, fishSide))) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ fishID, nrow = 3) +
  labs(title = "Évolution du temps passé (par poisson)",
       x = "Période (part)", y = "Temps (s)",
       color = "Côté") +
  theme_minimal()


# Summary plot
ggplot(emm_df, aes(x = part, y = emmean,
                   color = fishSide, group = fishSide)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.1) +
  labs(title = "Effets modélisés : temps estimé par période et côté",
       x = "Période (part)", y = "Temps estimé (s)",
       color = "Côté de l’arène") +
  theme_minimal() +
  theme(text = element_text(size = 13))

