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
df_AHL <- read.csv("data_AHL_model.csv", header = TRUE, na.strings = "NA")
head(df_AHL)

# 4.2. Data operations
df_AHL_long <- df_AHL %>%
  mutate(nocue = total_sum - sum) %>%               # 1. Créer la variable pour nocue et renommer la colonne directement "nocue"
  pivot_longer(
    cols = c(sum, nocue), # On pivote "sum" (= cue) et "nocue"
    names_to = "compartment",                            # 2. Indique la nouvelle colonne pour le type de compartiment
    values_to = "time"                                   # 3. Nouvelle colonne pour le temps associé
  ) %>%
  mutate(
    compartment = ifelse(as.character(compartment) == "sum", "cue", "nocue")  # 4. Renomme les types de compartiment
  ) %>%
  select(fish, minute, compartment, part, time, total_sum, ratioside)  # 5. Réorganise les colonnes

df_AHL_long <- df_AHL_long %>%
  mutate(ratioside = 100 * time / total_sum)

df_AHL <- df_AHL_long

df_AHL$fishSide <- df_AHL$compartment
df_AHL$fishID <- as.factor(df_AHL$fish)

# 4.3. Variables numeric/factor
df_AHL$part <- factor(df_AHL$part, levels = c("Before", "After"))
df_AHL$fishSide <- factor(df_AHL$fishSide, levels = c("nocue", "cue"))

# 4.4. Centrer la vriable minute (pour garder minute et part dans le modèle)
df_AHL$minute_centered <- ave(df_AHL$minute, df_AHL$part, FUN = function(x) x - mean(x))


#### 5. Modeling ####
model_list <- list(
  model1 <- lmer(time ~ part * fishSide * minute_centered + (1 | fishID), data = df_AHL),
  model2 <- lmer(time ~ part * fishSide + minute_centered + (1 | fishID), data = df_AHL),
  model3 <- lmer(time ~ part + fishSide * minute_centered + (1 | fishID), data = df_AHL),
  model4 <- lmer(time ~ part * minute_centered + fishSide + (1 | fishID), data = df_AHL),
  model5 <- lmer(time ~ part + minute_centered + fishSide + (1 | fishID), data = df_AHL)
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

library(emmeans)

emmeans(model, ~ part * fishSide * minute_centered,
        at = list(minute_centered = c(min(df_AHL$minute_centered),
                                      max(df_AHL$minute_centered))))

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
emm <- emmeans(model, ~ part * fishSide * minute_centered)
pairs(emm, simple = "fishSide")  # compare cue vs nocue dans chaque période
pairs(emm, simple = "part")      # compare before vs after pour chaque côté
pairs(emm, simple = "minute_centered")      # compare before vs after pour chaque côté

#### 6. Old Plots ####
emm_df_AHL <- as.data.frame(emm)
head(emm_df_AHL)

# Barplot : axe x = part, couleur = fishSide, y = temps
ggplot(emm_df_AHL, aes(x = part, y = emmean,
                       fill = fishSide, group = fishSide)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                position = position_dodge(0.8), width = 0.2) +
  labs(title = "Temps estimé par période et par côté",
       x = "Période (part)", y = "Temps estimé (s)",
       fill = "Côté de l’arène") +
  theme_minimal()

# Plot individuels par poisson
ggplot(df_AHL, aes(x = part, y = time,
                   color = fishSide, group = interaction(fishID, fishSide))) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ fishID, nrow = 3) +
  labs(title = "Évolution du temps passé (par poisson)",
       x = "Période (part)", y = "Temps (s)",
       color = "Côté") +
  theme_minimal()


# Summary plot
ggplot(emm_df_AHL, aes(x = part, y = emmean,
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

#### 7. New Plot ####

library(ggeffects)
library(tibble)
library(dplyr)
library(ggplot2)

# Génère les prédictions à partir du modèle
pred <- ggpredict(model, terms = c("minute_centered", "fishSide", "part"))

# Transforme en data.frame
pred_df_AHL <- as.data.frame(pred)

# Renomme les colonnes pour clarté
pred_df_AHL <- pred_df_AHL %>%
  rename(
    minute_centered = x,
    predicted = predicted,
    fishSide = group,
    part = facet
  )

# Reorder des facteurs
pred_df_AHL$part <- factor(pred_df_AHL$part, levels = c("Before", "After"))
pred_df_AHL$fishSide <- factor(pred_df_AHL$fishSide, levels = c("nocue", "cue"))

# Ajoute la moyenne des minutes pour repasser sur l'échelle "minute" réelle
mean_minutes_df_AHL <- data.frame(
  part = c("Before", "After"),
  mean_minute = c(5.5, 15.5)
)

# Fusion + recalcul de la minute réelle
pred_df_AHL <- pred_df_AHL %>%
  left_join(mean_minutes_df_AHL, by = "part") %>%
  mutate(minute = minute_centered + mean_minute)

# Ajouter une colonne pour % de temps dans pred_df_AHL et df_AHL
df_AHL <- df_AHL %>% mutate(time_pct = 100 * time / 60)
pred_df_AHL <- pred_df_AHL %>%
  mutate(predicted_pct = 100 * predicted / 60,
         conf.low_pct = 100 * conf.low / 60,
         conf.high_pct = 100 * conf.high / 60)

# Recode des noms pour l'affichage
df_AHL <- df_AHL %>%
  mutate(compartment = dplyr::recode(as.character(compartment),
                                     "cue" = "C8AHL", "nocue" = "DMSO"))
pred_df_AHL <- pred_df_AHL %>%
  mutate(fishSide = dplyr::recode(as.character(fishSide),
                                  "cue" = "C8AHL", "nocue" = "DMSO"))

# Couleurs manuelles (utilisées pour points, ribbon)
# Nouvelle palette saumon pour C8AHL
color_palette <- c("C8AHL" = "#FA8072",   # ligne saumon clair
                   "DMSO" = "grey50")    # gris inchangé

fill_palette <- c("C8AHL" = "#FAD6C4",    # fill saumon pâle
                  "DMSO" = "#DCDCDC")    # fill DMSO inchangé

border_palette <- c("C8AHL" = "#E07B5F",  # bordure saumon foncée
                    "DMSO" = "#a9a9a9")  # bordure DMSO inchangée

# Graphique final : effet modélisé du temps selon la minute, par part et côté
ggplot() +
  # Points avec contours
  geom_jitter(data = df_AHL,
              aes(x = minute, y = time_pct, fill = compartment, color = compartment),
              shape = 21, alpha = 0.6, stroke = 1, size = 2,
              width = 0.05, height = 0.01) +
  
  # Lignes prédictions C8AHL 0-10 min pointillées rouges
  geom_line(data = pred_df_AHL %>% filter(fishSide == "C8AHL", minute >= 0, minute <= 10),
            aes(x = minute, y = predicted_pct, color = fishSide, group = interaction(fishSide, part)),
            linewidth = 1.2, linetype = "dashed") +
  
  # Lignes prédictions C8AHL 11-20 min pleines rouges
  geom_line(data = pred_df_AHL %>% filter(fishSide == "C8AHL", minute > 10, minute <= 20),
            aes(x = minute, y = predicted_pct, color = fishSide, group = interaction(fishSide, part)),
            linewidth = 1.2, linetype = "solid") +
  
  # Lignes prédictions DMSO 0-10 min pointillées grises
  geom_line(data = pred_df_AHL %>% filter(fishSide == "DMSO", minute >= 0, minute <= 10),
            aes(x = minute, y = predicted_pct, group = interaction(fishSide, part)),
            linewidth = 1.2, linetype = "dashed", color = "grey50") +
  
  # Lignes prédictions DMSO 11-20 min pleines grises
  geom_line(data = pred_df_AHL %>% filter(fishSide == "DMSO", minute > 10, minute <= 20),
            aes(x = minute, y = predicted_pct, group = interaction(fishSide, part)),
            linewidth = 1.2, linetype = "solid", color = "grey50") +
  
  # Rubans de confiance (inchangés)
  geom_ribbon(data = pred_df_AHL,
              aes(x = minute, ymin = conf.low_pct, ymax = conf.high_pct,
                  fill = fishSide, color = fishSide, group = interaction(fishSide, part)),
              alpha = 0.2, color = NA) +
  
  # Bordures hautes des ribbons
  geom_line(data = pred_df_AHL, aes(x = minute, y = conf.high_pct, color = fishSide, group = interaction(fishSide, part)),
            linewidth = 0.4, linetype = "solid", alpha = 0.8) +
  
  # Bordures basses des ribbons
  geom_line(data = pred_df_AHL, aes(x = minute, y = conf.low_pct, color = fishSide, group = interaction(fishSide, part)),
            linewidth = 0.4, linetype = "solid", alpha = 0.8) +
  
  labs(title = "Effets modélisés du temps passé (par minute, côté et période)",
       x = "Temps (min)", y = "Temps passé dans la condition (%)", color = "Condition", fill = "Condition") +
  
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = fill_palette) +
  
  scale_y_continuous(
    limits = c(-1, 101),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0", "25", "50", "75", "100")
  ) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

