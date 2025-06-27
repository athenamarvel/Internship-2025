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

df <- read.csv("data_speed_AHL_DMSO.csv", header = TRUE, na.strings = "NA")
head(df)

df_ahl <- df %>%
  filter(treatment == "C8AHL")

df_ahl$treatment <- as.factor(df_ahl$treatment)
df_ahl$treatment <- droplevels(df_ahl$treatment)
levels(df_ahl$treatment)

df_ahl$fishID <- as.factor(df_ahl$fish)

df_ahl <- df %>% filter(treatment == "C8AHL")
df_ahl$treatment <- droplevels(as.factor(df_ahl$treatment))
df_ahl$fishID <- as.factor(df_ahl$fish)
df_ahl$part <- factor(df_ahl$part, levels = c("Before", "After"))
df_ahl$minute_centered <- ave(df_ahl$minute, df_ahl$part, FUN = function(x) x - mean(x))


#### 5. Modeling ####
model_list <- list(
  model1 = lmer(m1 ~ part * minute_centered + (1 | fishID), data = df_ahl),
  model2 = lmer(m1 ~ part + minute_centered + (1 | fishID), data = df_ahl)
)


# Run the comparison
aic_values <- sapply(model_list, AIC)
print(aic_values)

# Trouver le modèle avec le plus petit AIC
best_model_name <- names(which.min(aic_values))
best_model <- model_list[[best_model_name]]
cat("Best model based on AIC is:", best_model_name, "\n")


model <- best_model
isSingular(best_model, tol = 1e-4) # Un modèle est dit singulier lorsque la structure de ses effets aléatoires est trop complexe 

# Model checkings
# 1. Model Summary
summary(model) # extraire les p-values de là
anova(model)


# 2. Check residuals for normality & homoscedasticity
plot(model)  
qqnorm(resid(model)) # indique la présence de quelques valeurs extrêmes ou asymétrie légère = léger problème de normalité dans les résidus, surtout dans les valeurs élevées.
qqline(resid(model))
hist(resid(model), breaks = 30)

# 4. Calculate R² (marginal and conditional)
r.squaredGLMM(model)

# 5. Check multicollinearity (VIF)
vif(model)

# 6. Check assumptions
check_model(model)

# 7. Check fixed effects
emm <- emmeans(model, ~ part * minute_centered) 
pairs(emm, simple = "part")      # compare before vs after pour chaque côté
pairs(emm, simple = "minute_centered")      # compare before vs after pour chaque côté
pairs(emm, simple = "m1") 

mm <- emmeans(model, ~ part * minute_centered, 
              at = list(minute_centered = c(-4.5, 0, 4.5)))
pairs(emm, simple = "part")

emm_minute <- emmeans(model, ~ minute_centered | part, 
                      at = list(minute_centered = c(-4.5, 0, 4.5))) #compare minute_centered = -4.5 vs 0, 0 vs 4.5, etc.

pairs(emm_minute, reverse = TRUE)
emm
emm_minute


#### 6. Old Plots ####
emm_df <- as.data.frame(emm)
head(emm_df) # moyenne prédite de la variable réponse (m1) pour un certain niveau des variables fixes du modèle, en tenant compte des autres variables du modèle.
# L’emmean de 12.29 correspond à la moyenne estimée uniquement pour le groupe Before au temps centré égal à 0. Ce n’est pas la moyenne de toutes les minutes, mais celle à ce temps précis.

# Barplot : axe x = part, couleur = treatment, y = temps
emm_df_sub <- subset(emm_df, minute_centered == 0)

ggplot(emm_df_sub, aes(x = part, y = emmean, fill = part)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                position = position_dodge(0.8), width = 0.2) +
  theme_minimal() +
  labs(y = "Moyenne estimée (emmean)", x = "Part")
#labs(title = "Temps estimé par période et par côté",
#x = "Période (part)", y = "Temps estimé (s)",
#fill = "Côté de l’arène")

# Plot individuels par poisson
ggplot(df_ahl, aes(x = part, y = m1,
               color = part, group = interaction(fishID, part))) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ fishID, nrow = 3) +
  labs(title = "Évolution de la vitesse (par poisson)",
       x = "Période (part)", y = "m1",
       color = "Traitement") +
  theme_minimal()


# Summary plot
ggplot(emm_df, aes(x = part, y = emmean,
                   color = part, group = part)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.2) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.1) +
  labs(title = "Effets modélisés : vitesse de nage par période",
       x = "Période de l'expérience", y = "vitesse de nage des poissons (mm/s)",
       color = "période") +
  theme_minimal() +
  theme(text = element_text(size = 13))

#### 7. New Plot ####

library(ggeffects)
library(tibble)
library(dplyr)
library(ggplot2)

# Génère les prédictions à partir du modèle
pred <- ggpredict(model, terms = c("minute_centered", "part"))

# Transforme en data.frame
pred_df <- as.data.frame(pred)

# Renomme les colonnes pour clarté
pred_df <- pred_df %>%
  rename(
    minute_centered = x,
    predicted = predicted,
    part = group   # ici ggpredict met "group" pour part
  )

# Reorder des facteurs
pred_df$part <- factor(pred_df$part, levels = c("Before", "After"))

# Ajoute la moyenne des minutes pour repasser sur l'échelle "minute" réelle
mean_minutes_df <- data.frame(
  part = c("Before", "After"),
  mean_minute = c(5.5, 15.5)
)

# Fusion + recalcul de la minute réelle
pred_df <- pred_df %>%
  left_join(mean_minutes_df, by = "part") %>%
  mutate(minute = minute_centered + mean_minute)

pred_df <- pred_df %>%
  mutate(minute = ceiling(minute))

# Filtrer ton df original pour garder uniquement C8AHL
df_sub <- df %>% filter(treatment == "C8AHL")

# Couleurs simples par part (Before / After)
color_palette <- c("Before" = "#E75480",  # rose rouge vif (raspberry)
                    "After"  = "#2cccca")  # turquoise

fill_palette <- c("Before" = "#F8AFC1",  # rose pâle (rosé clair)
                  "After"  = "#9ee9e8")  # turquoise clair
# Graphique final : effet modélisé du temps selon la minute, par part
ggplot() +
  # Points avec contours
  geom_jitter(data = df_sub,
              aes(x = minute, y = m1, fill = part, color = part),
              shape = 21, alpha = 0.6, stroke = 1, size = 2,
              width = 0.05, height = 0.01) +
  
  # Lignes des prédictions
  geom_line(data = pred_df,
            aes(x = minute, y = predicted,
                color = part, group = part),
            linewidth = 1.2) +
  
  # Rubans de confiance
  geom_ribbon(data = pred_df,
              aes(x = minute, ymin = conf.low, ymax = conf.high,
                  fill = part, group = part),
              alpha = 0.2) +
  
  labs(title = "Effets modélisés de la vitesse (par minute et période)",
       x = "Temps (min)", y = "Temps passé dans la condition (%)",
       color = "Période", fill = "Période") +
  
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = fill_palette) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )


# graphe avec uniquement effet traitement

# Génère les prédictions à partir du modèle
pred <- ggpredict(model, terms = c("treatment"))

# Transforme en data.frame
pred_df <- as.data.frame(pred)

# Renomme les colonnes pour clarté
pred_df <- pred_df %>%
  rename(
    treatment = x,
    predicted = predicted,
    minute = group
  )

# Reorder des facteurs
pred_df$part <- factor(pred_df$part, levels = c("Before", "After"))
pred_df$treatment <- factor(pred_df$treatment, levels = c("C8AHL", "C8DMSO"))

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

# Couleurs manuelles
color_palette <- c("C8AHL" = "#2cccca", "DMSO" = "black")       # lignes, points
fill_palette <- c("C8AHL" = "#9ee9e8", "DMSO" = "#DCDCDC")      # ribbon, fill
border_palette <- c("C8AHL" = "#5acdcc", "DMSO" = "#a9a9a9")  # Couleurs intermédiaires


pred_df <- rbind(pred_df, pred_df)
pred_df$minute <- as.numeric(pred_df$minute)
pred_df[3,6] <- 20
pred_df[4,6] <- 20

# Graphique final : effet modélisé du temps selon la minute, par part et côté
ggplot() +
  # Points avec contours
  geom_jitter(data = df,
              aes(x = minute, y = m1, fill = treatment, color = treatment),
              shape = 21, alpha = 0.6, stroke = 1, size = 2,
              width = 0.05, height = 0.01) +
  
  # Lignes des prédictions
  geom_line(data = pred_df,
            aes(x = minute, y = predicted,
                color = treatment),
            linewidth = 1.2) +
  
  # Rubans de confiance
  geom_ribbon(data = pred_df,
              aes(x = minute, ymin = conf.low, ymax = conf.high,
                  fill = treatment, color = treatment),
              alpha = 0.2, color = NA) +
  
  # Bordures hautes des ribbons
  geom_line(data = pred_df, aes(x = minute, y = conf.high, color = treatment),
            linewidth = 0.4, linetype = "solid", alpha = 0.8) +
  
  # Bordures basses des ribbons
  geom_line(data = pred_df, aes(x = minute, y = conf.low, color = treatment),
            linewidth = 0.4, linetype = "solid", alpha = 0.8) +
  
  # facet_wrap(~ part) +
  labs(title = "Effets modélisés de la vitesse (par minute, côté et période)",
       x = "Temps (min)", y = "Temps passé dans la condition (%)", color = "Condition", fill = "Condition") +
  
  scale_color_manual(values = c(color_palette, border_palette)) +
  scale_fill_manual(values = fill_palette) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
