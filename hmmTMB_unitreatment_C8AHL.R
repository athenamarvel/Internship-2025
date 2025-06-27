##########################################################################################
#### Preamble ####
##########################################################################################

library("readxl")
library("dplyr")
library("janitor")
library("ggplot2")
library("patchwork")
library("hmmTMB")

# user defined functions for converting between beta shape parameters and mean/precision
shape_to_prec <- function(shape1, shape2) {
  mu <- shape1 / (shape1 + shape2)
  phi <- shape1 + shape2
  return(list(mu = mu, phi = phi))
}

prec_to_shape <- function(mu, phi) {
  shape1 <- mu * phi
  shape2 <- (1 - mu) * phi
  return(list(shape1 = shape1, shape2 = shape2))
}

# posterior predictive check
pp_check.hmm <- function(model, data = model$obs()$data(), draws = 10, cores = getOption("mc.cores", 2L)){
  
  mod_data <- model$obs()$data()
  resps <- names(model$obs()$dists())
  draws <- parallel::mclapply(1:draws,function(draw){
    model$simulate(n = NROW(mod_data),
                   data = data,
                   silent = TRUE) |>
      dplyr::mutate(.draw = draw)
    
  },mc.cores = cores) |>
    do.call("rbind",args = _) |>
    tidyr::pivot_longer(dplyr::all_of(resps),names_to = "variable",values_to = ".value")
  
  ggplot(draws, aes(x = .value)) +
    geom_density(aes(group = .draw),alpha = 0.3,col = "lightblue") +
    geom_density(data = mod_data |>
                   tidyr::pivot_longer(dplyr::all_of(resps),names_to = "variable",values_to = ".value")
                 ,col = "black") +
    facet_wrap(~variable,scales = "free") + theme_classic()
  
}

folder_path <- "/Users/athena.marvel/Documents/Stage M1 2025/R data/C8AHLdata"

# Lister tous les fichiers CSV dans le dossier
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Lire et combiner tous les fichiers CSV
dataF <- do.call(rbind, lapply(file_list, read.csv))

dataF$minute <- substr(dataF$time, 1, 2)

dataF$secondeTEST <- substr(dataF$time, 4, 5)
dataF$secsTEST <- substr(dataF$time, 7, 8)

dataF$minute<- as.numeric(dataF$minute)
dataF <- dataF %>%
  mutate(seconde = 60*as.numeric(minute) + as.numeric(secondeTEST) + as.numeric(secsTEST)*0.01)%>%
  mutate(minute = minute +1)
dataF <- dataF %>%
  rename(fish = fish_replicate)

dataF <- dataF %>%
  mutate(cue_x = ifelse(cue_side == "Right",olf_x,1-olf_x))%>%
  mutate(fishPositionLR = ifelse(olf_x < 0.5, 'Left', 'Right')) %>%  
  mutate(fishPositionCue = ifelse(cue_x < 0.5, 'control', 'cue')) %>%
  mutate(part = ifelse(seconde >= 600, "After","Before"))%>%
  mutate(part = factor(part, levels = c("Before", "After")))

dataF <- dataF %>%
  group_by(fish)%>%
  mutate(max_frame = max(frame))%>%
  mutate(timelag = ifelse(frame == 0, NA, 1200/max_frame))

dataF <- dataF %>% 
  mutate(seconde = as.numeric(seconde))

sub_data <- dataF %>% 
  mutate(time = floor(seconde) + 1) %>% 
  group_by(part, fish, time) %>% 
  summarise(
    cue_x = mean(cue_x, na.rm = TRUE),
    n_y = mean(olf_y, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  group_by(fish) %>% 
  arrange(time, .by_group = TRUE) %>%  # ensure time is in correct order
  ungroup() %>% 
  mutate(ID = factor(fish)) %>% 
  rename(secs = time) %>% 
  na.omit() %>% 
  group_by(fish) %>% 
  mutate(time = row_number()) %>% 
  ungroup() %>% 
  as.data.frame()

sub_data <- sub_data %>%
  mutate(fish = as.factor(fish))

ggplot(sub_data, aes(x = time, y = cue_x)) + # visualise time series
  geom_line(aes(col = part)) +
  facet_wrap(~fish) +
  theme_classic()

##########################################################################################
#### Prepare model ####
##########################################################################################

# prepare probability distributions of 3 n_mpx states (low, mid, high)
# beta distributions are a bit finicky to define in my opinion and consist of two shape parameters
# we can convert from mean and precision (a.k.a phi) to shape using the equations:
# shape1 = mean*phi
# shape2 = (1-mean)*phi
# phi <- 50 # phi influences how wide the distribution is. 50 is wider than 500
# phi_mid <- 500
# mean_low <- 0.25
# mean_mid <- 0.5
# mean_high <- 0.75
# shape1_0 <- c(mean_low * phi, mean_mid * phi_mid, mean_high * phi)
# shape2_0 <- c((1 - mean_low) * phi, (1 - mean_mid) * phi_mid, (1 - mean_high) * phi)

shape1_0 <- c(3,200, 10)
shape2_0 <- c(10,200, 3)

par0 <- list(cue_x = list(shape1 = shape1_0, shape2 = shape2_0)) # list of initial parameters
dists <- list(cue_x = "beta") # choose model distribution family

f <- list(cue_x = list(shape1 =  ~ 1 , shape2 =  ~ 1)) # define formulas for observation equation(s)

# create Observation object with data, distribution, number of states, initial values and formula
obs <- Observation$new(
  data =  sub_data,
  dists = dists,
  n_states = 3,
  par = par0,
  formulas = f
)

# create initial transition matrix. Assume likely to stay in current state (0.8) and
# unlikely transition to new state (0.1). Rows must sum to 1.
# if set 3->1 & 1->3 to 0, model doesn't fit properly
tpm0 <- matrix(c(0.8, 0.1, 0.1, 
                 0.1, 0.8, 0.1,
                 0.1, 0.1, 0.8),
               ncol = 3,
               byrow = TRUE)

# ensure that state values are fixed and shouldn't be estimated/changed
fixpar <- list(
  obs = c(
    "cue_x.shape1.state1.(Intercept)" = NA,
    "cue_x.shape1.state2.(Intercept)" = NA,
    "cue_x.shape1.state3.(Intercept)" = NA,
    "cue_x.shape2.state1.(Intercept)" = NA,
    "cue_x.shape2.state2.(Intercept)" = NA,
    "cue_x.shape2.state3.(Intercept)" = NA
  )
  # ,hid = c("S1>S3.(Intercept)"  = NA,
  #         "S3>S1.(Intercept)"  = NA) # here specific transition probabilities can be fixed
)

# create MarkovChain object (i.e. the transition matrix part)
hid <- MarkovChain$new(
  n_states = 3,
  formula = ~ part + s(fish,bs = "re"),
  # the formula to change transition matrix given covariates
  data = sub_data
  ,tpm = tpm0 # not necessarily needed in my experience
  ,initial_state = "stationary" # needed to improve convergence
)
# create HMM object
hmm_tmb <- HMM$new(obs = obs, hid = hid, fixpar = fixpar)

##########################################################################################
#### Fit model ####
##########################################################################################

# fit HMM object. I've chosen the "BFGS" optimiser as I think it's more effective
hmm_tmb$fit(silent = T, method = "BFGS") # silent can be set to TRUE to suppress readout

# we could also fit the same model Bayesian
# hmm_tmb$fit_stan(chains = 4,
#                  iter = 250,
#                  cores = 4,
#                  silent = FALSE)
#rstan::stan_trace(hmm_tmb$out_stan())

##########################################################################################
#### Extract model coefficients ####
##########################################################################################

# create QQ plot (indicates not great fit)
pr <- hmm_tmb$pseudores()
qqnorm(pr$cue_x)

# posterior predictive checks 
pp_check.hmm(hmm_tmb)

# extract the distribution of values that are classified as a "state".
# Should be the same as shape1_0 and shape2_0
state_means <- hmm_tmb$coeff_fe()[["obs"]] |> # extract coefficients
  as.data.frame() |> # convert from matrix to dataframe
  dplyr::rename(estimate = V1) |> # rename ambiguous column name
  dplyr::mutate(estimate = exp(estimate)) |> # exponentiate as estimate is made on log scale
  tibble::rownames_to_column(var = "variable") |> # add rownames as new column
  dplyr::mutate(variable = gsub("cue_x\\.|\\.\\(Intercept\\)", "", variable)) |> #tidy rownames
  tidyr::separate_wider_delim(
    cols = variable,
    delim = ".",
    names = c("variable", "state")
  ) |> # split rowname in to separate columns around the "." delimiter  
  tidyr::pivot_wider(id_cols = state, names_from = "variable", values_from = "estimate") |> # pivot the parameter names in to new columns
  dplyr::mutate(mean = shape_to_prec(shape1, shape2)[[1]],
                prec = shape_to_prec(shape1, shape2)[[2]]) |> # convert from shape to mean to double check it matches our initial aims
  dplyr::mutate(.value = list(rbeta(
    1000, shape1 = shape1, shape2 = shape2
  )), .by = state) |> # use shape parameters to simulate a beta distribution
  tidyr::unnest(.value) |>
  dplyr::mutate(.draw = 1:dplyr::n(), .by = state)

# plot state distributions
ggplot(state_means, aes(
  y = .value,
  group = state,
  col = state,
  fill = state
)) +
  ggdist::stat_slab(alpha = 0.5, normalize = "groups") +
  scale_color_manual(values = c("#40A24F", "#5C0FF7", "#554378"),
                     name = "State") +
  scale_fill_manual(values = c("#40A24F", "#5C0FF7", "#554378"),
                    name = "State") +
  xlab("Probability") + ylab("cue_x") +
  theme_classic()

# extract model classifications for each data point
tmb_viterbi <- sub_data |>
  dplyr::mutate(pred = hmm_tmb$viterbi())

ggplot(tmb_viterbi, aes(
  x = secs,
  y = cue_x,
  col = factor(pred),
  group = fish
)) +
  geom_line() +
  scale_color_manual(values = c("#40A24F", "#5C0FF7", "#554378"),
                     name = "State") +
  xlab("Time") +
  facet_wrap( ~ fish) +
  theme_classic() +
  labs(title= "Comportement des poissons exposés à une solution de menthe ")

##########################################################################################
#### Extract probability a fish is in each state ####
##########################################################################################

# 1. Modifier sub_data avant la prédiction
sub_data$part <- dplyr::recode(sub_data$part,
                               "Before" = "Avant",
                               "After" = "Après")

# 2. Prédiction avec les nouvelles valeurs
delta_preds <-  hmm_tmb$predict(what = "delta",
                                t = "all",
                                newdata = expand.grid("part" = unique(sub_data$part),
                                                      "fish" = unique(sub_data$fish)),
                                level = 0.95,
                                n_post = 1000,
                                return_post = TRUE)

# 3. Création de pred_data
pred_data <- lapply(seq_along(delta_preds$post), function(.draw) { 
  as.data.frame(delta_preds$post[[.draw]]) |>
    cbind(expand.grid("part" = unique(sub_data$part),
                      "fish" = unique(sub_data$fish))) |>
    dplyr::mutate(.draw = .draw)
}) |>
  do.call("rbind", args = _) |>
  tidyr::pivot_longer(cols = dplyr::contains("state"), names_to = "state", values_to = "prob") |>
  dplyr::group_by(part, state, fish) |>
  ggdist::median_qi(.exclude = c(".chain", ".iteration", ".draw", ".row")) |>
  dplyr::ungroup() |>
  dplyr::mutate(Probabilité = prob) |>
  dplyr::select(-prob)

# 4. Création de marg_pred
marg_pred <- lapply(seq_along(delta_preds$post), function(.draw) {
  as.data.frame(delta_preds$post[[.draw]]) |>
    cbind(expand.grid("part" = unique(sub_data$part),
                      "fish" = unique(sub_data$fish))) |>
    dplyr::mutate(.draw = .draw)
}) |>
  do.call("rbind", args = _) |>
  tidyr::pivot_longer(cols = dplyr::contains("state"), names_to = "state", values_to = "prob") |>
  dplyr::reframe(prob = mean(prob), .by = c(part, state, .draw)) |>
  dplyr::group_by(part, state) |>
  ggdist::median_qi(.exclude = c(".chain", ".iteration", ".draw", ".row")) |>
  dplyr::ungroup() |>
  dplyr::mutate(Probabilité = prob) |>
  dplyr::select(-prob)

# 5. Visualisation
state_labels <- c("state 1" = "Repoussé", "state 2" = "Indifférent", "state 3" = "Attiré")

ggplot(data = pred_data,
       aes(x = part, y = Probabilité)) +
  geom_point(size = 2, alpha = 0.5, position = position_dodge(width = 0.25)) +
  geom_line(aes(group = interaction(fish,state)), alpha = 0.25, position = position_dodge(width = 0.25)) +
  geom_pointrange(data = marg_pred, aes(ymin = .lower, ymax = .upper, color = part),
                  size = 1, position = position_dodge(width = 0.25)) +
  facet_wrap(~ state, labeller = as_labeller(state_labels)) +
  theme_bw() +
  labs(title = "Effet d'une solution d'AHL sur les états comportementaux des poissons",
       x = "Partie",
       y = "Probabilité")


##########################################################################################
#### Extract transition matrix ####
##########################################################################################

# extract estimated transition matrix across all time points and for each "part" separately
tpm_preds <- hmm_tmb$predict(
  what = "tpm",
  # "tpm" = transition matrix
  t = "all",
  # which time points
  newdata = expand.grid("part" = unique(sub_data$part),
                        "fish" = unique(sub_data$fish)),
  # this will need to change to match the covariates
  level = 0.95,
  # confidence interval
  n_post = 1000 # number of simulations to bootstrap the confidence interval
  , return_post = TRUE)

marg_df <- lapply(seq_along(tpm_preds$post),function(.draw){
  as.data.frame.table(tpm_preds$post[[.draw]])  |>
    dplyr::rename(from = Var1,
                  to = Var2,
                  covariate = Var3) |> # rename ambiguous column names
    dplyr::mutate(dplyr::across(from:to,  ~ gsub(" ", "", .x))) |> # remove white space
    dplyr::left_join(expand.grid("part" = unique(sub_data$part),
                                 "fish" = unique(sub_data$fish)) |>
                       dplyr::mutate(covariate = 1:dplyr::n()) |>
                       dplyr::mutate(covariate = provideDimnames(as.matrix(covariate), sep = "", base = list(LETTERS), unique = TRUE) |>
                                       rownames()), 
                     by = "covariate") |>
    dplyr::mutate(.draw = .draw)}) |>
  do.call("rbind",args = _) |> # relabel covariate to match the true covariate values 
  dplyr::reframe(Freq = mean(Freq), .by = c(part, from, to, .draw)) |>
  dplyr::group_by(part, from, to) |>
  ggdist::median_qi(.exclude = c(".chain", ".iteration", ".draw", ".row")) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    of_interest = dplyr::if_else(
      from == "state1" &
        to == "state3" | from == "state3" & to == "state1",
      TRUE,
      FALSE
    ))


ggplot(marg_df, aes(x = part, y = Freq)) +
  geom_rect(
    aes(fill = of_interest),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.01
  ) +
  geom_line(
    aes(
      x  = as.numeric(part),
      y = Freq,
      group = interaction(from, to)
    ),
    colour = "black",
    linewidth = 0.5 , position = position_dodge(width = 0.25)) +
  ggdist::geom_pointinterval(aes(ymin = .lower, ymax = .upper),fill = "black", position = position_dodge(width = 0.25)) +
  xlab("Part") +
  ylab("Probability") +
  scale_fill_manual(values = c("white", "#5C0FF7"), guide = "none") +
  scale_linetype_manual(values = c("solid","dashed"), name = "Treatment") +
  guides(linetype = guide_legend(override.aes = list(linetype = 0))) +
  facet_wrap( ~ paste(from, to, sep = "->"), scales = "free") +
  theme_classic()

