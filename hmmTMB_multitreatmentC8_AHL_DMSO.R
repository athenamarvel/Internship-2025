##########################################################################################
#### Preamble ####
##########################################################################################
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


raw_data <- read.csv("/Users/athena.marvel/Documents/Stage M1 2025/R script/dataF_AHL_DMSO.csv")


sub_data <- raw_data |> # example subset of 5 fish
  dplyr::mutate(time = floor(secs) + 1) |>
  dplyr::reframe(n_cuex = mean(cue_x),
                 n_y = mean(olf_y)
                 ,.by = c(treatment, part, fish, time)) |>
  dplyr::group_by(fish) |>
  dplyr::arrange(time, .by_group = TRUE) |> # ensure time is in correct order
  dplyr::ungroup() |>
  dplyr::mutate(ID = factor(fish)) |> 
  dplyr::rename(secs = time) |>
  na.omit() |>
  dplyr::mutate(time = 1:dplyr::n(),.by =fish) |>
  as.data.frame() # hmmTMB doesn't like tibbles

ggplot(sub_data, aes(x = secs, y = n_cuex)) + # visualise time series
  geom_line(aes(col = part)) +
  facet_grid(fish ~ treatment) +
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

par0 <- list(n_cuex = list(shape1 = shape1_0, shape2 = shape2_0)) # list of initial parameters
dists <- list(n_cuex = "beta") # choose model distribution family

f <- list(n_cuex = list(shape1 =  ~ 1 , shape2 =  ~ 1)) # define formulas for observation equation(s)

# create Observation object with data, distribution, number of states, initial values and formula
obs <- Observation$new(
  data = sub_data,
  dists = dists,
  n_states = 3,
  par = par0,
  formulas = f
)

# create initial transition matrix. Assume likely to stay in current state (0.8) and
# unlikely transition to new state (0.1). Rows must sum to 1
# if set 3->1 & 1->3 to 0, model doesn't fit properly
tpm0 <- matrix(c(0.8, 0.2, 0.0, 
                 0.1, 0.8, 0.1,
                 0.0, 0.2, 0.8),
               ncol = 3,
               byrow = TRUE)

# ensure that state values are fixed and shouldn't be estimated/changed
fixpar <- list(
  obs = c(
    "n_cuex.shape1.state1.(Intercept)" = NA,
    "n_cuex.shape1.state2.(Intercept)" = NA,
    "n_cuex.shape1.state3.(Intercept)" = NA,
    "n_cuex.shape2.state1.(Intercept)" = NA,
    "n_cuex.shape2.state2.(Intercept)" = NA,
    "n_cuex.shape2.state3.(Intercept)" = NA
  )
  ,hid = c("S1>S3.(Intercept)"  = NA,
          "S3>S1.(Intercept)"  = NA) # here specific transition probabilities can be fixed
)

# create MarkovChain object (i.e. the transition matrix part)
hid <- MarkovChain$new(
  n_states = 3,
  formula = ~ part*treatment + s(fish,bs = "re"),
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

qqnorm(pr$n_cuex)

# posterior predictive checks 
pp_check.hmm(hmm_tmb)

# extract the distribution of values that are classified as a "state".
# Should be the same as shape1_0 and shape2_0
state_means <- hmm_tmb$coeff_fe()[["obs"]] |> # extract coefficients
  as.data.frame() |> # convert from matrix to dataframe
  dplyr::rename(estimate = V1) |> # rename ambiguous column name
  dplyr::mutate(estimate = exp(estimate)) |> # exponentiate as estimate is made on log scale
  tibble::rownames_to_column(var = "variable") |> # add rownames as new column
  dplyr::mutate(variable = gsub("n_cuex\\.|\\.\\(Intercept\\)", "", variable)) |> #tidy rownames
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
  xlab("Probability") + ylab("n_cuex") +
  theme_classic()

# extract model classifications for each data point
tmb_viterbi <- sub_data |>
  dplyr::mutate(pred = hmm_tmb$viterbi())

ggplot(tmb_viterbi, aes(
  x = secs,
  y = n_cuex,
  col = factor(pred),
  group = fish
)) +
  geom_line() +
  scale_color_manual(values = c("#40A24F", "#5C0FF7", "#554378"),
                     name = "State") +
  xlab("Time") +
  facet_wrap(treatment~ fish) +
  theme_classic() 

##########################################################################################
#### Extract probability a fish is in each state ####
##########################################################################################

# Prédictions du modèle

library(ggdist)

# Prédictions du modèle
delta_preds <- hmm_tmb$predict(
  what = "delta",
  t = "all",
  newdata = expand.grid(
    "part" = unique(sub_data$part),
    "treatment" = unique(sub_data$treatment),
    "fish" = unique(sub_data$fish)
  ),
  level = 0.95,
  n_post = 1000,
  return_post = TRUE
)

# pred_data
list_pred_data <- lapply(seq_along(delta_preds$post), function(.draw) {
  as.data.frame(delta_preds$post[[.draw]]) %>%
    cbind(expand.grid(
      "part" = unique(sub_data$part),
      "treatment" = unique(sub_data$treatment),
      "fish" = unique(sub_data$fish)
    )) %>%
    mutate(.draw = .draw)
})
# Dans pred_data
pred_data <- do.call("rbind", list_pred_data) %>%
  pivot_longer(cols = contains("state"), names_to = "state", values_to = "prob") %>%
  group_by(part, treatment, state, fish) %>%
  ggdist::median_qi(.exclude = c(".chain", ".iteration", ".draw", ".row")) %>%
  ungroup() %>%
  rename(
    partie = part,
    traitement = treatment,
    probabilité = prob
  )

# Marg_pred avec les mêmes renommages
marg_pred <- do.call("rbind", list_marg_pred) %>%
  pivot_longer(cols = contains("state"), names_to = "state", values_to = "prob") %>%
  reframe(prob = mean(prob), .by = c(part, treatment, state, .draw)) %>%
  group_by(part, treatment, state) %>%
  ggdist::median_qi(.exclude = c(".chain", ".iteration", ".draw", ".row")) %>%
  ungroup() %>%
  rename(
    partie = part,
    traitement = treatment,
    probabilité = prob
  )

# Le plot avec les bons noms
ggplot(data = pred_data,
       aes(x = partie, y = probabilité,
           fill = traitement, color = traitement)) +
  geom_point(size = 2, alpha = 0.5, position = position_dodge(width = 0.25)) +
  geom_line(aes(group = interaction(fish, traitement, state)),
            alpha = 0.25, position = position_dodge(width = 0.25)) +
  geom_pointrange(data = marg_pred,
                  aes(ymin = .lower, ymax = .upper, fill = traitement, color = traitement),
                  size = 1, position = position_dodge(width = 0.25)) +
  scale_color_manual(values = c("#879BFF", "#626780")) +
  facet_wrap(~ state, labeller = as_labeller(c(
    "state 1" = "Repoussé",
    "state 2" = "Indifférent",
    "state 3" = "Attiré"
  ))) +
  theme_bw()


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
                        "treatment" = unique(sub_data$treatment), 
                        "fish" = unique(sub_data$fish)) |>
    dplyr::mutate(int = interaction(part,treatment)),
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
                                 "treatment" = unique(sub_data$treatment), 
                                 "fish" = unique(sub_data$fish)) |>
                       dplyr::mutate(covariate = 1:dplyr::n()) |>
                       dplyr::mutate(covariate = provideDimnames(as.matrix(covariate), sep = "", base = list(LETTERS), unique = TRUE) |>
                                       rownames()), 
                     by = "covariate") |>
    dplyr::mutate(.draw = .draw)}) |>
  do.call("rbind",args = _) |> # relabel covariate to match the true covariate values 
  dplyr::reframe(Freq = mean(Freq), .by = c(part, treatment, from, to, .draw)) |>
  dplyr::group_by(part, treatment, from, to) |>
  ggdist::median_qi(.exclude = c(".chain", ".iteration", ".draw", ".row")) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    of_interest = dplyr::if_else(
      from == "state1" &
        to == "state3" | from == "state3" & to == "state1",
      TRUE,
      FALSE
    ))

ggplot(marg_df, aes(x = part, y = Freq,shape = treatment)) +
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
      group = interaction(from, to, treatment),
      linetype = treatment
    ),
    colour = "black",
    linewidth = 0.5 , position = position_dodge(width = 0.25)) +
  ggdist::geom_pointinterval(aes(ymin = .lower, ymax = .upper),fill = "black", position = position_dodge(width = 0.25)) +
  xlab("Part") +
  ylab("Probability") +
  scale_fill_manual(values = c("white", "#5C0FF7"), guide = "none") +
  scale_shape_manual(values = c(21,1), name = "Treatment") +
  scale_linetype_manual(values = c("solid","dashed"), name = "Treatment") +
  guides(linetype = guide_legend(override.aes = list(linetype = 0))) +
  facet_wrap( ~ paste(from, to, sep = "->"), scales = "free") +
  theme_classic()
 
