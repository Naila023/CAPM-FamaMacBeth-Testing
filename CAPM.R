# -----------------------------------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------------------------------
library(tidyverse)
library(broom)
library(purrr)
library(systemfit)   # For SURE
library(gmm)         # For GMM
library(aod)         # For Wald test
library(GRS.test)    # For GRS test
library(lubridate)   # For date handling
library(car)         # For linearHypothesis

# -----------------------------------------------------------------------------
# Load the Dataset
# -----------------------------------------------------------------------------
df <- read_csv("USstocks_balanced.csv")

# -----------------------------------------------------------------------------
# Define Base Periods and Number of Periods (6 Periods)
# -----------------------------------------------------------------------------
base_formation_start  <- 1980
base_formation_end    <- 1986
base_estimation_start <- 1987
base_estimation_end   <- 1991
base_testing_start    <- 1992
base_testing_end      <- 1995

n_periods <- 6

# -----------------------------------------------------------------------------
# Fama–MacBeth Procedure: Portfolio Formation, Estimation, Testing, and Summary
# -----------------------------------------------------------------------------
results <- list()

for(i in 0:(n_periods - 1)){
  
  # Define Period Boundaries for the Current Iteration
  formation_start  <- base_formation_start + 4 * i
  formation_end    <- base_formation_end + 4 * i
  estimation_start <- base_estimation_start + 4 * i
  estimation_end   <- base_estimation_end + 4 * i
  testing_start    <- base_testing_start + 4 * i
  testing_end      <- base_testing_end + 4 * i
  
  # Portfolio Formation: Estimate each firm's beta in the formation period and rank into 20 portfolios
  pf_groups <- df %>% 
    filter(between(year, formation_start, formation_end)) %>% 
    group_by(permno) %>% 
    nest() %>% 
    mutate(betaest = map(data, ~ coef(lm(ri ~ rm, data = .x))[2])) %>% 
    unnest(betaest) %>% 
    ungroup() %>% 
    mutate(pfid = ntile(betaest, 20)) %>% 
    select(permno, betaest, pfid)
  
  # Merge Portfolio Assignments with Data for Estimation and Testing Periods
  df_withpfid <- df %>% 
    filter(between(year, estimation_start, testing_end)) %>% 
    left_join(pf_groups, by = "permno")
  
  # Portfolio Estimation: For each firm, re-estimate beta and compute residual std. error (sigma)
  n_years <- testing_end - testing_start + 1
  beta_sigma_period <- data.frame()
  for(j in 0:(n_years - 1)){
    beta_sigma <- df_withpfid %>% 
      filter(year <= (estimation_end + j)) %>% 
      group_by(permno) %>% 
      nest() %>% 
      mutate(
        betas = map(data, ~ coef(lm(ri ~ rm, data = .x))[2]),
        sigmas = map(data, ~ sigma(lm(ri ~ rm, data = .x))),
        yeartesting = testing_start + j
      ) %>% 
      unnest(c(betas, sigmas)) %>% 
      ungroup() %>% 
      select(permno, betas, sigmas, yeartesting)
    beta_sigma_period <- bind_rows(beta_sigma_period, beta_sigma)
  }
  beta_sigma_period <- rename(beta_sigma_period, year = yeartesting)
  
  # Merge Beta and Sigma Estimates with Testing Period Data
  df_testingperiod <- df_withpfid %>% 
    filter(between(year, testing_start, testing_end)) %>% 
    left_join(beta_sigma_period, by = c("permno", "year"))
  
  # Portfolio Testing: Compute Monthly Portfolio Averages
  testing_period <- df_testingperiod %>% 
    group_by(pfid, year, month) %>% 
    summarise(
      rp = mean(ri, na.rm = TRUE),   # Average portfolio return
      betap = mean(betas, na.rm = TRUE),
      betap_squared = mean(betas^2, na.rm = TRUE),
      sigmap = mean(sigmas, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Run Month-by-Month Cross-Sectional Regressions:
  # Regress lead(rp) on portfolio averages: betap, betap_squared, and sigmap.
  testing_period_gammas <- testing_period %>% 
    group_by(year, month) %>% 
    nest(data = !c(year, month)) %>% 
    mutate(
      estim = map(data, ~ lm(lead(rp) ~ betap + betap_squared + sigmap, data = .x)),
      estim = map(estim, tidy)
    ) %>% 
    unnest(estim) %>% 
    select(year, month, term, estimate) %>% 
    ungroup()
  
  # Summarize Cross-Sectional Regression Results
  summary_results <- testing_period_gammas %>% 
    group_by(term) %>% 
    summarise(
      mean_est = mean(estimate, na.rm = TRUE),
      sd_est   = sd(estimate, na.rm = TRUE),
      n        = n(),
      t_stat   = mean_est / (sd_est / sqrt(n))
    ) %>% 
    ungroup()
  
  # Save results for the current period
  results[[as.character(i+1)]] <- list(
    period = i + 1,
    formation = c(formation_start, formation_end),
    estimation = c(estimation_start, estimation_end),
    testing = c(testing_start, testing_end),
    df_testingperiod = df_testingperiod,
    gammas = testing_period_gammas,
    summary = summary_results
  )
}

# Combine and Print Summary Statistics Across All Periods
all_summary <- map_df(results, ~ .x$summary, .id = "period")
print(all_summary, n = 24)

# Hypothesis Testing on the Cross-Sectional Regression Estimates
hypothesis_tests <- map(results, function(res) {
  gammas <- res$gammas
  # Test for non-linearity (beta^2 term)
  gamma2 <- gammas %>% filter(term == "betap_squared") %>% pull(estimate)
  t_test_gamma2 <- t.test(gamma2, mu = 0)
  # Test for additional risk (sigma)
  gamma3 <- gammas %>% filter(term == "sigmap") %>% pull(estimate)
  t_test_gamma3 <- t.test(gamma3, mu = 0)
  # Test for positive return-risk tradeoff (beta coefficient)
  gamma1 <- gammas %>% filter(term == "betap") %>% pull(estimate)
  gamma0 <- gammas %>% filter(term == "(Intercept)") %>% pull(estimate)
  rm_val <- res$df_testingperiod %>% 
    filter(permno == unique(res$df_testingperiod$permno)[1]) %>% 
    pull(rm)
  rm_r0 <- rm_val - gamma0
  t_test_gamma1 <- t.test(gamma1, rm_r0)
  # Sharpe-Lintner test (comparing intercept to rf)
  rf_val <- res$df_testingperiod %>% 
    filter(permno == unique(res$df_testingperiod$permno)[1]) %>% 
    pull(rf)
  t_test_gamma0 <- t.test(gamma0, rf_val)
  # Market efficiency tests via Ljung-Box on each coefficient
  lb_gamma1 <- Box.test(gamma1, lag = 1, type = "Ljung-Box")
  lb_gamma2 <- Box.test(gamma2, lag = 1, type = "Ljung-Box")
  lb_gamma3 <- Box.test(gamma3, lag = 1, type = "Ljung-Box")
  
  list(
    "Linearity Test" = t_test_gamma2,
    "No Systematic Non-Beta Risk" = t_test_gamma3,
    "Positive Return-Risk Tradeoff" = t_test_gamma1,
    "Sharpe-Lintner Test" = t_test_gamma0,
    "Market Efficiency Test - 1 (Ljung-Box)" = lb_gamma1,
    "Market Efficiency Test - 2 (Ljung-Box)" = lb_gamma2,
    "Market Efficiency Test - 3 (Ljung-Box)" = lb_gamma3
  )
})
print(hypothesis_tests)

# -----------------------------------------------------------------------------
# Modern Tests of Factor Models using Portfolio Excess Returns
# -----------------------------------------------------------------------------

# Compute Excess Portfolio Returns (EPR) and Excess Market Returns (EMR)
excess_returns_results <- list()
for (i in 1:n_periods) {
  testing_period_data <- results[[as.character(i)]]$df_testingperiod
  excess_returns <- testing_period_data %>%
    group_by(pfid, year, month) %>%
    summarise(
      rp = mean(ri, na.rm = TRUE),      # Portfolio return
      erp = rp - mean(rf, na.rm = TRUE),  # Excess Portfolio Return
      erm = mean(rm, na.rm = TRUE) - mean(rf, na.rm = TRUE),  # Excess Market Return
      .groups = "drop"
    )
  excess_returns_results[[as.character(i)]] <- excess_returns
}

all_excess_returns <- bind_rows(excess_returns_results, .id = "period") %>%
  mutate(date = lubridate::ym(paste0(year, "-", month)))

# -----------------------------------------------------------------------------
# LS Approach: Estimate SCL for each portfolio
# -----------------------------------------------------------------------------
ls_test_modern <- all_excess_returns %>% 
  group_by(pfid) %>% 
  nest(data = !pfid) %>% 
  mutate(estim = map(data, ~ lm(erp ~ erm, data = .x)),
         estim = map(estim, tidy)) %>% 
  unnest(estim) %>% 
  select(-data) %>% 
  arrange(pfid) %>% 
  filter(term == "(Intercept)")

knitr::kable(ls_test_modern, digits = 4)

# -----------------------------------------------------------------------------
# SURE Approach - Independent Equations
# -----------------------------------------------------------------------------
sure_independent <- all_excess_returns %>% 
  group_by(pfid) %>% 
  nest(data = !pfid) %>% 
  mutate(estim = map(data, ~lm(erp ~ erm, data = .x)),
         estim = map(estim, tidy)) %>% 
  unnest(estim) %>% 
  select(-data) %>% 
  ungroup() %>% 
  filter(term == "(Intercept)") %>% 
  summarise(wald = sum(statistic^2))
print(sure_independent)
print(qchisq(p = 0.05, df = 20, lower.tail = FALSE))
print(pchisq(2.13, df = 20, lower.tail = FALSE))

# -----------------------------------------------------------------------------
# SURE Approach - Dependent Equations
# -----------------------------------------------------------------------------
df_final_wide <- all_excess_returns %>%
  select(pfid, date, erp, erm) %>%
  pivot_wider(names_from = pfid, values_from = erp, names_glue = "erp_{pfid}") %>%
  arrange(date)
print(head(df_final_wide))

system_equations <- map(1:20, ~paste0("erp_", .x, " ~ erm"))
system_equations <- map(system_equations, as.formula)
print(system_equations)

results_sure_modern <- systemfit(system_equations, method = "SUR", data = df_final_wide)
summary(results_sure_modern)

waldtest_sure_modern <- aod::wald.test(
  Sigma = vcov(results_sure_modern),
  b = coef(results_sure_modern),
  Terms = seq(1, 39, by = 2)
)
print(waldtest_sure_modern)
print(coef(results_sure_modern))
print(coef(results_sure_modern)[1])
print(coef(results_sure_modern)[2])
print(coef(results_sure_modern)[5])
print(seq(1, 39, by = 2))

# -----------------------------------------------------------------------------
# Convert Excess Portfolio Returns into Matrix Form
# -----------------------------------------------------------------------------
erp_matrix <- df_final_wide %>% 
  select(contains("erp")) %>%  
  as.matrix()
erm_matrix <- df_final_wide %>% 
  select(erm) %>%  
  as.matrix()

# -----------------------------------------------------------------------------
# GMM Approach
# -----------------------------------------------------------------------------
results_gmm_modern <- gmm(erp_matrix ~ erm_matrix,  
                          x = cbind(1, erm_matrix),  
                          vcov = "HAC",  
                          type = "twoStep")
summary(results_gmm_modern)
n_coef <- length(coef(results_gmm_modern))
print(paste("Number of estimated coefficients:", n_coef))
R <- diag(n_coef)
c_vec <- rep(0, n_coef)
if (nrow(R) != length(coef(results_gmm_modern))) {
  stop("Error: Dimensions of R and coefficient vector do not match!")
}
waldtest_gmm_modern <- aod::wald.test(
  Sigma = vcov(results_gmm_modern),  
  b = coef(results_gmm_modern),  
  Terms = 1:n_coef  
)
print(waldtest_gmm_modern)
print(coef(results_gmm_modern))
print(coef(results_gmm_modern)[1])
print(coef(results_gmm_modern)[2])
print(coef(results_gmm_modern)[5])
print(seq(1, n_coef, by = 2))

# Test for Restrictions on the Intercepts using linearHypothesis()
R_intercept <- diag(n_coef)
c_intercept <- rep(0, n_coef)
if (ncol(R_intercept) != length(coef(results_gmm_modern))) {
  stop("Error: Dimensions of R_intercept and coefficient vector do not match!")
}
linear_test_result <- car::linearHypothesis(results_gmm_modern, R_intercept, c_intercept, test = "Chisq")
print(linear_test_result)

# -----------------------------------------------------------------------------
# J-Test for Over identification Restrictions
# -----------------------------------------------------------------------------
results_gmm_jtest <- gmm(erp_matrix ~ erm_matrix - 1, x = cbind(1, erm_matrix))
j_test_result <- specTest(results_gmm_jtest)
print(j_test_result)

# -----------------------------------------------------------------------------
# GRS Test
# -----------------------------------------------------------------------------
grstest_results <- GRS.test(erp_matrix, erm_matrix)
print(grstest_results$GRS.stat)
print(grstest_results$GRS.pval)
