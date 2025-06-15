
                            
################################################################################
## Monte Carlo: Bias and MSE
################################################################################
seed_set <- 35151
amostras <- c(15, 30, 60, 120, 240, 480, 760)

simulate_monte_carlo_estimators <- function(
    models = c("gamma", "inv_gamma", "weibull", "inv_weibull"),
    ns = amostras,
    mu_vals = c(2),
    sigma_vals = c(1),
    delta = 1.5,
    B = 10000,
    alpha2 = 1/100, beta2 = 1/100,
    p = 1
) {
  results <- list()

  for (model in models) {
    for (mu in mu_vals) {
      for (sigma in sigma_vals) {
        for (n in ns) {

          est_mu <- numeric(B)
          est_sigma <- numeric(B)

          for (b in 1:B) {
            x <- tryCatch({
              rgeneric_expfam(n = n, model = model, mu = mu, sigma = sigma, delta = delta, p = p)
            }, error = function(e) return(NULL))

            if (is.null(x)) next

            est <- tryCatch({
              estimate_from_table1_model_analytic(x, model = model, delta = delta,
                                                  alpha2 = alpha2, beta2 = beta2, p = p)
            }, error = function(e) return(list(mu = NA, sigma = NA)))

            est_mu[b] <- est$mu
            est_sigma[b] <- est$sigma
          }

          df <- data.frame(
            model = model,
            n = n,
            mu = mu,
            sigma = sigma,
            bias_mu = abs((mean(est_mu, na.rm = TRUE) - mu)/mu),
            bias_sigma = abs((mean(est_sigma, na.rm = TRUE) - sigma)/sigma),
            mse_mu = mean((est_mu - mu)^2, na.rm = TRUE),
            mse_sigma = mean((est_sigma - sigma)^2, na.rm = TRUE)
          )

          results[[length(results) + 1]] <- df
        }
      }
    }
  }

  final_results <- do.call(rbind, results)
  return(final_results)
}

set.seed(seed_set)
results <- simulate_monte_carlo_estimators()
model_names <- c(
  gamma = "Gamma",
  inv_gamma = "Inverse Gamma",
  weibull = "Weibull",
  inv_weibull = "Inverse Weibull"
)

results <- results %>%
  mutate(Model = model_names[model]) %>%
  mutate(Model = factor(Model, levels = unique(model_names)))

linetypes <- c("Gamma" = "solid", "Inverse Gamma" = "dashed",
               "Weibull" = "dotdash", "Inverse Weibull" = "twodash")

colors <- c("Gamma" = "#1b9e77", "Inverse Gamma" = "#d95f02",
            "Weibull" = "#7570b3", "Inverse Weibull" = "#e7298a")

base_font_size <- 13

p1 <- ggplot(results, aes(x = n, y = bias_mu, color = Model, linetype = Model)) +
  geom_line(size = 1) + geom_point() +
  labs(title = expression("Relative Bias of the Estimator of " * mu),
       x = "Sample Size (n)", y = "Relative Bias") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  theme_classic(base_size = base_font_size) +
  theme(legend.title = element_text(size = base_font_size),
        legend.text = element_text(size = base_font_size),
        axis.title = element_text(size = base_font_size),
        axis.text = element_text(size = base_font_size))
ggsave("bias_mu.eps", plot = p1, device = "eps", width = 6, height = 4)

p2 <- ggplot(results, aes(x = n, y = mse_mu, color = Model, linetype = Model)) +
  geom_line(size = 1) + geom_point() +
  labs(title = expression("MSE of the Estimator of " * mu),
       x = "Sample Size (n)", y = "MSE") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  theme_classic(base_size = base_font_size) +
  theme(legend.title = element_text(size = base_font_size),
        legend.text = element_text(size = base_font_size),
        axis.title = element_text(size = base_font_size),
        axis.text = element_text(size = base_font_size))
ggsave("mse_mu.eps", plot = p2, device = "eps", width = 6, height = 4)

p3 <- ggplot(results, aes(x = n, y = bias_sigma, color = Model, linetype = Model)) +
  geom_line(size = 1) + geom_point() +
  labs(title = expression("Relative Bias of the Estimator of " * sigma),
       x = "Sample Size (n)", y = "Relative Bias") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  theme_classic(base_size = base_font_size) +
  theme(legend.title = element_text(size = base_font_size),
        legend.text = element_text(size = base_font_size),
        axis.title = element_text(size = base_font_size),
        axis.text = element_text(size = base_font_size))
ggsave("bias_sigma.eps", plot = p3, device = "eps", width = 6, height = 4)

p4 <- ggplot(results, aes(x = n, y = mse_sigma, color = Model, linetype = Model)) +
  geom_line(size = 1) + geom_point() +
  labs(title = expression("MSE of the Estimator of " * sigma),
       x = "Sample Size (n)", y = "MSE") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  theme_classic(base_size = base_font_size) +
  theme(legend.title = element_text(size = base_font_size),
        legend.text = element_text(size = base_font_size),
        axis.title = element_text(size = base_font_size),
        axis.text = element_text(size = base_font_size))
ggsave("mse_sigma.eps", plot = p4, device = "eps", width = 6, height = 4)


                            
################################################################################
## Monte Carlo: Compare Analytic vs MAP vs ML
################################################################################

simulate_comparison_mle_vs_custom <- function(
    model = c("gamma"),
    n     = 100,
    mu    = 2,
    sigma = 1,
    delta = 1.5,
    p     = 1,
    B     = 10000
) {
  model <- match.arg(model)
  # Pre-allocate
  res_analytic_mu    <- numeric(B)
  res_analytic_sigma <- numeric(B)
  res_map_mu         <- numeric(B)
  res_map_sigma      <- numeric(B)
  res_mle_mu         <- numeric(B)
  res_mle_sigma      <- numeric(B)
  
  for (b in seq_len(B)) {
    X <- rgeneric_expfam(n, model, mu, sigma, delta, p)
    
    # 1) Proposed analytic
    est_a <- estimate_from_table1_model_analytic(
      X, model, delta, p, alpha2 = 1/100, beta2 = 1/100)
    res_analytic_mu[b]    <- est_a$mu
    res_analytic_sigma[b] <- est_a$sigma
    
    # 2) Numeric MAP
    est_m <- estimate_map_numerical(
      X,
      T_fun         = function(x) x,
      T_prime       = function(x) 1,
      T_double_prime= function(x) 0,
      alpha1        = 1/100, beta1 = 1/100,
      alpha2        = 1/100, beta2 = 1/100,
      init          = c(mu = est_a$mu, sigma = est_a$sigma)
    )
    res_map_mu[b]    <- est_m$mu
    res_map_sigma[b] <- est_m$sigma
    
    # 3) Standard MLE
    est_l <- mle_gamma_mu_sigma(X, init = c(est_a$mu, est_a$sigma))
    res_mle_mu[b]    <- est_l$mu_hat
    res_mle_sigma[b] <- est_l$sigma_hat
  }
  
  data.frame(
    bias_mu_analytic    = abs(mean(res_analytic_mu - mu)/mu),
    bias_sigma_analytic = abs(mean(res_analytic_sigma - sigma)/sigma),
    mse_mu_analytic     = mean((res_analytic_mu - mu)^2),
    mse_sigma_analytic  = mean((res_analytic_sigma - sigma)^2),
    bias_mu_map         = abs(mean(res_map_mu - mu)/mu),
    bias_sigma_map      = abs(mean(res_map_sigma - sigma)/sigma),
    mse_mu_map          = mean((res_map_mu - mu)^2),
    mse_sigma_map       = mean((res_map_sigma - sigma)^2),
    bias_mu_mle         = abs(mean(res_mle_mu - mu)/mu),
    bias_sigma_mle      = abs(mean(res_mle_sigma - sigma)/sigma),
    mse_mu_mle          = mean((res_mle_mu - mu)^2),
    mse_sigma_mle       = mean((res_mle_sigma - sigma)^2)
  )
}

################################################################################
## Run full comparison over a vector of sample sizes and save plots
################################################################################
run_comparison <- function(sim_sizes,
                           model = "gamma",
                           mu    = 2,
                           sigma = 1,
                           delta = 1.5,
                           p     = 1,
                           B     = 10000,
                           seed  = 123) {
  set.seed(seed)
  library(tidyr)
  
  results_list <- lapply(sim_sizes, function(n) {
    df <- simulate_comparison_mle_vs_custom(model, n, mu, sigma, delta, p, B)
    df$n <- n
    df
  })
  comp_df <- bind_rows(results_list)
  
  # Reshape to long format
  long_df <- comp_df %>%
    pivot_longer(
      cols       = -n,
      names_to   = c("metric","param","est"),
      names_pattern = "(bias|mse)_(mu|sigma)_(analytic|map|mle)",
      values_to  = "value"
    ) %>%
    mutate(
      Metric    = toupper(metric),
      Estimator = recode(est,
                         analytic = "Proposed",
                         map      = "MAP",
                         mle      = "ML"),
      Parameter = ifelse(param=="mu","μ","σ")
    )
  
  # Plot settings
  colors <- c(Proposed = "#1b9e77", MAP = "#d95f02", ML = "#7570b3")
  lines  <- c(Proposed = "solid",    MAP = "dashed",    ML = "dotdash")
  
  # Relative bias
  p_bias <- ggplot(filter(long_df, Metric=="BIAS"),
                   aes(n, value, color=Estimator, linetype=Estimator)) +
    geom_line(size=1) + geom_point() +
    facet_wrap(~Parameter, scales="free_y") +
    scale_color_manual(values=colors) +
    scale_linetype_manual(values=lines) +
    labs(x="Sample Size (n)",
         y="Relative Bias",
         title="Relative Bias: Proposed vs MAP vs ML") +
    theme_classic(base_size=13) +
    theme(legend.position="top", legend.title=element_blank())
  
  # MSE
  p_mse <- ggplot(filter(long_df, Metric=="MSE"),
                  aes(n, value, color=Estimator, linetype=Estimator)) +
    geom_line(size=1) + geom_point() +
    facet_wrap(~Parameter, scales="free_y") +
    scale_color_manual(values=colors) +
    scale_linetype_manual(values=lines) +
    labs(x="Sample Size (n)",
         y="MSE",
         title="MSE: Proposed vs MAP vs ML") +
    theme_classic(base_size=13) +
    theme(legend.position="top", legend.title=element_blank())
  
  # Save
  ggsave("bias_comparison.eps", p_bias, width=7, height=4, device="eps")
  ggsave("mse_comparison.eps",  p_mse,  width=7, height=4, device="eps")
  
  invisible(list(bias_plot = p_bias, mse_plot = p_mse, data = comp_df))
}

