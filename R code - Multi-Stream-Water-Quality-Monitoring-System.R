# =========================
# Step 0: Libraries
# =========================
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(foreach)
library(doParallel)
library(openxlsx)
library(tibble)

# =========================
# Step 1: Load and preprocess
# =========================
df_raw <- read_csv("brisbane_water_quality.csv")

# Column names (as in your file)
col_timestamp      <- "Timestamp"
col_record_number  <- "Record number"
col_pH             <- "pH"
col_do_sat         <- "Dissolved Oxygen (%Saturation)"
col_turbidity      <- "Turbidity"
col_salinity       <- "Salinity"
col_speed          <- "Average Water Speed"
col_direction      <- "Average Water Direction"

# Clean and keep only needed cols; use all_of() to avoid tidyselect warnings
df_cleaned <- df_raw %>%
  mutate(Timestamp = ymd_hms(.data[[col_timestamp]], tz = "UTC")) %>%
  select(
    Timestamp,
    all_of(col_record_number),
    all_of(col_pH),
    all_of(col_do_sat),
    all_of(col_turbidity),
    all_of(col_salinity),
    all_of(col_speed),
    all_of(col_direction)
  ) %>%
  drop_na()

# Split 2023 (Phase I) and 2024 (Phase II)
df_2023 <- df_cleaned %>% filter(year(Timestamp) == 2023)
df_2024 <- df_cleaned %>% filter(year(Timestamp) == 2024)
#View(df_2024)

cat("\n[Step 1] Data loaded and cleaned.\n")
cat("Phase I (2023) records:", nrow(df_2023), "\n")
cat("Phase II (2024) records:", nrow(df_2024), "\n")

# =========================
# Step 2: 1st/99th percentile bands from 2023
# =========================
calculate_percentiles <- function(df) {
  df %>%
    summarise(
      pH_1st  = quantile(.data[[col_pH]],      0.01, na.rm = TRUE),
      pH_99th = quantile(.data[[col_pH]],      0.99, na.rm = TRUE),
      
      DO_Sat_1st  = quantile(.data[[col_do_sat]], 0.01, na.rm = TRUE),
      DO_Sat_99th = quantile(.data[[col_do_sat]], 0.99, na.rm = TRUE),
      
      Turbidity_1st  = quantile(.data[[col_turbidity]], 0.01, na.rm = TRUE),
      Turbidity_99th = quantile(.data[[col_turbidity]], 0.99, na.rm = TRUE),
      
      Salinity_1st  = quantile(.data[[col_salinity]], 0.01, na.rm = TRUE),
      Salinity_99th = quantile(.data[[col_salinity]], 0.99, na.rm = TRUE),
      
      Speed_1st  = quantile(.data[[col_speed]], 0.01, na.rm = TRUE),
      Speed_99th = quantile(.data[[col_speed]], 0.99, na.rm = TRUE),
      
      Direction_1st  = quantile(.data[[col_direction]], 0.01, na.rm = TRUE),
      Direction_99th = quantile(.data[[col_direction]], 0.99, na.rm = TRUE)
    )
}

percentiles_2023 <- calculate_percentiles(df_2023)

# Binary mapping: 1 if outside band, 0 otherwise
convert_to_binary <- function(x, min_val, max_val) {
  ifelse(x < min_val | x > max_val, 1L, 0L)
}

# Build 2024 binary indicators
df_2024_binary <- df_2024 %>%
  mutate(
    x_pH        = convert_to_binary(.data[[col_pH]],        percentiles_2023$pH_1st,       percentiles_2023$pH_99th),
    x_DO_Sat    = convert_to_binary(.data[[col_do_sat]],    percentiles_2023$DO_Sat_1st,   percentiles_2023$DO_Sat_99th),
    x_Turbidity = convert_to_binary(.data[[col_turbidity]], percentiles_2023$Turbidity_1st,percentiles_2023$Turbidity_99th),
    x_Salinity  = convert_to_binary(.data[[col_salinity]],  percentiles_2023$Salinity_1st, percentiles_2023$Salinity_99th),
    x_Speed     = convert_to_binary(.data[[col_speed]],     percentiles_2023$Speed_1st,    percentiles_2023$Speed_99th),
    x_Direction = convert_to_binary(.data[[col_direction]], percentiles_2023$Direction_1st,percentiles_2023$Direction_99th)
  )

cat("\n[Step 2] 2024 binaries created from 2023 1st/99th percentiles.\n")

# =========================
# Step 3: Estimate p0 from 2024 binaries; build C_j, Q_t, W_t
# =========================
# p0 per stream (proportion of 1s)
p0_tbl <- df_2024_binary %>%
  summarise(
    p0_pH        = mean(x_pH,        na.rm = TRUE),
    p0_DO_Sat    = mean(x_DO_Sat,    na.rm = TRUE),
    p0_Turbidity = mean(x_Turbidity, na.rm = TRUE),
    p0_Salinity  = mean(x_Salinity,  na.rm = TRUE),
    p0_Speed     = mean(x_Speed,     na.rm = TRUE),
    p0_Direction = mean(x_Direction, na.rm = TRUE)
  )
print(p0_tbl)

# fixed order vector (pH, DO, Turbidity, Salinity, Speed, Direction)
p0_vec <- as.numeric(p0_tbl[1, c("p0_pH","p0_DO_Sat","p0_Turbidity","p0_Salinity","p0_Speed","p0_Direction")])

# C_j, Q_t, W_t with μ = sum p0_i and σ² = sum p0_i(1-p0_i)
ensure_W_cols <- function(df, p0_vec, id_col_name = col_record_number) {
  Cj <- with(df, x_pH + x_DO_Sat + x_Turbidity + x_Salinity + x_Speed + x_Direction)
  mu     <- sum(p0_vec)
  sigma2 <- sum(p0_vec * (1 - p0_vec))
  t_idx  <- seq_along(Cj)
  Qt     <- cumsum(Cj)
  Wt     <- (Qt - mu * t_idx) / sqrt(t_idx * sigma2)
  
  df %>%
    mutate(
      C_j = Cj,
      Q_t = Qt,
      W_t = Wt
    ) %>%
    # keep Timestamp and Record number for reporting
    rename(`Record number` = all_of(id_col_name))
}

df_2024_w <- ensure_W_cols(df_2024_binary, p0_vec)

cat("\n[Step 3] Built C_j, Q_t, W_t from df_2024_binary using estimated p0.\n")
#view(df_2024_w)
# =========================
# Step 4: CSB-EWMA on W_t and signal detection
# =========================
# Correct EWMA recursion (equivalent to closed form with r0=0)
ewma_from_W <- function(W, lambda) {
  n <- length(W)
  r <- numeric(n)
  prev <- 0
  one_m <- 1 - lambda
  for (t in seq_len(n)) {
    prev <- lambda * W[t] + one_m * prev
    r[t] <- prev
  }
  r
}

# Control limits per your simplified spec (±L). If you later want exact Var[r_t] prefix, we can plug it in.
compute_limits <- function(L, n) {
  list(UCL = rep(L, n), LCL = rep(-L, n))
}

# One (λ, L) run returning a tidy frame for that combo
run_csb_ewma_one <- function(lambda, L, dfW) {
  n  <- nrow(dfW)
  rt <- ewma_from_W(dfW$W_t, lambda)
  lim <- compute_limits(L, n)
  sig_up    <- rt >  lim$UCL
  sig_down  <- rt <  lim$LCL
  sig_any   <- sig_up | sig_down
  first_idx <- if (any(sig_any)) which(sig_any)[1] else NA_integer_
  first_ts  <- if (!is.na(first_idx)) dfW$Timestamp[first_idx] else as.POSIXct(NA, tz = "UTC")
  
  tibble(
    Timestamp      = dfW$Timestamp,
    `Record number`= dfW$`Record number`,
    C_j            = dfW$C_j,
    Q_t            = dfW$Q_t,
    W_t            = dfW$W_t,
    r_t            = rt,
    UCL            = lim$UCL,
    LCL            = lim$LCL,
    signal_up      = sig_up,
    signal_down    = sig_down,
    signal_any     = sig_any,
    lambda         = lambda,
    L              = L,
    first_signal_index = first_idx,
    first_signal_time  = first_ts
  )
}

# Parameter grid
lambda_set <- c(0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
L_set      <- seq(1.5, 2.5, by = 0.1)
grid <- tidyr::expand_grid(lambda = lambda_set, L = L_set)

# =========================
# Step 5: Parallel compute (no workbook IO in workers)
# =========================
num_cores <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# Export only what workers need
parallel::clusterExport(
  cl,
  varlist = c("df_2024_w", "ewma_from_W", "compute_limits", "run_csb_ewma_one"),
  envir = environment()
)

t0 <- Sys.time()

res_list <- foreach(i = seq_len(nrow(grid)), .packages = c("dplyr","tibble","tidyr")) %dopar% {
  lam <- grid$lambda[i]
  LL  <- grid$L[i]
  run_csb_ewma_one(lam, LL, df_2024_w)
}

parallel::stopCluster(cl)

t1 <- Sys.time()
elapsed_minutes <- as.numeric(difftime(t1, t0, units = "mins"))
cat(sprintf("\n[Step 5] Parallel CSB-EWMA finished in %.2f minutes\n", elapsed_minutes))

# =========================
# Step 6: Write per-combo sheets and a summary
# =========================
wb <- createWorkbook()
for (i in seq_along(res_list)) {
  lam <- grid$lambda[i]
  LL  <- grid$L[i]
  sheet_name <- paste0("lambda_", lam, "_L_", LL)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, res_list[[i]])
}
saveWorkbook(wb, "csb_ewma_detection_results_by_combo.xlsx", overwrite = TRUE)
cat("Saved per-combo results to csb_ewma_detection_results_by_combo.xlsx\n")

summary_tbl <- dplyr::bind_rows(lapply(seq_along(res_list), function(i) {
  lam <- grid$lambda[i]
  LL  <- grid$L[i]
  df  <- res_list[[i]]
  tibble(
    lambda = lam,
    L = LL,
    n_obs = nrow(df),
    n_signals = sum(df$signal_any),
    first_signal_time = if (any(df$signal_any)) df$first_signal_time[which(df$signal_any)[1]] else as.POSIXct(NA, tz = "UTC")
  )
}))
readr::write_csv(summary_tbl, "csb_ewma_detection_summary.csv")
cat("Wrote summary to csb_ewma_detection_summary.csv\n")
cat(sprintf("Total elapsed time: %.2f minutes\n", elapsed_minutes))




















#########################################################


# r

# ==== prerequisites ====
library(dplyr)
library(tidyr)
library(lubridate)
library(readr)
library(ggplot2)
library(openxlsx)
library(doParallel)
library(foreach)

# ==== Step 1: data already cleaned/split and df_2024_binary exists ====
# Assumes df_2024_binary has: Timestamp, `Record number`,
# x_pH, x_DO_Sat, x_Turbidity, x_Salinity, x_Speed, x_Direction

# ---- Step 2: compute Ct, Qt, Wt with p0 from 2024 binary ----
df_2024_bin2 <- df_2024_binary %>%
  mutate(
    C_t = x_pH + x_DO_Sat + x_Turbidity + x_Salinity + x_Speed + x_Direction
  )

# p0_i from 2024 binary (proportion of 1s for each stream)
p0_vec <- colMeans(df_2024_bin2[, c("x_pH","x_DO_Sat","x_Turbidity",
                                    "x_Salinity","x_Speed","x_Direction")])
p0_total  <- sum(p0_vec)
sigma2_total <- sum(p0_vec * (1 - p0_vec))

# guardrails
stopifnot(is.finite(p0_total), is.finite(sigma2_total), sigma2_total > 0)

n <- nrow(df_2024_bin2)
t_idx <- seq_len(n)

Q_t <- cumsum(df_2024_bin2$C_t)
W_t <- (Q_t - p0_total * t_idx) / sqrt(t_idx * sigma2_total)

# quick sanity checks (these should look ~0 mean, ~1 sd, and not huge range)
cat("W_t mean=", mean(W_t), " sd=", sd(W_t),
    " min=", min(W_t), " max=", max(W_t), "\n")

df_2024_w <- df_2024_bin2 %>%
  mutate(W_t = W_t) %>%
  select(Timestamp, `Record number`, C_t, W_t)

# ---- Step 3: EWMA with steady-state limits ----
ewma_run <- function(W, lambda) {
  r <- numeric(length(W))
  prev <- 0
  one_m <- 1 - lambda
  for (i in seq_along(W)) {
    prev <- lambda * W[i] + one_m * prev
    r[i] <- prev
  }
  r
}

csb_ewma_once <- function(dfW, lambda, L) {
  r <- ewma_run(dfW$W_t, lambda)
  var_factor <- lambda / (2 - lambda)  # steady-state Var(r_t) when input Var(W)=1
  UCL <-  L * sqrt(var_factor)
  LCL <- -L * sqrt(var_factor)
  
  tibble::tibble(
    Timestamp     = dfW$Timestamp,
    RecordNumber  = dfW$`Record number`,
    W_t           = dfW$W_t,
    r_t           = r,
    UCL           = UCL,
    LCL           = LCL,
    lambda        = lambda,
    L             = L
  )
}

# ---- Step 4: run grid in parallel but write workbook sequentially ----
lambda_set <- c(0.1, 0.3, 0.4, 0.5, 0.6, 0.9)
L_set      <- seq(1.5, 2.5, by = 0.5)

num_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(num_cores)
registerDoParallel(cl)

start_time <- Sys.time()

# compute results in parallel and bring back as a list
res_list <- foreach(ll = lambda_set, .packages = c("dplyr","tibble")) %:%
  foreach(LL = L_set, .packages = c("dplyr","tibble")) %dopar% {
    csb_ewma_once(df_2024_w, ll, LL)
  }

stopCluster(cl)

# flatten to named list for writing
named_res <- list()
for (i in seq_along(lambda_set)) {
  for (j in seq_along(L_set)) {
    nm <- sprintf("lambda_%s_L_%s", format(lambda_set[i]), format(L_set[j]))
    named_res[[nm]] <- res_list[[i]][[j]]
  }
}

# write workbook sequentially (thread-safe)
wb <- createWorkbook()
for (nm in names(named_res)) {
  addWorksheet(wb, nm)
  writeData(wb, nm, named_res[[nm]])
}
saveWorkbook(wb, "csb_ewma_results_by_sheet.xlsx", overwrite = TRUE)

end_time <- Sys.time()
cat("Workbook written. Total minutes:", round(as.numeric(difftime(end_time, start_time, units = "mins")), 2), "\n")

# ---- Step 5: plot a few representative charts ----
plot_cc <- function(res_df) {
  ggplot(res_df, aes(Timestamp, r_t)) +
    geom_line(linewidth = 0.6) +
    geom_hline(yintercept = unique(res_df$UCL), linetype = "dashed", color = "red") +
    geom_hline(yintercept = unique(res_df$LCL), linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    labs(
      title = paste0("CSB-EWMA (λ=", res_df$lambda[1], ", L=", res_df$L[1], ")"),
      y = "EWMA statistic r_t", x = "Timestamp"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# example: plot two combinations
p1 <- plot_cc(named_res[["lambda_0.3_L_2"]])
p2 <- plot_cc(named_res[["lambda_0.6_L_2.5"]])
print(p1); print(p2)

# plot for lambda=0.9, L=1.5 
p3 <- plot_cc(named_res[["lambda_0.9_L_1.5"]])
p4 <- plot_cc(named_res[["lambda_0.1_L_2.5"]])
print(p3)
print(p4)





ggsave("plot_lambda_0.3_L_2.png", p1, width = 10, height = 5, dpi = 300)
ggsave("plot_lambda_0.6_L_2.5.png", p2, width = 10, height = 5, dpi = 300)
ggsave("plot_lambda_0.9_L_1.5.png", p3, width = 10, height = 5, dpi = 300)
ggsave("plot_lambda_0.1_L_2.5.png", p4, width = 10, height = 5, dpi = 300)





# =========================
# Step 7: Recompute r_t with steady-state variance adjustment and plot control charts
# =========================

library(ggplot2)

# Function to recompute EWMA and proper control limits
run_csb_ewma_varscaled <- function(lambda, L, dfW) {
  n  <- nrow(dfW)
  W  <- dfW$W_t
  r  <- numeric(n)
  prev <- 0
  one_m <- 1 - lambda
  for (t in seq_len(n)) {
    prev <- lambda * W[t] + one_m * prev
    r[t] <- prev
  }
  
  # steady-state EWMA variance factor
  var_factor <- lambda / (2 - lambda)
  UCL <-  L * sqrt(var_factor)
  LCL <- -L * sqrt(var_factor)
  
  tibble::tibble(
    Timestamp = dfW$Timestamp,
    RecordNumber = dfW$`Record number`,
    W_t = W,
    r_t = r,
    UCL = UCL,
    LCL = LCL,
    lambda = lambda,
    L = L
  )
}

# Choose a few representative parameter sets for visualization
plot_sets <- list(
  list(lambda = 0.3, L = 2.0),
  list(lambda = 0.3, L = 2.5),
  list(lambda = 0.6, L = 2.0),
  list(lambda = 0.6, L = 2.5)
)

# Compute and store the adjusted results
plot_results <- lapply(plot_sets, function(p) {
  run_csb_ewma_varscaled(p$lambda, p$L, df_2024_w)
})

# Plot function
plot_control_chart <- function(res) {
  ggplot(res, aes(x = Timestamp, y = r_t)) +
    geom_line(color = "steelblue", size = 0.6) +
    geom_hline(aes(yintercept = UCL), color = "red", linetype = "dashed") +
    geom_hline(aes(yintercept = LCL), color = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "black", size = 0.3) +
    labs(
      title = paste0("CSB-EWMA Control Chart (λ=", res$lambda[1],
                     ", L=", res$L[1], ")"),
      y = "EWMA Statistic (r_t)",
      x = "Timestamp"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Display each chart
for (res in plot_results) {
  print(plot_control_chart(res))
}

# Optionally, save them as PNGs
for (res in plot_results) {
  fname <- paste0("control_chart_lambda_", res$lambda[1], "_L_", res$L[1], ".png")
  ggsave(filename = fname, plot = plot_control_chart(res),
         width = 10, height = 5, dpi = 300)
}
cat("\n[Step 7] Control charts with variance-adjusted limits saved as PNGs.\n")