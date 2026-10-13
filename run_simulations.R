##############################################################################
# Replication Script: Monte Carlo Simulations and Mosaic Plot Generation
# Paper: "ggskewboxplots: An R Package for Outlier Detection in Skewed Data"
# Author: Mustafa Cavus
#
# Description:
#   This self-contained script reproduces the full Monte Carlo simulation study
#   and generates all 84 mosaic plot artifacts (21 swamping + 63 masking) using
#   the Skewed Exponential Power Distribution (SEPD) across a 49x49 grid.
#
# Requirements:
#   - R (>= 4.0.0)
#   - Packages: ggskewboxplots, viridisLite, parallel
##############################################################################

# ============================================================================
# 0. Dependencies and Parallel Configuration
# ============================================================================

required_packages <- c("ggskewboxplots", "viridisLite")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(ggskewboxplots)
library(viridisLite)
library(parallel)

# Detect CPU cores (Linux/macOS support mclapply; on Windows fallback to 1)
is_windows <- .Platform$OS.type == "windows"
N_CORES <- if (is_windows) 1 else max(1, parallel::detectCores() - 1)
cat(sprintf("[Setup] Detected %d available cores. Using %d parallel worker(s).\n", 
            parallel::detectCores(), N_CORES))

# Output directories
dir.create("figures_png", showWarnings = FALSE)
dir.create("simulation_rds", showWarnings = FALSE)

# ============================================================================
# 1. Skewed Exponential Power Distribution (SEPD) Random Generator
# ============================================================================

rsepd <- function(n, mu = 0, sigma = 1, alpha = 0.5, p = 2) {
  kp <- 1 / (2 * gamma(1 + 1/p) * p^(1/p))
  u  <- runif(n)
  w  <- rgamma(n, shape = 1/p, scale = 1)
  y  <- ifelse(u > alpha, 1 - alpha, -alpha) * (w^(1/p)) / gamma(1 + 1/p)
  y  <- y / kp
  A  <- 2 * p^(1/p) * ((1 - alpha)^2 - alpha^2) * gamma(2/p) / gamma(1/p)
  A2 <- (2 * p^(1/p))^2 * ((1 - alpha)^3 + alpha^3) * gamma(3/p) / gamma(1/p)
  B  <- sqrt(A2 - A^2)
  mu + sigma * (y - A) / B
}

# ============================================================================
# 2. Performance Metrics (Swamping & Masking)
# ============================================================================

# Swamping: Uncontaminated sample, proportion of false positives outside [ymin, ymax]
mc_swamping <- function(siz, nsim, mu, sigma, alpha, p, method, k = 1.5) {
  ysim <- matrix(
    rsepd(siz * nsim, mu = mu, sigma = sigma, alpha = alpha, p = p),
    nrow = nsim, ncol = siz
  )
  res <- apply(ysim, 1, function(x) {
    st <- compute_skew_stats(x, method = method, k = k)
    (sum(x < st$ymin) + sum(x > st$ymax)) / length(x)
  })
  mean(res)
}

# Masking: Synthetic contamination via order statistics multiplied by c_mult
# Masking Error = 1 - detection_rate (proportion of undetected genuine outliers)
mc_masking <- function(siz, nsim, mu, sigma, alpha, p, method,
                       contam_rate = 0.05, c_mult = 2, k = 1.5) {
  n_each   <- max(1, floor(siz * contam_rate / 2))
  detected <- numeric(nsim)
  
  for (m in seq_len(nsim)) {
    x <- sort(rsepd(siz, mu = mu, sigma = sigma, alpha = alpha, p = p))
    # Contaminate lower and upper extreme tails
    x[seq_len(n_each)]         <- x[seq_len(n_each)] * c_mult
    x[(siz - n_each + 1):siz] <- x[(siz - n_each + 1):siz] * c_mult
    contam_idx <- c(seq_len(n_each), (siz - n_each + 1):siz)
    
    st <- tryCatch(
      compute_skew_stats(x, method = method, k = k),
      error = function(e) list(ymin = -Inf, ymax = Inf)
    )
    detected[m] <- sum(
      x[contam_idx] < st$ymin | x[contam_idx] > st$ymax
    ) / (2 * n_each)
  }
  mean(detected)
}

# ============================================================================
# 3. Parameter Grid and Mosaic Plot Visualizer
# ============================================================================

make_grid <- function(n = 49) {
  if (n %% 2 == 0) n <- n + 1
  alfa <- (0:(n - 1)) / (n - 1)
  alfa[alfa == 0] <- 1e-12
  alfa[alfa == 1] <- 1 - 1e-12
  j      <- (log10(50) / log10(2))^(1 / (n / 2 - 0.5))
  p_grid <- c(1, 2^(1 / (j^(((n - 3) / 2):1))), 2, 2^(j^((1:(n / 2 - 0.5)))))
  expand.grid(alpha = alfa, p = p_grid)
}

mosaic_draw <- function(scores, n = 49, label = FALSE) {
  if (n %% 2 == 0) n <- n + 1
  mu <- 0; sigma <- 1; s <- 3.5
  ymax_plot <- c(0, 1 / sigma)
  
  alfa <- (0:(n - 1)) / (n - 1)
  alfa[alfa == 0] <- 1e-12
  alfa[alfa == 1] <- 1 - 1e-12
  j      <- (log10(50) / log10(2))^(1 / (n / 2 - 0.5))
  p_grid <- c(1, 2^(1 / (j^(((n - 3) / 2):1))), 2, 2^(j^((1:(n / 2 - 0.5)))))
  params <- expand.grid(alpha = alfa, p = p_grid)
  
  cex1 <- c(1.2, 1, 0.8, 0.7, 0.6)
  old.par <- par(
    mfrow = c(n, n), xaxt = "n", yaxt = "n",
    mar = c(0, 0, 0, 0),
    cex = cex1[findInterval(n, c(1, 5, 10, 15, 25, 9999))],
    xaxs = "i", yaxs = "i"
  )
  
  palette <- viridisLite::viridis(11, direction = 1)
  breaks  <- c(0, 0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065, 0.075, 0.085, 0.095, 1)
  
  for (i in seq_len(nrow(params))) {
    a  <- params$alpha[i]
    pv <- params$p[i]
    sc <- scores[i]
    plot(0, 0, xlim = c(mu - s * sigma, mu + s * sigma),
         ylim = ymax_plot, type = "n", axes = FALSE)
    id   <- findInterval(sc, breaks, all.inside = TRUE)
    colo <- palette[id]
    rect(mu - s * sigma, 0, mu + s * sigma, ymax_plot[2], col = colo, border = NA)
    if (a == 0.5 && pv == 2) box(lwd = 2, col = "white") else box()
    if (label) text(mu, ymax_plot[2] * 0.5, sc,
                    col = ifelse(id <= 5, "white", "black"))
  }
  par(old.par)
}

# ============================================================================
# 4. Simulation Execution Loop
# ============================================================================

N_GRID      <- 49
NSIM        <- 10000
SIZES       <- c(20, 50, 100)
CONTAM_RATE <- 0.05
C_LEVELS    <- c(2, 4, 6)

methods <- c("tukey", "kimber", "hubert", "adil", "babura", "walker", "junsawang")

# File naming alignment with the paper manuscript
paper_names <- list(
  tukey     = list(swamp = "swamping_tukey",     mask = "Classic"),
  kimber    = list(swamp = "swamping_kimber",    mask = "Kimber"),
  hubert    = list(swamp = "swamping_hubert",    mask = "Hubert"),
  adil      = list(swamp = "swamping_adil",      mask = "Adil"),
  babura    = list(swamp = "swamping_babura",    mask = "Babura"),
  walker    = list(swamp = "swamping_walker",    mask = "Walker"),
  junsawang = list(swamp = "swamping_junsawang", mask = "Jusawang")
)

params  <- make_grid(N_GRID)
n_cells <- nrow(params)
total_steps <- length(methods) * length(SIZES) * (1 + length(C_LEVELS))
current_step <- 0

start_all <- Sys.time()
cat(sprintf("[Simulation Start] Total %d parameter cells per mosaic | %d iterations\n", n_cells, NSIM))

for (m in methods) {
  cat(sprintf("\n>>> Method: %s <<<\n", m))
  
  for (siz in SIZES) {
    cat(sprintf("  -> Sample size n = %d\n", siz))
    
    # ----------------------------------------------------
    # A. Swamping Simulation
    # ----------------------------------------------------
    t0 <- Sys.time()
    cat("     [1/4] Computing Swamping...")
    res_sw <- if (is_windows) {
      lapply(seq_len(n_cells), function(i) {
        mc_swamping(siz, NSIM, 0, 1, params$alpha[i], params$p[i], method = m)
      })
    } else {
      mclapply(seq_len(n_cells), function(i) {
        mc_swamping(siz, NSIM, 0, 1, params$alpha[i], params$p[i], method = m)
      }, mc.cores = N_CORES)
    }
    
    scores_sw <- unlist(res_sw)
    saveRDS(scores_sw, file = sprintf("simulation_rds/swamping_%s_n%d.rds", m, siz))
    
    png_sw <- sprintf("figures_png/%s_n%d.png", paper_names[[m]]$swamp, siz)
    png(png_sw, width = 500, height = 500, res = 96)
    mosaic_draw(scores_sw, n = N_GRID, label = FALSE)
    dev.off()
    
    current_step <- current_step + 1
    cat(sprintf(" Done (%.1f min) -> Saved: %s [%d/%d]\n",
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                png_sw, current_step, total_steps))
    
    # ----------------------------------------------------
    # B. Masking Simulation across c in {2, 4, 6}
    # ----------------------------------------------------
    for (cv in C_LEVELS) {
      t0 <- Sys.time()
      cat(sprintf("     Computing Masking (c = %d)...", cv))
      res_m <- if (is_windows) {
        lapply(seq_len(n_cells), function(i) {
          det <- mc_masking(siz, NSIM, 0, 1, params$alpha[i], params$p[i],
                            method = m, contam_rate = CONTAM_RATE, c_mult = cv)
          1 - det
        })
      } else {
        mclapply(seq_len(n_cells), function(i) {
          det <- mc_masking(siz, NSIM, 0, 1, params$alpha[i], params$p[i],
                            method = m, contam_rate = CONTAM_RATE, c_mult = cv)
          1 - det
        }, mc.cores = N_CORES)
      }
      
      scores_m <- unlist(res_m)
      saveRDS(scores_m, file = sprintf("simulation_rds/masking_c%d_%s_n%d.rds", cv, m, siz))
      
      png_m <- sprintf("figures_png/masking_c%d_%s_n%d.png", cv, paper_names[[m]]$mask, siz)
      png(png_m, width = 500, height = 500, res = 96)
      mosaic_draw(scores_m, n = N_GRID, label = FALSE)
      dev.off()
      
      current_step <- current_step + 1
      cat(sprintf(" Done (%.1f min) -> Saved: %s [%d/%d]\n",
                  as.numeric(difftime(Sys.time(), t0, units = "mins")),
                  png_m, current_step, total_steps))
    }
  }
}

cat(sprintf("\n======================================================\n"))
cat(sprintf("[Completed] Total time: %.2f hours.\n", 
            as.numeric(difftime(Sys.time(), start_all, units = "hours"))))
cat(sprintf("Generated 84 mosaic figures in 'figures_png/' and RDS data in 'simulation_rds/'.\n"))
