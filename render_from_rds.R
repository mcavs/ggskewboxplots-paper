##############################################################################
# Quick Visualization Re-rendering Script from Cached RDS Scores
# Paper: "ggskewboxplots: An R Package for Outlier Detection in Skewed Data"
#
# Description:
#   If the simulation scores (*.rds) have already been generated, this script
#   re-draws all 84 mosaic PNG plots instantly (in ~5 seconds) without
#   running any computationally intensive Monte Carlo iterations.
##############################################################################

if (!requireNamespace("viridisLite", quietly = TRUE)) {
  install.packages("viridisLite", repos = "https://cloud.r-project.org")
}

library(viridisLite)

dir.create("figures_png", showWarnings = FALSE)

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

methods <- c("tukey", "kimber", "hubert", "adil", "babura", "walker", "junsawang")
SIZES   <- c(20, 50, 100)
C_LEVELS <- c(2, 4, 6)

paper_names <- list(
  tukey     = list(swamp = "swamping_tukey",     mask = "Classic"),
  kimber    = list(swamp = "swamping_kimber",    mask = "Kimber"),
  hubert    = list(swamp = "swamping_hubert",    mask = "Hubert"),
  adil      = list(swamp = "swamping_adil",      mask = "Adil"),
  babura    = list(swamp = "swamping_babura",    mask = "Babura"),
  walker    = list(swamp = "swamping_walker",    mask = "Walker"),
  junsawang = list(swamp = "swamping_junsawang", mask = "Jusawang")
)

rds_dir <- if (dir.exists("simulation_rds")) "simulation_rds" else "../codes/rds_scores"

cat(sprintf("[Render] Reading cached scores from: %s\n", rds_dir))

for (m in methods) {
  for (siz in SIZES) {
    # Swamping
    f_sw_rds <- file.path(rds_dir, sprintf("swamping_%s_n%d.rds", m, siz))
    if (file.exists(f_sw_rds)) {
      scores <- readRDS(f_sw_rds)
      png_out <- sprintf("figures_png/%s_n%d.png", paper_names[[m]]$swamp, siz)
      png(png_out, width = 500, height = 500, res = 96)
      mosaic_draw(scores, n = 49, label = FALSE)
      dev.off()
    }
    
    # Masking
    for (cv in C_LEVELS) {
      f_m_rds <- file.path(rds_dir, sprintf("masking_c%d_%s_n%d.rds", cv, m, siz))
      if (file.exists(f_m_rds)) {
        scores <- readRDS(f_m_rds)
        png_out <- sprintf("figures_png/masking_c%d_%s_n%d.png", cv, paper_names[[m]]$mask, siz)
        png(png_out, width = 500, height = 500, res = 96)
        mosaic_draw(scores, n = 49, label = FALSE)
        dev.off()
      }
    }
  }
}

cat("[Render Complete] All available mosaic plots regenerated in 'figures_png/'.\n")
