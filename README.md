# ggskewboxplots: Enhanced Boxplots for Skewed Data in R

### Replication Materials for Simulation Studies

**Author:** Mustafa Cavus  
**Software Package:** `ggskewboxplots` (CRAN: https://cran.r-project.org/package=ggskewboxplots)

---

## 1. Overview

This directory provides the complete, self-contained replication workflow required to reproduce all numerical results, Monte Carlo simulations, and mosaic figures presented in the paper:
- **Swamping Evaluation:** False positive rates across uncontaminated distributions (Figure 6).
- **Masking Evaluation (Mild, $c = 2$):** Outlier masking errors under mild contamination (Figure 7).
- **Contamination Severity Evaluation ($c \in \{2, 4, 6\}$):** Masking errors under varying outlier extremity for sample sizes $n = 20$ (Figure 8), $n = 50$ (Appendix Figure A1), and $n = 100$ (Appendix Figure A2).

A total of **84 unique mosaic plots** ($7 \text{ methods} \times 3 \text{ sample sizes} \times 4 \text{ scenarios}$) are produced.

---

## 2. Directory Structure

```text
replication_package/
├── README.md               # Complete replication documentation (this file)
├── run_simulations.R       # Standalone Monte Carlo simulation script (generates RDS & PNGs)
├── render_from_rds.R       # Instant plotting script from precomputed RDS score vectors
├── simulation_rds/         # Directory storing raw Monte Carlo error rates (*.rds)
└── figures_png/            # Generated high-resolution mosaic PNG files (84 figures)
```

---

## 3. Computational Requirements & Dependencies

### System Requirements
- **R Version:** $\ge 4.0.0$
- **Operating System:** macOS, Linux, or Windows (parallel execution uses `parallel::mclapply` on Unix-like systems and single-core fallback on Windows).
- **Recommended Hardware:** Multi-core processor (4–8+ cores) and at least 8 GB of RAM.

### Required R Packages
Install the required CRAN packages prior to running the scripts:

```r
install.packages(c("ggskewboxplots", "viridisLite", "parallel"))
```

---

## 4. Replication Instructions

### Option A: Complete Monte Carlo Simulation from Scratch (Full Run)
To re-run the entire Monte Carlo experiment ($10{,}000$ replications per cell over a $49 \times 49 = 2{,}401$ grid of skewness $\alpha$ and kurtosis $p$ parameters):

```bash
# Navigate to the replication directory
cd /Users/mustafacavus/Documents/ggskewboxplots/replication_package

# Execute the simulation script in background or terminal
Rscript run_simulations.R
```

**Runtime Expectation:**
- Depending on the number of available CPU cores (typically 7–8 cores), the full simulation across all 7 methods, 3 sample sizes, and 4 scenarios takes approximately **12 to 18 hours**.
- As each step finishes, the raw numeric matrix is saved to `simulation_rds/` and the resulting mosaic PNG is saved to `figures_png/`.

---

### Option B: Fast Reproduction / Plot Re-rendering (Using Cached Data)
If you already have the precomputed `.rds` files generated or copied from the repository (`../codes/rds_scores/`), you can instantly regenerate all 84 publication-quality mosaic graphics without re-running the heavy simulations:

```bash
cd /Users/mustafacavus/Documents/ggskewboxplots/replication_package

Rscript render_from_rds.R
```

**Runtime Expectation:** ~5 to 10 seconds.

---

## 5. Parameter Setup & Method Specifications

| Parameter | Values / Descriptions |
| :--- | :--- |
| **Methods Evaluated ($7$)** | Tukey (`classic`), Kimber (`kimber`), Hubert (`hubert`), Adil (`adil`), Babura (`babura`), Walker (`walker`), Junsawang (`junsawang`) |
| **Sample Sizes ($n$)** | $n \in \{20, 50, 100\}$ |
| **Mosaic Grid Size** | $49 \times 49$ parameter pairs ($\alpha \in (0, 1)$, $p \in (0.5, 50)$) |
| **Monte Carlo Replications** | $10{,}000$ replications per parameter cell |
| **Contamination Fraction** | $\varepsilon = 0.05$ ($k = \max(1, \lfloor n \varepsilon / 2 \rfloor)$ in each tail) |
| **Contamination Extremity ($c$)**| Mild ($c = 2$), Moderate ($c = 4$), Severe ($c = 6$) |
| **Color Scale** | Viridis palette ($0\%$ dark purple to $\ge 10\%$ bright yellow) |

---

## 6. Verification and Matching with the Manuscript

Generated files in `figures_png/` map directly to the figures in the LaTeX manuscript (`template.tex`):

- **Figure 6 (Swamping):** `swamping_{method}_n{20,50,100}.png`
- **Figure 7 (Masking $c=2$):** `masking_c2_{Method}_n{20,50,100}.png`
- **Figure 8 (Masking $n=20$ by $c$):** `masking_c{2,4,6}_{Method}_n20.png`
- **Appendix Figure A1 ($n=50$):** `masking_c{2,4,6}_{Method}_n50.png`
- **Appendix Figure A2 ($n=100$):** `masking_c{2,4,6}_{Method}_n100.png`

---

## 7. Citation

If using this replication package or the `ggskewboxplots` package, please cite:

```bibtex
@article{cavus2025ggskewboxplots,
  title     = {ggskewboxplots: Enhanced Boxplots for Skewed Data in R},
  author    = {Cavus, Mustafa},
  year      = {2025}
}
```
