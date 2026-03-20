# EIT-Metabolic-Syndrome

**Metabolic Syndrome and Regional Lung Physiology in ARDS: A Bayesian Analysis of Electrical Impedance Tomography Data**

T. Gaulton, MD MSc — Massachusetts General Hospital / Harvard Medical School

---

## Overview

This repository contains the analysis code for a study examining how metabolic syndrome affects regional ventilation-perfusion distributions in patients with acute respiratory distress syndrome (ARDS), measured by electrical impedance tomography (EIT). The analysis uses Bayesian regression to evaluate whether metabolic syndrome modifies regional lung physiology through mechanical and/or vascular pathways.

**Key questions:**

- Does metabolic syndrome reduce dorsal (dependent) ventilation?
- Does metabolic syndrome impair gravitational perfusion redistribution?
- Do combined ventilation and perfusion abnormalities increase dorsal shunt?

---

## Repository Structure

```
EIT-Metabolic-Syndrome/
├── analysis_main.R            # Complete analysis script
├── data/                      # Input data (not included; see Data Availability)
├── models/
│   ├── bayesian/              # Primary brms model objects (auto-created on run)
│   ├── sensitivity/           # Interaction model
│   └── gradients/             # Gradient and decomposition models
├── results/
│   ├── tables/                # Table 1 and Table 2 (CSV)
│   └── data/                  # Model summaries, diagnostics, complete workspace
├── figures/                   # All figures (PDF and PNG)
└── submission_figures/        # Final submission-ready figures
```

---

## Data Availability

Patient-level data are not included in this repository. The analysis requires two Excel files placed in the `data/` directory:

- `data/Phenotypes_distributions_demo.xlsx` — demographic and clinical variables
- `data/Phenotypes_distributions_vq.xlsx` — regional EIT ventilation-perfusion data

The datasets can be made available upon reasonable request to the corresponding author, subject to data use agreement and IRB requirements.

---

## Requirements

**R version:** 4.3.0 or later recommended

**Packages:**

| Package    | Purpose                          |
|------------|----------------------------------|
| tidyverse  | Data wrangling and visualization |
| brms       | Bayesian multilevel regression   |
| bayesplot  | MCMC diagnostics                 |
| patchwork  | Figure composition               |
| readxl     | Reading input Excel files        |

Install all packages with:

```r
install.packages(c("tidyverse", "brms", "bayesplot", "patchwork", "readxl"))
```

Stan (via `rstan` or `cmdstanr`) is required as a backend for `brms`. See the [brms installation guide](https://github.com/paul-buerkner/brms#how-to-use-brms) for details.

---

## Usage

1. Clone the repository:

```bash
git clone https://github.com/[username]/EIT-Metabolic-Syndrome.git
cd EIT-Metabolic-Syndrome
```

2. Place data files in the `data/` directory (see Data Availability above).

3. Open `EIT-Metabolic-Syndrome.Rproj` in RStudio, then run:

```r
source("analysis_main.R")
```

The script will create all output directories, fit Bayesian models (expect 10–20 minutes on first run), and save all results, figures, and tables. Fitted models are cached — subsequent runs are fast unless the model specification changes.

To reload results without refitting models:

```r
load("results/data/complete_analysis.RData")
```

---

## Session Information

A `session_info.txt` file is written to `results/data/` each time the script runs, recording the R version and package versions used.

---

## Citation

If you use this code, please cite the associated manuscript (citation to be added upon publication).

---

## License

MIT License. See `LICENSE` for details.

---

## Contact

Timothy Gaulton, MD MSc
Division of Critical Care, Department of Anesthesia, Critical Care and Pain Medicine
Massachusetts General Hospital / Harvard Medical School
