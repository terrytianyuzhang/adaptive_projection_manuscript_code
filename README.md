# Adaptive Projection Manuscript Code

This repository contains the code used in the paper titled **"[Insert full paper title here]"**, currently under review. The code implements statistical procedures and generates all the figures and tables included in the manuscript.

## 📁 Repository Structure

```
.
├── data/                            # Processed or raw data (if shareable)
├── results/                         # Output folder for generated plots/results
├── scripts/
│   ├── fig1_generate_projection.R   # Code for Figure 1
│   ├── fig2_power_comparison.R      # Code for Figure 2
│   ├── fig3_real_data_analysis.R    # Code for Figure 3
│   ├── ...
│   └── helper_functions.R           # Utility functions used across scripts
├── install_packages.R               # Installs required R packages
├── README.md                        # This file
└── session_info.txt                 # R session info for reproducibility
```

## 🛠️ Setup Instructions

1. **Clone the repository**:

```bash
git clone https://github.com/your_username/adaptive_projection_manuscript_code.git
cd adaptive_projection_manuscript_code
```

2. **Install required R packages**:

From within R:

```r
source("install_packages.R")
```

This script installs all necessary R packages (e.g., `glmnet`, `ggplot2`, `data.table`, `Seurat`, etc.).

3. **Run analysis scripts**:

Each script in `scripts/` is responsible for generating specific figures in the paper. Run them sequentially to reproduce the results.

| Figure | Script                          | Description                                |
|--------|----------------------------------|--------------------------------------------|
| Fig 1  | `scripts/fig1_generate_projection.R` | Synthetic example for projection directions |
| Fig 2  | `scripts/fig2_power_comparison.R`     | Power comparison under various alternatives |
| Fig 3  | `scripts/fig3_real_data_analysis.R`  | Real data analysis (e.g., CommonMind)       |
| ...    | ...                              | ...                                        |

All plots and intermediate results will be saved in the `results/` folder.

## 🧬 Data

- Some datasets used in this manuscript (e.g., **CommonMind**) are large or subject to sharing restrictions.
- Where possible, we provide:
  - Preprocessed versions (in `data/`)
  - Instructions or links for public access
- For restricted datasets, please follow the instructions in the manuscript or contact the data owners.

## 📜 Reproducibility

- A snapshot of the R session (`session_info.txt`) is provided for reproducibility.
- Random seeds are fixed within each script to ensure consistent results.

## 📧 Contact

For questions or issues, please contact:

**Tianyu Zhang**  
Department of Statistics and Applied Probability  
University of California, Santa Barbara  
📧 your.email@ucsb.edu
