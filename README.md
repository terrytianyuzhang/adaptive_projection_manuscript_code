# Adaptive Projection Manuscript Code

This repository contains the code used in the paper titled **"Adaptive Projected Two-Sample Comparisons for Single-Cell Expression Data"**, currently under review. The code implements statistical procedures and generates all the figures included in the manuscript.

## 📁 Repository Structure

```
├── README.md
└── code_submitted/
    ├── approximate_orthogonality/     # Simulation for Approximate Orthogonality section
    ├── Biostatistics_revision/       # Sparse-covariance simulations for the Biostatistics revision
    ├── cleary_data_mean_comparison/   # Perturb-seq data analysis [1]
    ├── code_paper/                    # Core functions and utilities
    ├── jinhong_deviance/              # Application to a Lupus study [2]
    ├── main_simulation/               # Type-I error and power assessment
    └── try_Cleary_data/               # Preprocessed data for Perturb-seq analysis [1]
```


## References

[1] Yao, D., Binan, L., Bezney, J., Simonton, B., Freedman, J., Frangieh, C. J., Dey, K., Geiger-Schuller, K., Eraslan, B., Gusev, A., et al. (2024). Scalable genetic screening for regulatory circuits using compressed Perturb-seq. *Nature Biotechnology*, 42(8), 1282–1295.

[2] Perez, R. K., Gordon, M. G., Subramaniam, M., Kim, M. C., Hartoularos, G. C., Targ, S., Sun, Y., Ogorodnikov, A., Bueno, R., Lu, A., et al. (2022). Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. *Science*, 376(6589), eabf1970.
