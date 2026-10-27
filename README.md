# TLSLeafArea-Quercus

R pipeline for whole-plant leaf area and Leaf Area Index (LAI) estimation from terrestrial laser scanning (TLS) point clouds, using Random Forest wood–leaf classification and the Maximum Unit Normal Component (MUNC) area estimation method.

Developed and tested on *Quercus frainetto* Ten. individuals in a declining Mediterranean oak stand in Basilicata (southern Italy).

**Associated publication:**

Giuzio E (2026). A reproducible terrestrial laser scanning workflow for whole-plant leaf area and leaf area index estimation in declining Mediterranean oaks. *iForest – Biogeosciences and Forestry*. ms26/5291.

---

## Workflow overview

```
TLS point clouds (.las)
        │
        ▼
[Section 3] RF classification ──►  Wood–leaf separation (Random Forest)
                                    Output: *_RF.las  (RF.class, RF.prob)
        │
        ▼
[Section 4] LOTO validation   ──►  Leave-One-Tree-Out cross-validation
                                    Output: LOTO_results.csv, metrics plot
        │
        ▼
[Section 5] Leaf area         ──►  MUNC estimation + uncertainty
                                    Output: RF_summary.csv, diagnostic plots
        │
        ▼
[Section 6] LAI               ──►  Crown-projected LAI via ITSMe
                                    Output: LAI_summary.csv, comparison plots
```

All steps are contained in the single script **`tls_leaf_area_pipeline.R`**.  
Run them in order by calling the corresponding batch functions (see Section 7 of the script for the exact call sequence).

---

## Script structure

| Section | Content |
|---------|---------|
| 0 | Package installation and loading |
| 1 | `params` list — **edit paths and co-registration errors here** |
| 2 | Helper functions (feature extraction, MUNC, error propagation) |
| 3 | RF training, classification, tiled classification for large clouds |
| 4 | Leave-One-Tree-Out (LOTO) cross-validation |
| 5 | Leaf area estimation (MUNC) and uncertainty quantification |
| 6 | Crown-projected LAI (ITSMe) |
| 7 | Step-by-step run instructions |

---

## Features used for classification (26 per point)

- Normalised height (Z)
- RGB colour components (R, G, B) and Excess Green Index (ExG)
- Laser intensity, return number, number of returns, return ratio
- Multiscale geometric features (linearity, planarity, sphericity, curvature, local point density) at four radii: 0.05, 0.10, 0.20, 0.40 m

---

## Uncertainty quantification

**Classification uncertainty** — asymmetric range [−FP, +FN] derived from posterior RF probabilities (RF.prob): false positives overestimate leaf area, false negatives underestimate it.

**Geometric uncertainty** — σ_geo propagated in quadrature from three components:
1. Co-registration error (from Cyclone REGISTER 360 PLUS reports, entered manually in `params`)
2. Voxelisation positional uncertainty
3. Normal estimation error (PCA-based)

Final credible range per tree: [area − FP − 1.96 σ_geo, area + FN + 1.96 σ_geo]

---

## Requirements

R ≥ 4.1.0

| Package | Purpose |
|---------|---------|
| lidR | LAS file I/O and processing |
| ranger | Random Forest implementation |
| caret | Cross-validation |
| FNN | k-nearest neighbour search |
| ggplot2 | Diagnostic plots |
| dplyr | Data manipulation |
| tidyr | Data reshaping |
| viridis | Colour scales |
| progress | Progress bar |
| tools | File path utilities |
| ITSMe | Crown projected area (GitHub only) |

Install CRAN packages:

```r
install.packages(c("lidR", "ranger", "caret", "FNN", "ggplot2",
                   "dplyr", "tidyr", "viridis", "progress", "tools"))
```

Install ITSMe from GitHub:

```r
# install.packages("devtools")
devtools::install_github("lmterryn/ITSMe")
```

---

## Usage

1. Open `tls_leaf_area_pipeline.R` and edit the `params` list in **Section 1**:
   - `input_dir` / `output_dir` — working directories
   - `training_wood_las` / `training_leaf_las` — manually labelled training files
   - `coreg_table` — co-registration errors from your scan registration report
2. Run **Section 3** (`run_classification_batch()`) to train the RF and classify all trees.
3. *(Optional)* Run **Section 4** (`run_loto()`) to validate the classifier with Leave-One-Tree-Out CV. Requires annotated wood/leaf files for at least two trees.
4. Run **Section 5** (`run_leaf_area_batch()`) to estimate leaf area and propagate uncertainties.
5. Run **Section 6** (`run_lai_batch()`) to compute crown-projected LAI.

> **Note:** Before applying to a new species or scanning season, retrain the RF classifier on representative training data, particularly where the spectral contrast between leaf and woody tissue differs from that of deciduous oaks during the leaf-on season.

---

## Data

Point cloud data and training files are not included in this repository due to file size constraints. The trained RF model (`modello_RF.rds`) is available upon request to the corresponding author (emanuele.giuzio@unibas.it).

---

## License

This project is licensed under the MIT License.

---

## Contact

Emanuele Giuzio  
University of Basilicata, Department of Agricultural, Forestry, Food and Environmental Sciences (DAFE)  
emanuele.giuzio@unibas.it
