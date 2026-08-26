# TLSLeafArea-Quercus

R scripts for whole-plant leaf area and Leaf Area Index (LAI) estimation from terrestrial laser scanning (TLS) point clouds, using Random Forest wood-leaf classification and the Maximum Unit Normal Component (MUNC) area estimation method.

Developed and tested on *Quercus frainetto* Ten. individuals in a declining Mediterranean oak stand in Basilicata (southern Italy).

Associated publication:
> Giuzio E (2026). A reproducible terrestrial laser scanning workflow for whole-plant leaf area and leaf area index estimation in declining Mediterranean oaks. *iForest – Biogeosciences and Forestry*. ms26/5291.

---

## Workflow overview

```
TLS point clouds (.las)
        │
        ▼
rf_classify.R  ──►  Wood-leaf separation (Random Forest)
                     Output: *_RF.las (RF.class, RF.prob fields)
        │
        ▼
rf_report.R    ──►  Leaf area estimation (MUNC method)
                     Uncertainty quantification (FP, FN, σ_geo)
                     Output: summary table, diagnostic plots
```

---

## Scripts

### `rf_classify.R`
Trains a Random Forest classifier to separate leaf and wood points in TLS point clouds.

**Input:**
- `wood_training.las` — manually labelled wood points
- `leaf_training.las` — manually labelled leaf points
- Individual tree `.las` files to classify

**Features used (26 per point):**
- Normalised height (Z)
- RGB colour components (R, G, B)
- Excess Green Index (ExG)
- Multiscale geometric features (linearity, planarity, sphericity, curvature, point density) at four radii: 0.05, 0.10, 0.20, 0.40 m

**Output per tree:**
- `*_RF.las` with added fields:
  - `RF.class`: binary classification (1 = leaf, 2 = wood)
  - `RF.prob`: posterior probability of belonging to the leaf class

**Output batch:**
- `00_RF_summary.csv`
- `00_RF_feature_importance.png`
- `00_RF_vs_CANUPO.png` (if CANUPO classification is available)

---

### `rf_report.R`
Estimates whole-plant leaf area and LAI from classified point clouds using the MUNC method, and quantifies estimation uncertainty.

**Input:**
- `*_RF.las` files produced by `rf_classify.R`
- Co-registration error values (from Cyclone REGISTER 360 PLUS reports, entered manually)
- Projected crown area (PCA) per tree (from ITSMe R package)

**Uncertainty quantification:**
- Classification uncertainty: asymmetric range [−FP, +FN] based on posterior RF.prob
- Geometric uncertainty: σ_geo propagated in quadrature from co-registration error, voxelisation positional uncertainty, and normal estimation error

**Output per tree:**
- `*_report.txt`: leaf area, LAI, uncertainty components
- Diagnostic plots (top view, side view, vertical profile)

**Output batch:**
- `00_RF_summary.csv`: complete results table
- `00_compare_area_RF.png`: leaf area comparison across individuals

---

## Requirements

- R version [X.X.X]
- Required packages:

| Package  | Purpose                        |
|----------|--------------------------------|
| lidR     | LAS file I/O and processing    |
| ranger   | Random Forest implementation   |
| caret    | Cross-validation               |
| FNN      | k-nearest neighbour search     |
| ggplot2  | Diagnostic plots               |
| dplyr    | Data manipulation              |
| tidyr    | Data reshaping                 |
| viridis  | Colour scales                  |
| progress | Progress bar                   |
| tools    | File path utilities            |

Install all packages with:
```r
install.packages(c("lidR", "ranger", "caret", "FNN", "ggplot2",
                   "dplyr", "tidyr", "viridis", "progress", "tools"))
```

---

## Usage

1. Edit the `params` list at the top of each script to set input/output paths.
2. Run `rf_classify.R` to train the RF classifier and classify all trees.
3. Run `rf_report.R` to estimate leaf area and quantify uncertainty.

> **Note:** Before applying to a new species or scanning season, the RF classifier should be retrained or re-calibrated on representative training data, particularly where the spectral contrast between leaf and woody tissue differs from that of deciduous oaks during the leaf-on season.

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
