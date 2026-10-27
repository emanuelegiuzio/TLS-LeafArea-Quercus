# TLS LEAF AREA PIPELINE
# Random Forest wood/leaf classification, Leave-One-Tree-Out (LOTO) cross-validation,
# MUNC leaf area estimation with error propagation, and LAI calculation
# for TLS point clouds acquired with Leica BLK360.
#
# Workflow:
#   1. Train RF classifier on annotated wood/leaf point clouds
#   2. Classify all trees (batch); optionally recover remaining trees
#      using a previously saved model
#   3. Leave-One-Tree-Out (LOTO) cross-validation at tree level
#   4. Estimate leaf area (MUNC method) with geometric error propagation
#      and FP/FN classification uncertainty range
#   5. Compute crown-projected LAI via ITSMe
#
# Reference for ITSMe: Terryn et al. (2023), Methods Ecol Evol 14:231-241.
#   doi:10.1111/2041-210X.14026

# PACKAGES

packages <- c("lidR", "ranger", "caret", "FNN", "ggplot2",
              "dplyr", "tidyr", "viridis", "tools", "progress", "ITSMe")

# ITSMe is on GitHub, not CRAN
if (!"ITSMe" %in% installed.packages()[, "Package"]) {
  if (!"devtools" %in% installed.packages()[, "Package"])
    install.packages("devtools")
  devtools::install_github("lmterryn/ITSMe")
}

missing_pkg <- packages[packages != "ITSMe" &
                          !packages %in% installed.packages()[, "Package"]]
if (length(missing_pkg) > 0)
  install.packages(missing_pkg, dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

# 1. USER PARAMETERS  <-- edit these paths before running

params <- list(

  # --- Training files (annotated wood and leaf point clouds) ---
  training_wood_path  = "path/to/wood_training.las",
  training_leaf_path  = "path/to/leaf_training.las",

  # --- Per-tree annotated files for LOTO validation (wood_XXXX.las / leaf_XXXX.las) ---
  training_dir = "path/to/training_set",

  # --- Trees included in LOTO (must have matching wood_*.las and leaf_*.las) ---
  loto_tree_ids = c("id1", "id2", "id3", "id4"),

  # --- Input folder: one .las per tree to classify ---
  input_dir  = "path/to/classified",

  # --- Output folder (classified LAS, CSVs, plots) ---
  output_dir = "path/to/output_RF",

  # --- Geometric feature scales (m): from fine leaf detail to branch diameter ---
  scale_geom = c(0.05, 0.10, 0.20, 0.40),

  # --- Radiometric features to include ---
  use_rgb        = TRUE,   # R, G, B (normalized)
  use_intensity  = TRUE,   # laser return intensity
  use_return_num = TRUE,   # Return Number, Number of Returns, ratio
  use_height     = TRUE,   # normalized Z
  use_density    = TRUE,   # local point density

  # --- Random Forest ---
  n_trees       = 500,     # trees in forest
  n_fold_cv     = 5,       # folds for point-level cross-validation
  max_train_pts = 50000,   # max points per class (down-samples if needed)
  max_test_pts  = 50000,   # max points per class in LOTO test set
  seed          = 42,

  # --- Classification threshold: P(leaf) >= threshold -> leaf ---
  leaf_threshold = 0.50,

  # --- Class codes ---
  class_leaf = 1L,
  class_wood = 2L,

  # --- MUNC parameters (calibrated for Leica BLK360) ---
  k_normals  = 20,
  voxel_size = 0.007,      # m

  # --- ITSMe concavity parameter for projected_area_pc ---
  concavity  = 2
)

# Co-registration error per tree (Cyclone REGISTER 360 PLUS group error, m)
coreg <- data.frame(
  tree_name   = c("names of the tree files"),
  n_setup     = c(2L, 2L, 3L, 3L, 3L, 3L, 4L, 2L, 2L, 2L),
  err_coreg_m = c(0.023, 0.008, 0.008, 0.008,
                  0.009, 0.009, 0.007, 0.008,
                  0.008, 0.008),
  stringsAsFactors = FALSE
)

if (!dir.exists(params$output_dir))
  dir.create(params$output_dir, recursive = TRUE)

# 2. HELPER FUNCTIONS

# -2.1 Co-registration error lookup 
get_coreg_err <- function(tree_name, coreg_table) {
  idx <- match(tree_name, coreg_table$tree_name)
  if (is.na(idx)) {
    warning("Co-registration error not found for: ", tree_name,
            " -> using default 0.009 m")
    return(0.009)
  }
  coreg_table$err_coreg_m[idx]
}

# -2.2 Geometric feature extraction (multi-scale PCA)
extract_features <- function(las, scale_geom, params, label = "") {
  cat(sprintf("  Extracting features%s (%d points)...\n",
              label, npoints(las)))

  df  <- las@data
  xyz <- cbind(df$X, df$Y, df$Z)
  n   <- nrow(df)

  feat <- data.frame(row.names = seq_len(n))

  if (params$use_height)
    feat$Z_norm <- df$Z

  if (params$use_intensity && "Intensity" %in% names(df))
    feat$intensity <- as.numeric(df$Intensity)

  rn_field  <- intersect(c("ReturnNumber",    "Return.Number",   "Return_Number"),  names(df))
  nor_field <- intersect(c("NumberOfReturns", "Number.of.Returns","Number_of_Returns","NumberReturns"), names(df))

  if (params$use_return_num && length(rn_field)  > 0) feat$return_num <- as.numeric(df[[rn_field[1]]])
  if (params$use_return_num && length(nor_field) > 0) feat$n_returns  <- as.numeric(df[[nor_field[1]]])
  if (!is.null(feat$return_num) && !is.null(feat$n_returns))
    feat$return_ratio <- feat$return_num / (feat$n_returns + 1e-10)

  if (params$use_rgb) {
    if (all(c("R","G","B") %in% names(df))) {
      mx <- if (max(df$R, na.rm = TRUE) > 255) 65535 else 255
      feat$R   <- as.numeric(df$R) / mx
      feat$G   <- as.numeric(df$G) / mx
      feat$B   <- as.numeric(df$B) / mx
      feat$ExG <- 2 * feat$G - feat$R - feat$B   # Excess Green index
    } else {
      cat("  WARN: RGB not found in LAS\n")
    }
  }

  # Multi-scale geometric features (Demantke et al. 2011)
  cat("  Multi-scale geometric features")
  for (r in scale_geom) {
    cat(sprintf(" [%.2fm]", r))
    k   <- min(30, n - 1)
    knn <- get.knn(xyz, k = k)

    lin_v  <- numeric(n)
    plan_v <- numeric(n)
    sfer_v <- numeric(n)
    curv_v <- numeric(n)
    dens_v <- numeric(n)

    pb <- progress_bar$new(total = n, format = "  [:bar] :percent",
                           clear = TRUE, width = 40)
    for (i in seq_len(n)) {
      pb$tick()
      dists <- knn$nn.dist[i, ]
      idx   <- knn$nn.index[i, dists <= r]

      if (length(idx) < 3) { lin_v[i] <- 1; next }

      pts    <- xyz[idx, , drop = FALSE]
      centro <- colMeans(pts)
      D      <- sweep(pts, 2, centro)
      eig    <- tryCatch(eigen(crossprod(D) / nrow(D), symmetric = TRUE),
                         error = function(e) NULL)
      if (is.null(eig)) { lin_v[i] <- 1; next }

      lam <- sort(eig$values, decreasing = TRUE)
      lam <- pmax(lam, 0)
      tot <- sum(lam) + 1e-10

      lin_v[i]  <- (lam[1] - lam[2]) / tot
      plan_v[i] <- (lam[2] - lam[3]) / tot
      sfer_v[i] <- lam[3] / tot
      curv_v[i] <- lam[3] / (lam[1] + 1e-10)
      dens_v[i] <- length(idx) / (4/3 * pi * r^3)
    }

    sufx <- gsub("\\.", "", sprintf("_r%03d", round(r * 100)))
    feat[[paste0("lin",  sufx)]] <- lin_v
    feat[[paste0("plan", sufx)]] <- plan_v
    feat[[paste0("sfer", sufx)]] <- sfer_v
    feat[[paste0("curv", sufx)]] <- curv_v
    if (params$use_density)
      feat[[paste0("dens", sufx)]] <- dens_v
  }
  cat("\n")
  feat
}

# -2.3 MUNC helper functions 

compute_normals <- function(xyz, k = 20) {
  n   <- nrow(xyz)
  N   <- matrix(0, n, 3)
  knn <- get.knn(xyz, k = k)
  pb  <- progress_bar$new(total = n, format = "  [:bar] :percent ETA: :eta",
                           clear = FALSE)
  for (i in seq_len(n)) {
    pb$tick()
    vicini <- xyz[knn$nn.index[i, ], , drop = FALSE]
    centro <- colMeans(vicini)
    D      <- sweep(vicini, 2, centro)
    eig    <- tryCatch(eigen(crossprod(D) / k, symmetric = TRUE),
                       error = function(e) NULL)
    if (!is.null(eig)) N[i, ] <- eig$vectors[, 3]
  }
  N
}

voxelise <- function(xyz, normals, weights, voxel) {
  # Reduce point cloud to one representative per voxel (Fan et al. 2025 eq. 7)
  ix  <- floor(xyz[, 1] / voxel)
  iy  <- floor(xyz[, 2] / voxel)
  iz  <- floor(xyz[, 3] / voxel)
  key <- paste(ix, iy, iz, sep = "_")

  voxel_keys <- unique(key)
  nv <- length(voxel_keys)
  N_out <- matrix(0, nv, 3)
  p_out <- if (!is.null(weights)) numeric(nv) else NULL

  for (k in seq_along(voxel_keys)) {
    ii    <- which(key == voxel_keys[k])
    N_avg <- colMeans(normals[ii, , drop = FALSE])
    len   <- sqrt(sum(N_avg^2)) + 1e-10
    N_out[k, ] <- N_avg / len
    if (!is.null(weights)) p_out[k] <- mean(weights[ii], na.rm = TRUE)
  }
  list(normals = N_out, weights = p_out, n_voxel = nv)
}

munc_area <- function(xyz, normals, voxel, weights = NULL) {
  vox   <- voxelise(xyz, normals, weights, voxel)
  N     <- vox$normals
  p     <- vox$weights
  len   <- sqrt(rowSums(N^2)) + 1e-10
  N_hat <- N / len
  munc  <- apply(abs(N_hat), 1, max)
  munc[munc < 1e-6] <- 1e-6
  da    <- (voxel^2) / munc
  if (!is.null(p)) da_w <- da * pmin(pmax(p, 0), 1) else da_w <- da
  n   <- length(da_w)
  tot <- sum(da_w)
  se  <- sd(da_w) * sqrt(n)
  list(area_m2  = tot,
       se_m2    = se,
       cv_pct   = round(100 * se / (tot + 1e-10), 2),
       n_voxel  = vox$n_voxel,
       n_points = n)
}

fp_fn_range <- function(xyz_leaf, norm_leaf, prob_leaf,
                         xyz_wood, norm_wood, prob_wood, voxel) {
  # FP: leaf points likely misclassified as wood (weight = P(wood) = 1-P(leaf))
  w_fp   <- 1 - pmin(pmax(prob_leaf, 0), 1)
  vox_fp <- voxelise(xyz_leaf, norm_leaf, w_fp, voxel)
  mfp    <- apply(abs(vox_fp$normals /
                      (sqrt(rowSums(vox_fp$normals^2)) + 1e-10)), 1, max)
  mfp[mfp < 1e-6] <- 1e-6
  fp_m2  <- sum((voxel^2) / mfp * pmin(pmax(vox_fp$weights, 0), 1))

  # FN: wood points likely misclassified as leaf (weight = P(leaf) stored in prob_wood)
  fn_m2 <- 0
  if (!is.null(xyz_wood) && nrow(xyz_wood) >= 10) {
    w_fn   <- pmin(pmax(prob_wood, 0), 1)
    vox_fn <- voxelise(xyz_wood, norm_wood, w_fn, voxel)
    mfn    <- apply(abs(vox_fn$normals /
                        (sqrt(rowSums(vox_fn$normals^2)) + 1e-10)), 1, max)
    mfn[mfn < 1e-6] <- 1e-6
    fn_m2  <- sum((voxel^2) / mfn * pmin(pmax(vox_fn$weights, 0), 1))
  }
  list(fp_m2 = fp_m2, fn_m2 = fn_m2)
}

# -2.4 Geometric error propagation 
# Three independent components (co-registration, voxelization, normal estimation)
# combined in quadrature. See rf_report.R inline comments for derivation.
propagate_errors <- function(area_m2, n_voxel, voxel_size, k_normals,
                              err_coreg_m, sigma_scanner = 0.006) {
  r_local <- voxel_size * sqrt(k_normals)   # local geometric scale

  sigma_coreg   <- (area_m2 / sqrt(n_voxel)) * (err_coreg_m / r_local)
  sigma_voxel   <- (area_m2 / sqrt(n_voxel)) * ((voxel_size / (2 * sqrt(3))) / r_local)
  delta_norm    <- min(sigma_scanner / r_local, pi / 6)
  sigma_normals <- (area_m2 / sqrt(n_voxel)) * sin(delta_norm)

  sigma_geo <- sqrt(sigma_coreg^2 + sigma_voxel^2 + sigma_normals^2)

  list(sigma_coreg_m2   = round(sigma_coreg,   4),
       sigma_voxel_m2   = round(sigma_voxel,   4),
       sigma_normals_m2 = round(sigma_normals, 4),
       sigma_geo_m2     = round(sigma_geo,     4),
       ic95_geo_m2      = round(1.96 * sigma_geo, 4))
}

# 3. RF TRAINING AND BATCH CLASSIFICATION

load_training_set <- function(params) {
  cat("=== LOADING TRAINING SET ===\n")
  las_wood <- readLAS(params$training_wood_path, select = "*")
  las_leaf <- readLAS(params$training_leaf_path, select = "*")
  cat(sprintf("  Wood:  %d points\n", npoints(las_wood)))
  cat(sprintf("  Leaf:  %d points\n", npoints(las_leaf)))

  set.seed(params$seed)
  n_max <- params$max_train_pts
  if (npoints(las_wood) > n_max) las_wood <- las_wood[sample(npoints(las_wood), n_max)]
  if (npoints(las_leaf) > n_max) las_leaf <- las_leaf[sample(npoints(las_leaf), n_max)]

  feat_w <- extract_features(las_wood, params$scale_geom, params, " (wood)")
  feat_l <- extract_features(las_leaf, params$scale_geom, params, " (leaf)")

  feat_w$class <- "wood";  feat_l$class <- "leaf"
  df <- bind_rows(feat_w, feat_l)
  df$class <- factor(df$class, levels = c("leaf", "wood"))

  # Remove zero-variance or all-NA columns
  df <- df[, colSums(is.na(df)) == 0, drop = FALSE]
  df <- df[, apply(df[, -ncol(df), drop = FALSE], 2,
                   function(x) is.character(x) || var(as.numeric(x), na.rm = TRUE) > 1e-10),
           drop = FALSE]

  cat(sprintf("  Training set: %d points x %d features\n", nrow(df), ncol(df) - 1))
  df
}

train_RF <- function(df, params) {
  cat("=== TRAINING RANDOM FOREST ===\n")
  set.seed(params$seed)
  X <- df[, names(df) != "class", drop = FALSE]
  y <- df$class

  # Point-level k-fold CV (internal consistency check)
  cat(sprintf("  %d-fold cross-validation...\n", params$n_fold_cv))
  folds  <- createFolds(y, k = params$n_fold_cv, list = TRUE)
  acc_cv <- numeric(params$n_fold_cv)
  ba_cv  <- numeric(params$n_fold_cv)

  for (f in seq_along(folds)) {
    idx_val <- folds[[f]]
    idx_trn <- setdiff(seq_len(nrow(df)), idx_val)
    mod_f   <- ranger(class ~ ., data = df[idx_trn, ],
                      num.trees = 200, probability = TRUE, seed = params$seed)
    prob_v  <- predict(mod_f, data = df[idx_val, ])$predictions[, "leaf"]
    pred_v  <- factor(ifelse(prob_v >= params$leaf_threshold, "leaf", "wood"),
                      levels = c("leaf", "wood"))
    cm      <- confusionMatrix(pred_v, y[idx_val], positive = "leaf")
    acc_cv[f] <- cm$overall["Accuracy"]
    ba_cv[f]  <- cm$byClass["Balanced Accuracy"]
    cat(sprintf("    Fold %d: Accuracy=%.3f | Balanced Accuracy=%.3f\n",
                f, acc_cv[f], ba_cv[f]))
  }
  cat(sprintf("  CV mean: Accuracy=%.3f±%.3f | BA=%.3f±%.3f\n\n",
              mean(acc_cv), sd(acc_cv), mean(ba_cv), sd(ba_cv)))

  # Final model on full training set
  cat(sprintf("  Training final model (%d trees)...\n", params$n_trees))
  model <- ranger(class ~ ., data = df,
                  num.trees   = params$n_trees,
                  probability = TRUE,
                  importance  = "permutation",
                  seed        = params$seed)
  cat("  Done.\n\n")

  list(model   = model,
       acc_cv  = acc_cv,
       ba_cv   = ba_cv,
       feature = names(X),
       n_train = nrow(df))
}

classify_las <- function(path, rf_obj, params, outdir) {
  tree_name <- file_path_sans_ext(basename(path))
  cat(sprintf("\n--- %s ---\n", tree_name))
  t0 <- proc.time()

  las <- readLAS(path, select = "*")
  if (is.null(las) || npoints(las) == 0) { cat("  ERROR: empty file\n"); return(NULL) }

  feat <- extract_features(las, params$scale_geom, params)
  feat_common <- intersect(rf_obj$feature, names(feat))
  feat_missing <- setdiff(rf_obj$feature, names(feat))
  if (length(feat_missing) > 0)
    cat("  WARN: missing features:", paste(feat_missing, collapse = ", "), "\n")

  feat_pred <- feat[, feat_common, drop = FALSE]
  for (col in names(feat_pred))
    feat_pred[[col]][is.na(feat_pred[[col]])] <- median(feat_pred[[col]], na.rm = TRUE)

  cat("  RF classification...\n")
  pred      <- predict(rf_obj$model, data = feat_pred)
  prob_leaf <- pred$predictions[, "leaf"]
  class_rf  <- ifelse(prob_leaf >= params$leaf_threshold, params$class_leaf, params$class_wood)

  n_leaf <- sum(class_rf == params$class_leaf)
  n_wood <- sum(class_rf == params$class_wood)
  cat(sprintf("  Leaf: %d (%.1f%%) | Wood: %d (%.1f%%)\n",
              n_leaf, 100 * n_leaf / npoints(las),
              n_wood, 100 * n_wood / npoints(las)))

  las <- add_lasattribute(las, as.integer(class_rf),      "RF.class", "RF classification: 1=leaf 2=wood")
  las <- add_lasattribute(las, as.numeric(round(prob_leaf, 4)), "RF.prob",  "RF P(leaf)")

  # Optional: compare with CANUPO classification if present
  concordance <- NA_real_
  if ("CANUPO.class" %in% names(las@data)) {
    canupo      <- as.integer(las@data[["CANUPO.class"]])
    concordance <- round(100 * mean(class_rf == canupo, na.rm = TRUE), 1)
    cat(sprintf("  RF vs CANUPO agreement: %.1f%%\n", concordance))
  }

  out_path <- file.path(outdir, paste0(tree_name, "_RF.las"))
  writeLAS(las, out_path)
  cat(sprintf("  Saved: %s\n", basename(out_path)))
  cat(sprintf("  Done in %.1f s\n", (proc.time() - t0)["elapsed"]))

  list(tree_name       = tree_name,
       n_pts           = npoints(las),
       n_leaf          = n_leaf,
       n_wood          = n_wood,
       pct_leaf        = round(100 * n_leaf / npoints(las), 1),
       prob_mean       = round(mean(prob_leaf[class_rf == params$class_leaf]), 3),
       concordance_pct = concordance,
       elapsed_sec     = round((proc.time() - t0)["elapsed"], 1),
       status          = "OK")
}

# Tiled classification for large point clouds (reduces RAM usage)
classify_las_tiled <- function(path, rf_obj, params, outdir,
                                tile_size = 10, overlap = 0.5) {
  tree_name <- file_path_sans_ext(basename(path))
  cat(sprintf("\n--- %s (TILED mode) ---\n", tree_name))
  t0 <- proc.time()

  las <- readLAS(path, select = "*")
  if (is.null(las) || npoints(las) == 0) { cat("  ERROR: empty file\n"); return(NULL) }
  cat(sprintf("  Total points: %d\n", npoints(las)))

  xmin <- min(las@data$X); xmax <- max(las@data$X)
  ymin <- min(las@data$Y); ymax <- max(las@data$Y)
  xs   <- seq(xmin, xmax, by = tile_size); if (tail(xs, 1) < xmax) xs <- c(xs, xmax)
  ys   <- seq(ymin, ymax, by = tile_size); if (tail(ys, 1) < ymax) ys <- c(ys, ymax)

  n_pts     <- npoints(las)
  class_rf  <- integer(n_pts)
  prob_leaf <- numeric(n_pts)
  count     <- integer(n_pts)

  tile_idx <- 0
  for (i in seq_len(length(xs) - 1)) {
    for (j in seq_len(length(ys) - 1)) {
      tile_idx <- tile_idx + 1
      x1 <- xs[i] - overlap;  x2 <- xs[i+1] + overlap
      y1 <- ys[j] - overlap;  y2 <- ys[j+1] + overlap
      idx_tile <- which(las@data$X >= x1 & las@data$X <= x2 &
                        las@data$Y >= y1 & las@data$Y <= y2)
      if (length(idx_tile) < 10) next
      cat(sprintf("  Tile %d: %d points\n", tile_idx, length(idx_tile)))

      las_tile   <- las[idx_tile]
      feat_tile  <- extract_features(las_tile, params$scale_geom, params)
      feat_c     <- intersect(rf_obj$feature, names(feat_tile))
      feat_pred  <- feat_tile[, feat_c, drop = FALSE]
      for (col in names(feat_pred))
        feat_pred[[col]][is.na(feat_pred[[col]])] <- median(feat_pred[[col]], na.rm = TRUE)

      pred_t  <- predict(rf_obj$model, data = feat_pred)
      prob_i  <- pred_t$predictions[, "leaf"]
      class_i <- ifelse(prob_i >= params$leaf_threshold, params$class_leaf, params$class_wood)

      prob_leaf[idx_tile] <- prob_leaf[idx_tile] + prob_i
      class_rf[idx_tile]  <- class_rf[idx_tile]  + as.integer(class_i)
      count[idx_tile]     <- count[idx_tile] + 1L
      rm(las_tile, feat_tile, feat_pred, pred_t, prob_i, class_i); gc()
    }
  }

  # Average probabilities for points covered by multiple tiles (overlap zones)
  proc <- count > 0
  prob_leaf[proc] <- prob_leaf[proc] / count[proc]
  class_rf[proc]  <- ifelse(prob_leaf[proc] >= params$leaf_threshold,
                             params$class_leaf, params$class_wood)
  if (any(!proc)) {
    cat(sprintf("  WARN: %d uncovered points -> assigned to wood\n", sum(!proc)))
    class_rf[!proc] <- params$class_wood
  }

  n_leaf <- sum(class_rf == params$class_leaf)
  n_wood <- sum(class_rf == params$class_wood)
  cat(sprintf("  Leaf: %d (%.1f%%) | Wood: %d (%.1f%%)\n",
              n_leaf, 100 * n_leaf / n_pts, n_wood, 100 * n_wood / n_pts))

  las <- add_lasattribute(las, as.integer(class_rf),           "RF.class", "RF classification: 1=leaf 2=wood")
  las <- add_lasattribute(las, as.numeric(round(prob_leaf, 4)), "RF.prob",  "RF P(leaf)")

  concordance <- NA_real_
  if ("CANUPO.class" %in% names(las@data)) {
    canupo      <- as.integer(las@data[["CANUPO.class"]])
    concordance <- round(100 * mean(class_rf == canupo, na.rm = TRUE), 1)
    cat(sprintf("  RF vs CANUPO agreement: %.1f%%\n", concordance))
  }

  out_path <- file.path(outdir, paste0(tree_name, "_RF.las"))
  writeLAS(las, out_path)
  cat(sprintf("  Saved: %s  (%.1f s)\n", basename(out_path),
              (proc.time() - t0)["elapsed"]))

  list(tree_name       = tree_name,
       n_pts           = n_pts,
       n_leaf          = n_leaf,
       n_wood          = n_wood,
       pct_leaf        = round(100 * n_leaf / n_pts, 1),
       prob_mean       = NA_real_,
       concordance_pct = concordance,
       elapsed_sec     = round((proc.time() - t0)["elapsed"], 1),
       status          = "OK (tiled)")
}

run_classification_batch <- function(params) {
  # Check input files
  if (!file.exists(params$training_wood_path))
    stop("Training wood file not found: ", params$training_wood_path)
  if (!file.exists(params$training_leaf_path))
    stop("Training leaf file not found: ", params$training_leaf_path)
  if (!dir.exists(params$input_dir))
    stop("Input directory not found: ", params$input_dir)

  # Train
  df_train <- load_training_set(params)
  rf_obj   <- train_RF(df_train, params)

  # Save model for future use (e.g., classify additional trees without retraining)
  model_path <- file.path(params$output_dir, "modello_RF.rds")
  saveRDS(rf_obj, model_path)
  cat(sprintf("Model saved: %s\n\n", model_path))

  # Plot feature importance
  imp    <- sort(rf_obj$model$variable.importance, decreasing = TRUE)
  df_imp <- data.frame(feature = names(imp), importance = imp,
                       stringsAsFactors = FALSE)
  df_imp$feature <- factor(df_imp$feature, levels = rev(df_imp$feature))
  df_imp$type <- case_when(
    grepl("^lin|^plan|^sfer|^curv|^dens", df_imp$feature) ~ "Geometric",
    grepl("^R$|^G$|^B$|^ExG",            df_imp$feature) ~ "RGB",
    grepl("intensity",                    df_imp$feature) ~ "Intensity",
    grepl("return",                       df_imp$feature) ~ "Returns",
    grepl("Z_norm",                       df_imp$feature) ~ "Height",
    TRUE                                                   ~ "Other"
  )
  p_imp <- ggplot(df_imp, aes(x = feature, y = importance, fill = type)) +
    geom_col() + scale_fill_viridis(discrete = TRUE) + coord_flip() +
    labs(title = "Feature importance – Random Forest",
         subtitle = "Permutation importance",
         x = "", y = "Importance", fill = "Type") +
    theme_minimal(base_size = 11)
  ggsave(file.path(params$output_dir, "00_RF_feature_importance.png"),
         p_imp, width = 10, height = max(5, nrow(df_imp) * 0.25 + 2), dpi = 150)

  # Classify all trees
  files <- list.files(params$input_dir, pattern = "\\.las$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No .las files found in: ", params$input_dir)
  cat(sprintf("=== BATCH CLASSIFICATION (%d files) ===\n", length(files)))

  results <- list()
  t_batch <- proc.time()
  for (i in seq_along(files)) {
    cat(sprintf("\n[%d/%d] ", i, length(files)))
    # Use tiled mode for large clouds (> 1.5 M points)
    n_est <- tryCatch({
      h <- readLASheader(files[i])
      h@PHB[["Number of point records"]]
    }, error = function(e) 0)

    res <- tryCatch(
      if (n_est > 1500000)
        classify_las_tiled(files[i], rf_obj, params, params$output_dir)
      else
        classify_las(files[i], rf_obj, params, params$output_dir),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        list(tree_name = file_path_sans_ext(basename(files[i])),
             status    = paste("ERR:", conditionMessage(e)))
      }
    )
    if (!is.null(res)) results[[i]] <- res
    cat(rep("-", 60), "\n", sep = "")
  }

  summary_df <- bind_rows(lapply(results, as.data.frame))
  write.csv2(summary_df, file.path(params$output_dir, "00_RF_classify_summary.csv"),
             row.names = FALSE)

  elapsed_tot <- (proc.time() - t_batch)["elapsed"]
  cat(sprintf("\n=== CLASSIFICATION DONE in %.1f min ===\n", elapsed_tot / 60))
  cat(sprintf("Trees OK: %d / %d\n",
              sum(summary_df$status %in% c("OK","OK (tiled)"), na.rm = TRUE),
              nrow(summary_df)))

  invisible(list(rf = rf_obj, summary = summary_df))
}

# Convenience function: classify additional trees using a previously saved model
# (run when some trees were added after the initial batch)
recover_classification <- function(params, model_path = NULL) {
  if (is.null(model_path))
    model_path <- file.path(params$output_dir, "modello_RF.rds")
  if (!file.exists(model_path))
    stop("Model not found: ", model_path,
         "\nRun run_classification_batch() first.")

  cat("=== LOADING SAVED RF MODEL ===\n")
  rf_obj <- readRDS(model_path)
  cat(sprintf("  Model loaded: %d trees | %d features\n",
              rf_obj$model$num.trees, length(rf_obj$feature)))

  files <- list.files(params$input_dir, pattern = "\\.las$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No .las files found in: ", params$input_dir)
  cat(sprintf("=== RECOVERY BATCH (%d trees) ===\n", length(files)))

  summary_df <- data.frame()
  t_batch    <- proc.time()
  for (i in seq_along(files)) {
    cat(sprintf("\n[%d/%d] ", i, length(files)))
    n_est <- tryCatch({
      h <- readLASheader(files[i])
      h@PHB[["Number of point records"]]
    }, error = function(e) 0)
    riga <- tryCatch(
      if (n_est > 1500000)
        classify_las_tiled(files[i], rf_obj, params, params$output_dir)
      else
        classify_las(files[i], rf_obj, params, params$output_dir),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        data.frame(tree_name = file_path_sans_ext(basename(files[i])),
                   status    = paste("ERR:", conditionMessage(e)),
                   stringsAsFactors = FALSE)
      }
    )
    if (!is.null(riga)) summary_df <- bind_rows(summary_df, as.data.frame(riga))
    cat(rep("-", 60), "\n", sep = "")
  }
  write.csv2(summary_df, file.path(params$output_dir, "00_RF_recovery_summary.csv"),
             row.names = FALSE)
  cat(sprintf("\n=== RECOVERY DONE in %.1f min ===\n",
              (proc.time() - t_batch)["elapsed"] / 60))
  invisible(summary_df)
}

# 4. LEAVE-ONE-TREE-OUT (LOTO) CROSS-VALIDATION

load_tree_data <- function(tree_id, training_dir, params, max_pts = NULL) {
  cat(sprintf("\n  Loading %s...\n", tree_id))
  path_w <- file.path(training_dir, sprintf("wood_%s.las", tree_id))
  path_l <- file.path(training_dir, sprintf("leaf_%s.las", tree_id))
  if (!file.exists(path_w)) stop("Missing file: ", path_w)
  if (!file.exists(path_l)) stop("Missing file: ", path_l)

  las_w <- readLAS(path_w, select = "*")
  las_l <- readLAS(path_l, select = "*")
  cat(sprintf("    Wood:  %d points | Leaf: %d points\n",
              npoints(las_w), npoints(las_l)))

  n_max <- if (!is.null(max_pts)) max_pts else params$max_train_pts
  set.seed(params$seed)
  if (npoints(las_w) > n_max) las_w <- las_w[sample(npoints(las_w), n_max)]
  if (npoints(las_l) > n_max) las_l <- las_l[sample(npoints(las_l), n_max)]

  feat_w <- extract_features(las_w, params$scale_geom, params,
                              sprintf(" wood (%s)", tree_id))
  feat_l <- extract_features(las_l, params$scale_geom, params,
                              sprintf(" leaf (%s)", tree_id))

  feat_w$class <- "wood";  feat_l$class <- "leaf"
  df <- bind_rows(feat_w, feat_l)
  df$class   <- factor(df$class, levels = c("leaf", "wood"))
  df$tree_id <- tree_id
  df
}

run_loto <- function(params) {
  tree_ids <- params$loto_tree_ids
  n        <- length(tree_ids)

  # Verify all annotated files exist
  for (id in tree_ids) {
    pw <- file.path(params$training_dir, sprintf("wood_%s.las", id))
    pl <- file.path(params$training_dir, sprintf("leaf_%s.las", id))
    if (!file.exists(pw)) stop("Missing file: ", pw)
    if (!file.exists(pl)) stop("Missing file: ", pl)
  }
  cat(sprintf("=== LEAVE-ONE-TREE-OUT VALIDATION (n=%d trees) ===\n\n", n))

  # Load all trees
  cat("--- Loading all trees ---\n")
  data_list <- setNames(vector("list", n), tree_ids)
  for (id in tree_ids) {
    data_list[[id]] <- tryCatch(
      load_tree_data(id, params$training_dir, params),
      error = function(e) { cat(sprintf("ERROR loading %s: %s\n", id, e$message)); NULL }
    )
  }
  failed <- names(which(sapply(data_list, is.null)))
  if (length(failed) > 0) stop("Could not load: ", paste(failed, collapse = ", "))

  # Keep only features present in all trees
  feat_cols <- lapply(data_list, function(d)
    names(d)[!names(d) %in% c("class", "tree_id")])
  common_feats <- Reduce(intersect, feat_cols)
  cat(sprintf("\nCommon features: %d\n\n", length(common_feats)))

  results <- list()

  for (fold in seq_len(n)) {
    test_id   <- tree_ids[fold]
    train_ids <- tree_ids[-fold]

    cat(sprintf("\n========================================\n"))
    cat(sprintf("FOLD %d/%d  |  Training: %s  |  Test: %s\n",
                fold, n, paste(train_ids, collapse = "+"), test_id))
    cat(sprintf("========================================\n"))

    df_train <- bind_rows(lapply(train_ids, function(id)
      data_list[[id]][, c(common_feats, "class")]
    ))
    df_train$class <- factor(df_train$class, levels = c("leaf", "wood"))

    df_test  <- data_list[[test_id]][, c(common_feats, "class")]
    df_test$class <- factor(df_test$class, levels = c("leaf", "wood"))

    cat(sprintf("  Training: %d points (leaf=%d, wood=%d)\n",
                nrow(df_train), sum(df_train$class=="leaf"), sum(df_train$class=="wood")))
    cat(sprintf("  Test:     %d points (leaf=%d, wood=%d)\n",
                nrow(df_test),  sum(df_test$class=="leaf"),  sum(df_test$class=="wood")))

    cat(sprintf("  Training RF (%d trees)...\n", params$n_trees))
    set.seed(params$seed)
    mod <- ranger(class ~ ., data = df_train,
                  num.trees   = params$n_trees,
                  probability = TRUE,
                  seed        = params$seed)

    feat_test <- df_test[, common_feats, drop = FALSE]
    for (col in names(feat_test))
      feat_test[[col]][is.na(feat_test[[col]])] <- median(feat_test[[col]], na.rm = TRUE)

    pred        <- predict(mod, data = feat_test)
    prob_leaf   <- pred$predictions[, "leaf"]
    class_pred  <- factor(ifelse(prob_leaf >= params$leaf_threshold, "leaf", "wood"),
                          levels = c("leaf", "wood"))
    cm <- confusionMatrix(class_pred, df_test$class, positive = "leaf")

    BA   <- cm$byClass["Balanced Accuracy"]
    sens <- cm$byClass["Sensitivity"]
    spec <- cm$byClass["Specificity"]
    prec_l <- cm$byClass["Precision"]
    prec_w <- cm$byClass["Neg Pred Value"]
    acc    <- cm$overall["Accuracy"]

    cat(sprintf("\n  Results fold %d (%s):\n", fold, test_id))
    cat(sprintf("    Balanced Accuracy:  %.4f\n", BA))
    cat(sprintf("    Sensitivity (leaf): %.4f\n", sens))
    cat(sprintf("    Specificity (wood): %.4f\n", spec))
    cat(sprintf("    Precision (leaf):   %.4f\n", prec_l))
    cat(sprintf("    Precision (wood):   %.4f\n", prec_w))
    cat(sprintf("    Accuracy:           %.4f\n", acc))
    cat("\n  Confusion matrix:\n"); print(cm$table)

    results[[fold]] <- data.frame(
      fold              = fold,
      test_tree         = test_id,
      train_trees       = paste(train_ids, collapse = "+"),
      n_train           = nrow(df_train),
      n_test            = nrow(df_test),
      n_test_leaf       = sum(df_test$class == "leaf"),
      n_test_wood       = sum(df_test$class == "wood"),
      accuracy          = round(acc,    4),
      balanced_accuracy = round(BA,     4),
      sensitivity_leaf  = round(sens,   4),
      specificity_wood  = round(spec,   4),
      precision_leaf    = round(prec_l, 4),
      precision_wood    = round(prec_w, 4),
      TP = cm$table["leaf","leaf"],
      FP = cm$table["leaf","wood"],
      FN = cm$table["wood","leaf"],
      TN = cm$table["wood","wood"],
      stringsAsFactors = FALSE
    )
  }

  df_res <- bind_rows(results)

  # Summary
  metrics  <- c("balanced_accuracy","accuracy","sensitivity_leaf",
                "specificity_wood","precision_leaf","precision_wood")
  labels   <- c("Balanced Accuracy","Accuracy","Sensitivity (leaf)",
                "Specificity (wood)","Precision (leaf)","Precision (wood)")
  cat("\n\n==========================================\n")
  cat("LOTO SUMMARY\n")
  cat("==========================================\n")
  for (k in seq_along(metrics)) {
    v <- df_res[[metrics[k]]]
    cat(sprintf("%-22s %.4f ± %.4f  [%.4f - %.4f]\n",
                paste0(labels[k], ":"),
                mean(v, na.rm=TRUE), sd(v, na.rm=TRUE),
                min(v, na.rm=TRUE),  max(v, na.rm=TRUE)))
  }

  # Save outputs
  loto_output_dir <- file.path(params$output_dir, "LOTO_validation")
  if (!dir.exists(loto_output_dir)) dir.create(loto_output_dir, recursive=TRUE)

  write.csv2(df_res, file.path(loto_output_dir, "00_LOTO_results.csv"), row.names = FALSE)

  sink(file.path(loto_output_dir, "00_LOTO_summary.txt"))
  cat("=== LEAVE-ONE-TREE-OUT VALIDATION ===\n")
  cat(sprintf("Date: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("Trees: %s\n", paste(params$loto_tree_ids, collapse=", ")))
  cat(sprintf("n_trees RF: %d | threshold: %.2f\n\n", params$n_trees, params$leaf_threshold))
  cat("--- Per-fold metrics ---\n")
  print(df_res[, c("fold","test_tree","balanced_accuracy","sensitivity_leaf",
                   "specificity_wood","precision_leaf","precision_wood",
                   "TP","FP","FN","TN")], row.names=FALSE)
  cat("\n--- Summary (mean ± SD  [min - max]) ---\n")
  for (k in seq_along(metrics)) {
    v <- df_res[[metrics[k]]]
    cat(sprintf("%-22s %.4f ± %.4f  [%.4f - %.4f]\n",
                paste0(labels[k], ":"), mean(v, na.rm=TRUE), sd(v, na.rm=TRUE),
                min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
  }
  sink()

  # LOTO metrics plot
  df_plot <- df_res %>%
    select(fold, test_tree, balanced_accuracy, sensitivity_leaf,
           specificity_wood, precision_leaf, precision_wood) %>%
    pivot_longer(cols = -c(fold, test_tree), names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
      balanced_accuracy = "Balanced Accuracy",
      sensitivity_leaf  = "Sensitivity (leaf)",
      specificity_wood  = "Specificity (wood)",
      precision_leaf    = "Precision (leaf)",
      precision_wood    = "Precision (wood)"
    ))

  df_mean <- df_plot %>%
    group_by(metric) %>%
    summarise(mean_val = mean(value, na.rm=TRUE), .groups="drop")

  p_loto <- ggplot(df_plot, aes(x = test_tree, y = value, fill = test_tree)) +
    geom_col(width = 0.6) +
    geom_hline(data = df_mean, aes(yintercept = mean_val),
               linetype = "dashed", color = "gray40", linewidth = 0.7) +
    geom_text(data = df_mean,
              aes(x = Inf, y = mean_val, label = sprintf("mean=%.3f", mean_val)),
              hjust = 1.05, vjust = -0.4, size = 3, color = "gray30", inherit.aes = FALSE) +
    facet_wrap(~metric, ncol = 3, scales = "free_y") +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_number(accuracy = 0.01)) +
    scale_fill_brewer(palette = "Set2") +
    labs(title    = "Leave-One-Tree-Out Validation",
         subtitle = sprintf("RF classifier | %d annotated trees | %d trees RF",
                            length(params$loto_tree_ids), params$n_trees),
         x = "Held-out tree (test)", y = "Metric value", fill = "Held-out tree") +
    theme_minimal(base_size = 12) +
    theme(legend.position   = "none",
          strip.text        = element_text(face = "bold"),
          axis.text.x       = element_text(angle = 30, hjust = 1),
          panel.grid.major  = element_line(color = "gray90"),
          panel.grid.minor  = element_blank())

  ggsave(file.path(loto_output_dir, "00_LOTO_metrics.png"),
         p_loto, width = 12, height = 7, dpi = 150)

  cat(sprintf("\nLOTO outputs saved to: %s\n", loto_output_dir))
  invisible(df_res)
}

# 5. LEAF AREA ESTIMATION (MUNC) WITH ERROR PROPAGATION

process_classified_las <- function(path, params) {
  tree_name <- gsub("_RF$", "", file_path_sans_ext(basename(path)))
  cat(sprintf("\n--- %s ---\n", tree_name))
  t0 <- proc.time()

  outdir <- file.path(params$output_dir, tree_name)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  las <- readLAS(path, select = "*")
  if (is.null(las) || npoints(las) == 0) { cat("  ERROR: empty file\n"); return(NULL) }

  df    <- las@data
  xyz   <- cbind(df$X, df$Y, df$Z)
  n_pts <- npoints(las)

  # Read RF classification
  rf_field <- intersect(c("RF.class","RF_class"), names(df))
  if (length(rf_field) == 0)
    stop("RF.class field not found. Fields: ", paste(names(df), collapse=", "))
  class_rf <- as.integer(df[[rf_field[1]]])

  n_leaf <- sum(class_rf == params$class_leaf)
  n_wood <- sum(class_rf == params$class_wood)
  cat(sprintf("  Leaf: %d (%.1f%%) | Wood: %d (%.1f%%)\n",
              n_leaf, 100*n_leaf/n_pts, n_wood, 100*n_wood/n_pts))

  # RF probability (used as weight in MUNC)
  prob_field <- intersect(c("RF.prob","RF_prob"), names(df))
  rf_prob    <- if (length(prob_field) > 0) as.numeric(df[[prob_field[1]]]) else NULL

  # CANUPO agreement (if available)
  concordance <- NA_real_
  if ("CANUPO.class" %in% names(df)) {
    canupo      <- as.integer(df[["CANUPO.class"]])
    concordance <- round(100 * mean(class_rf == canupo, na.rm=TRUE), 1)
    cat(sprintf("  RF vs CANUPO agreement: %.1f%%\n", concordance))
  }

  idx_leaf <- which(class_rf == params$class_leaf)
  idx_wood <- which(class_rf == params$class_wood)
  if (length(idx_leaf) < 10) { cat("  WARN: fewer than 10 leaf points\n"); return(NULL) }

  prob_leaf <- if (!is.null(rf_prob)) rf_prob[idx_leaf] else NULL
  prob_wood <- if (!is.null(rf_prob)) rf_prob[idx_wood] else NULL
  prob_mean <- if (!is.null(prob_leaf)) round(mean(prob_leaf, na.rm=TRUE), 3) else NA

  # MUNC leaf area (weighted by RF.prob, after voxelisation)
  cat("  Computing normals (leaf)...\n")
  norm_leaf <- compute_normals(xyz[idx_leaf, ], k = params$k_normals)
  cat("  MUNC area (voxelised)...\n")
  munc <- munc_area(xyz[idx_leaf,], norm_leaf, params$voxel_size, weights = prob_leaf)
  cat(sprintf("  Voxels: %d / %d leaf points\n", munc$n_voxel, n_leaf))
  cat(sprintf("  Leaf area: %.4f m2\n", munc$area_m2))

  # FP/FN range (requires normals on wood points)
  fp_fn <- NULL
  if (!is.null(rf_prob) && length(idx_wood) >= 10) {
    cat("  Computing normals (wood, for FN estimate)...\n")
    norm_wood <- compute_normals(xyz[idx_wood,], k = params$k_normals)
    fp_fn <- fp_fn_range(xyz[idx_leaf,], norm_leaf, prob_leaf,
                          xyz[idx_wood,], norm_wood, prob_wood,
                          params$voxel_size)
    cat(sprintf("  FP (overestimate): %.4f m2\n", fp_fn$fp_m2))
    cat(sprintf("  FN (underestimate): %.4f m2\n", fp_fn$fn_m2))
  }

  # Geometric error propagation
  err_coreg <- get_coreg_err(tree_name, coreg)
  errors    <- propagate_errors(munc$area_m2, munc$n_voxel,
                                 params$voxel_size, params$k_normals, err_coreg)

  # Total range: geometric IC95% (symmetric) + FP/FN (asymmetric)
  if (!is.null(fp_fn)) {
    range_min <- munc$area_m2 - fp_fn$fp_m2 - errors$ic95_geo_m2
    range_max <- munc$area_m2 + fp_fn$fn_m2 + errors$ic95_geo_m2
  } else {
    range_min <- munc$area_m2 - errors$ic95_geo_m2
    range_max <- munc$area_m2 + errors$ic95_geo_m2
  }
  range_min <- max(0, range_min)

  cat(sprintf("  Errors: coreg=%.4f | voxel=%.4f | normals=%.4f | sigma_geo=%.4f m2\n",
              errors$sigma_coreg_m2, errors$sigma_voxel_m2,
              errors$sigma_normals_m2, errors$sigma_geo_m2))
  cat(sprintf("  Final range: [%.4f, %.4f] m2\n", range_min, range_max))
  cat(sprintf("  Done in %.1f s\n", (proc.time() - t0)["elapsed"]))

  data.frame(
    tree_name          = tree_name,
    n_pts              = n_pts,
    n_leaf             = n_leaf,
    n_wood             = n_wood,
    pct_leaf           = round(100 * n_leaf / n_pts, 1),
    prob_mean_leaf     = prob_mean,
    area_m2            = round(munc$area_m2, 4),
    se_geo_m2          = round(munc$se_m2,   4),
    cv_pct             = munc$cv_pct,
    n_voxel            = munc$n_voxel,
    fp_m2              = if (!is.null(fp_fn)) round(fp_fn$fp_m2, 4) else NA_real_,
    fn_m2              = if (!is.null(fp_fn)) round(fp_fn$fn_m2, 4) else NA_real_,
    sigma_coreg_m2     = errors$sigma_coreg_m2,
    sigma_voxel_m2     = errors$sigma_voxel_m2,
    sigma_normals_m2   = errors$sigma_normals_m2,
    sigma_geo_m2       = errors$sigma_geo_m2,
    err_coreg_input_m  = err_coreg,
    range_min_m2       = round(range_min, 4),
    range_max_m2       = round(range_max, 4),
    range_width_m2     = round(range_max - range_min, 4),
    concordance_pct    = concordance,
    elapsed_sec        = round((proc.time() - t0)["elapsed"], 1),
    status             = "OK",
    stringsAsFactors   = FALSE
  )
}

run_leaf_area_batch <- function(params) {
  files <- list.files(params$output_dir, pattern = "_RF\\.las$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0)
    stop("No *_RF.las files found in: ", params$output_dir)

  cat(sprintf("=== LEAF AREA BATCH (%d files) ===\n", length(files)))
  for (f in files) cat("  -", basename(f), "\n")
  cat(rep("=", 60), "\n", sep = "")

  summary_df <- data.frame()
  t_batch    <- proc.time()

  for (i in seq_along(files)) {
    cat(sprintf("\n[%d/%d] ", i, length(files)))
    riga <- tryCatch(
      process_classified_las(files[i], params),
      error = function(e) {
        cat("  ERROR:", conditionMessage(e), "\n")
        data.frame(
          tree_name = gsub("_RF$", "", file_path_sans_ext(basename(files[i]))),
          status    = paste("ERR:", conditionMessage(e)),
          stringsAsFactors = FALSE
        )
      }
    )
    if (!is.null(riga)) summary_df <- bind_rows(summary_df, riga)
    cat(rep("-", 60), "\n", sep = "")
  }

  write.csv2(summary_df,
             file.path(params$output_dir, "00_RF_summary.csv"),
             row.names = FALSE)

  # Summary plot: leaf area with error range per tree
  ok <- summary_df[!is.na(summary_df$range_min_m2) &
                     summary_df$status == "OK", ]
  if (nrow(ok) > 1) {
    p_range <- ggplot(ok, aes(x = reorder(tree_name, area_m2))) +
      geom_segment(aes(xend = reorder(tree_name, area_m2),
                       y = range_min_m2, yend = range_max_m2),
                   color = "gray60", linewidth = 3, alpha = 0.4) +
      geom_point(aes(y = area_m2), color = "#2d8a4e", size = 3.5) +
      geom_point(aes(y = range_min_m2), shape = "|", size = 5, color = "#E87722") +
      geom_point(aes(y = range_max_m2), shape = "|", size = 5, color = "#3B6FA0") +
      coord_flip() +
      labs(title    = "Leaf area per tree – full error range",
           subtitle = "green = estimate | orange = min (area - FP - 1.96*sigma_geo) | blue = max (area + FN + 1.96*sigma_geo)",
           x = "Tree", y = "Leaf area (m2)") +
      theme_minimal(base_size = 12)
    ggsave(file.path(params$output_dir, "00_compare_range_RF.png"),
           p_range, width = 10, height = 6, dpi = 150)
  }

  elapsed_tot <- (proc.time() - t_batch)["elapsed"]
  cat(sprintf("\n=== LEAF AREA DONE in %.1f min ===\n", elapsed_tot / 60))
  cat(sprintf("Trees OK: %d / %d\n",
              sum(summary_df$status == "OK", na.rm=TRUE), nrow(summary_df)))

  cols_print <- intersect(c("tree_name","pct_leaf","area_m2","sigma_geo_m2",
                             "fp_m2","fn_m2","range_min_m2","range_max_m2",
                             "range_width_m2","status"), names(summary_df))
  print(summary_df[, cols_print, drop=FALSE], row.names=FALSE)

  invisible(summary_df)
}

# 6. CROWN-PROJECTED LAI (ITSMe)
# LAI = MUNC-estimated leaf area / projected crown area (leaf points only)

compute_projected_area <- function(path, params, summary_rf) {
  tree_name <- gsub("_RF$", "", file_path_sans_ext(basename(path)))
  cat(sprintf("\n--- %s ---\n", tree_name))
  t0 <- proc.time()

  las <- readLAS(path, select = "*")
  if (is.null(las) || npoints(las) == 0) { cat("  ERROR: empty file\n"); return(NULL) }

  df    <- las@data
  xyz   <- data.frame(X = df$X, Y = df$Y, Z = df$Z)

  rf_field <- intersect(c("RF.class","RF_class"), names(df))
  if (length(rf_field) == 0)
    stop("RF.class field not found. Fields: ", paste(names(df), collapse=", "))
  class_rf <- as.integer(df[[rf_field[1]]])

  # Helper: extract numeric value from projected_area_pc output
  # (returns a 'units' object from sf::st_area or a list depending on ITSMe version)
  get_pa <- function(pc_df, concavity) {
    tryCatch({
      res <- projected_area_pc(pc = pc_df, concavity = concavity, plot = FALSE)
      val <- if (is.list(res)) res[[1]] else res
      as.numeric(val)
    }, error = function(e) {
      cat("  ERROR projected_area_pc:", e$message, "\n")
      NA_real_
    })
  }

  cat("  Projected area (all points)...\n")
  pa_all_m2  <- as.numeric(get_pa(xyz, params$concavity))
  cat(sprintf("  PA all:  %.4f m2\n", pa_all_m2))

  idx_leaf <- which(class_rf == params$class_leaf)
  pa_leaf_m2 <- NA_real_
  if (length(idx_leaf) >= 10) {
    cat("  Projected area (leaf points only)...\n")
    pa_leaf_m2 <- as.numeric(get_pa(xyz[idx_leaf, ], params$concavity))
    cat(sprintf("  PA leaf: %.4f m2\n", pa_leaf_m2))
  } else {
    cat("  WARN: too few leaf points for crown PA\n")
  }

  diff_m2  <- if (!is.na(pa_all_m2) && !is.na(pa_leaf_m2)) pa_all_m2 - pa_leaf_m2 else NA_real_
  diff_pct <- if (!is.na(diff_m2) && pa_all_m2 > 0) round(100 * diff_m2 / pa_all_m2, 1) else NA_real_

  # Retrieve leaf area from RF summary
  idx_rf    <- match(tree_name, summary_rf$tree_name)
  area_m2   <- as.numeric(if (!is.na(idx_rf)) summary_rf$area_m2[idx_rf]     else NA_real_)
  range_min <- as.numeric(if (!is.na(idx_rf)) summary_rf$range_min_m2[idx_rf] else NA_real_)
  range_max <- as.numeric(if (!is.na(idx_rf)) summary_rf$range_max_m2[idx_rf] else NA_real_)

  # Crown-projected LAI (leaf points only as denominator)
  lai        <- if (!is.na(area_m2) && !is.na(pa_leaf_m2) && pa_leaf_m2 > 0) round(area_m2 / pa_leaf_m2, 4) else NA_real_
  lai_min    <- if (!is.na(range_min) && !is.na(pa_leaf_m2) && pa_leaf_m2 > 0) round(range_min / pa_leaf_m2, 4) else NA_real_
  lai_max    <- if (!is.na(range_max) && !is.na(pa_leaf_m2) && pa_leaf_m2 > 0) round(range_max / pa_leaf_m2, 4) else NA_real_
  lai_all    <- if (!is.na(area_m2)   && !is.na(pa_all_m2)  && pa_all_m2 > 0)  round(area_m2  / pa_all_m2,  4) else NA_real_

  cat(sprintf("  Leaf area: %.4f m2\n", area_m2))
  cat(sprintf("  LAI (leaf PA): %.3f  [%.3f - %.3f]\n", lai, lai_min, lai_max))
  cat(sprintf("  LAI (all PA):  %.3f\n", lai_all))
  cat(sprintf("  Done in %.1f s\n", (proc.time() - t0)["elapsed"]))

  data.frame(
    tree_name        = tree_name,
    pa_all_m2        = round(pa_all_m2,  4),
    pa_leaf_m2       = round(pa_leaf_m2, 4),
    diff_pa_m2       = round(diff_m2,    4),
    diff_pa_pct      = diff_pct,
    leaf_area_m2     = round(area_m2,    4),
    leaf_area_min_m2 = round(range_min,  4),
    leaf_area_max_m2 = round(range_max,  4),
    lai              = lai,
    lai_min          = lai_min,
    lai_max          = lai_max,
    lai_all          = lai_all,
    elapsed_sec      = round((proc.time() - t0)["elapsed"], 1),
    status           = "OK",
    stringsAsFactors = FALSE
  )
}

run_lai_batch <- function(params) {
  # Load leaf area summary produced by run_leaf_area_batch()
  summary_rf_path <- file.path(params$output_dir, "00_RF_summary.csv")
  if (!file.exists(summary_rf_path))
    stop("00_RF_summary.csv not found in: ", params$output_dir,
         "\nRun run_leaf_area_batch() first.")

  summary_rf <- read.csv2(summary_rf_path, stringsAsFactors = FALSE)
  for (col in c("area_m2","range_min_m2","range_max_m2"))
    if (col %in% names(summary_rf))
      summary_rf[[col]] <- as.numeric(gsub(",", ".", summary_rf[[col]]))
  cat(sprintf("Loaded 00_RF_summary.csv: %d trees\n", nrow(summary_rf)))

  files <- list.files(params$output_dir, pattern = "_RF\\.las$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0)
    stop("No *_RF.las files found in: ", params$output_dir)

  cat(sprintf("=== LAI BATCH (%d files) ===\n", length(files)))
  lai_df  <- data.frame()
  t_batch <- proc.time()

  for (i in seq_along(files)) {
    cat(sprintf("\n[%d/%d] ", i, length(files)))
    riga <- tryCatch(
      compute_projected_area(files[i], params, summary_rf),
      error = function(e) {
        cat("  ERROR:", e$message, "\n")
        data.frame(tree_name = gsub("_RF$","",file_path_sans_ext(basename(files[i]))),
                   status = "ERROR", stringsAsFactors = FALSE)
      }
    )
    if (!is.null(riga)) lai_df <- bind_rows(lai_df, riga)
    cat(rep("-", 60), "\n", sep = "")
  }

  write.csv2(lai_df, file.path(params$output_dir, "00_LAI_summary.csv"),
             row.names = FALSE)

  ok <- lai_df[lai_df$status == "OK" & !is.na(lai_df$lai), ]
  if (nrow(ok) > 1) {

    # LAI per tree with uncertainty range
    p_lai <- ggplot(ok, aes(x = reorder(tree_name, lai))) +
      geom_segment(aes(xend = reorder(tree_name, lai),
                       y = lai_min, yend = lai_max),
                   color = "gray60", linewidth = 3, alpha = 0.4) +
      geom_point(aes(y = lai),     color = "#2d8a4e", size = 4) +
      geom_point(aes(y = lai_min), shape = "|", size = 5, color = "#E87722") +
      geom_point(aes(y = lai_max), shape = "|", size = 5, color = "#3B6FA0") +
      geom_point(aes(y = lai_all), color = "#8B4513", size = 2.5, shape = 17) +
      coord_flip() +
      labs(title    = "Crown-projected LAI per tree",
           subtitle = "green = LAI (leaf area / PA leaf) | brown triangle = LAI (leaf area / PA all)\norange/blue = uncertainty range from error propagation",
           x = "Tree", y = "LAI") +
      theme_minimal(base_size = 12)
    ggsave(file.path(params$output_dir, "00_LAI_range.png"),
           p_lai, width = 10, height = 6, dpi = 150)

    # Scatter: leaf area vs projected crown area
    p_scatter <- ggplot(ok, aes(x = pa_leaf_m2, y = leaf_area_m2, label = tree_name)) +
      geom_smooth(method = "lm", se = TRUE, color = "gray70", linetype = "dashed") +
      geom_point(aes(color = lai), size = 4) +
      geom_text(vjust = -0.8, size = 3.2) +
      scale_color_viridis_c(name = "LAI") +
      labs(title    = "Leaf area vs Projected crown area",
           subtitle = "Slope ≈ mean LAI | colour = individual LAI",
           x = "Projected crown area (m2)", y = "Leaf area (m2)") +
      theme_minimal(base_size = 12)
    ggsave(file.path(params$output_dir, "00_LAI_scatter.png"),
           p_scatter, width = 9, height = 7, dpi = 150)
  }

  elapsed_tot <- (proc.time() - t_batch)["elapsed"]
  cat(sprintf("\n=== LAI DONE in %.1f min ===\n", elapsed_tot / 60))
  cat(sprintf("Trees OK: %d / %d\n",
              sum(lai_df$status == "OK", na.rm=TRUE), nrow(lai_df)))

  cols_print <- intersect(c("tree_name","pa_all_m2","pa_leaf_m2","diff_pa_pct",
                             "leaf_area_m2","lai","lai_min","lai_max",
                             "lai_all","status"), names(lai_df))
  print(lai_df[, cols_print, drop=FALSE], row.names=FALSE)

  invisible(lai_df)
}

# =============================================================================
# 7. RUN PIPELINE
# Edit the params list in Section 1, then call these functions in order.
#
#   Step 1 – RF training + batch classification:
#     results <- run_classification_batch(params)
#
#   Step 1b – Classify additional trees with an existing model (optional):
#     recover_classification(params, model_path = "path/to/modello_RF.rds")
#
#   Step 2 – LOTO cross-validation:
#     loto_results <- run_loto(params)
#
#   Step 3 – Leaf area estimation (MUNC) + error propagation:
#     leaf_summary <- run_leaf_area_batch(params)
#
#   Step 4 – Crown-projected LAI:
#     lai_summary <- run_lai_batch(params)
# =============================================================================
