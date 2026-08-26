# =============================================================================
# RF_REPORT - Genera report e plot dai LAS già classificati con RF
# Legge i file *_RF.las nella cartella output e produce:
#   - report.txt per ogni albero
#   - plot diagnostici (top view, side view, profilo verticale)
#   - 00_RF_summary.csv completo
#   - 00_compare_area_RF.png
#   - 00_RF_feature_importance.png (se modello disponibile)
# =============================================================================

# =============================================================================
# 0. PACCHETTI
# =============================================================================
pacchetti <- c("lidR", "ranger", "ggplot2", "dplyr", "tidyr",
               "viridis", "tools", "progress", "FNN")

mancanti <- pacchetti[!pacchetti %in% installed.packages()[,"Package"]]
if (length(mancanti) > 0) {
  install.packages(mancanti, dependencies=TRUE)
}
invisible(lapply(pacchetti, library, character.only=TRUE))

# =============================================================================
# 1. PARAMETRI
# =============================================================================
params <- list(

  # Cartella con i *_RF.las classificati
  input_dir  = "C:/Users/emagi/OneDrive/Desktop/Petiole_XLA/TLS/Isolated_trees/SOR_per_classification/output_RF",

  # Cartella output report (puo coincidere con input_dir)
  output_dir = "C:/Users/emagi/OneDrive/Desktop/Petiole_XLA/TLS/Isolated_trees/SOR_per_classification/output_RF",

  # Modello RF salvato (per plot feature importance)
  modello_path = "C:/Users/emagi/OneDrive/Desktop/Petiole_XLA/TLS/Isolated_trees/SOR_per_classification/output_RF/modello_RF.rds",

  # Codici classe
  class_leaf = 1L,
  class_wood = 2L,

  # Calcolo area MUNC (per confronto con v7)
  k_normali  = 20,
  voxel_size = 0.007    # m - calibrato su BLK360
)

# =============================================================================
# 2. TABELLA CO-REGISTRAZIONE (Cyclone REGISTER 360 PLUS)
# Errore gruppo per ogni albero - inserire manualmente dai report PDF
# =============================================================================
coreg <- data.frame(
  tree_name   = c("SP_3957","SP_3958","SP_3961","SP_3962",
                  "SP_3969","SP_3970","SP_3975","SP_3976",
                  "SP_3982","SP_3983"),
  n_setup     = c(2L, 2L, 3L, 3L, 3L, 3L, 4L, 2L, 2L, 2L),
  err_coreg_m = c(0.023, 0.008, 0.008, 0.008,
                  0.009, 0.009, 0.007, 0.008,
                  0.008, 0.008),
  stringsAsFactors = FALSE
)
# Nota: SP_3957 ha errore 0.023 m (qualita bassa, solo 2 setup).
# Tutti gli altri rientrano nel range 0.007-0.009 m (qualita buona).

# Funzione per recuperare errore co-registrazione di un albero
get_err_coreg <- function(tree_name, coreg_table) {
  idx <- match(tree_name, coreg_table$tree_name)
  if (is.na(idx)) {
    warning("Errore co-registrazione non trovato per: ", tree_name,
            " -> uso valore medio 0.009 m")
    return(0.009)
  }
  coreg_table$err_coreg_m[idx]
}

# =============================================================================
# 3. PROPAGAZIONE ERRORI
# Stima analitica delle tre componenti di errore geometrico,
# poi propagazione in quadratura per ottenere sigma_geo totale.
#
# Le tre componenti:
#
# A) sigma_coreg (co-registrazione):
#    L'errore gruppo di Cyclone e' un RMS 3D della distanza residua
#    tra le scansioni sovrapposte. Per una superficie con normali
#    distribuite uniformemente, l'errore sull'area differenziale
#    da = voxel^2 / MUNC si propaga come:
#      sigma_coreg_area = area * (2/3) * (err_coreg / voxel_size)
#    Il fattore 2/3 deriva dalla media sull'angolo solido di superfici
#    con orientamento casuale (integrale di |cos(theta)| su sfera unitaria).
#
# B) sigma_voxel (voxelizzazione):
#    Il centro del voxel rappresentativo puo essere spostato fino a
#    voxel_size/2 rispetto alla superficie reale in direzione normale.
#    Questo introduce un'incertezza sull'orientamento della normale
#    locale, che si propaga su MUNC come:
#      sigma_voxel_area = area * (voxel_size / (2 * sqrt(3) * r_norm_medio))
#    dove r_norm_medio e' la distanza media al k-esimo vicino (proxy
#    della scala geometrica locale). Lo stimiamo come voxel_size * sqrt(k).
#
# C) sigma_normali (stima normali con PCA su k vicini):
#    L'incertezza sull'angolo della normale stimata da k punti rumorosi
#    e' approssimata da (Hoppe et al. 1992):
#      delta_theta ~ sigma_rumore / (lambda2 - lambda3)
#    dove sigma_rumore e' il rumore di posizione del singolo punto
#    (BLK360: ~0.006 m a 10 m) e (lambda2-lambda3) e' la separazione
#    degli autovalori (proxy della "planarità" locale).
#    In assenza di dati PCA per punto, approssimiamo come:
#      sigma_normali_area = area * sin(delta_theta_medio)
#    con delta_theta_medio stimato da sigma_rumore e k.
#
# =============================================================================
propaga_errori <- function(area_m2, munc_result, voxel_size, k_normali,
                            err_coreg_m,
                            sigma_rumore_scanner = 0.006) {
  # Propagazione corretta: ogni componente agisce sul singolo voxel,
  # i contributi dei voxel si sommano in quadratura (indipendenti).
  # Struttura comune: sigma = (area / sqrt(n_voxel)) * fattore
  # dove area/n_voxel = da_medio (contributo medio per voxel).
  #
  # Riferimento concettuale: Hoppe et al. 1992 per errore normali;
  # propagazione classica per misure indipendenti (Taylor 1997).

  n_voxel  <- munc_result$n_voxel
  da_medio <- area_m2 / (n_voxel + 1e-10)

  # Scala geometrica locale: distanza media ai k vicini su nuvola
  # voxelizzata, approssimata come voxel_size * sqrt(k_normali).
  # E' la lunghezza su cui la normale e' stimata -> denominatore
  # per convertire errore di posizione in errore angolare.
  r_locale <- voxel_size * sqrt(k_normali)

  # --- A) Co-registrazione ---
  # err_coreg_m e' l'errore di posizione 3D di ogni punto (RMS Cyclone).
  # Genera un errore angolare sulla normale: delta_theta_coreg = err_coreg / r_locale
  # Errore su da_i = da_i * delta_theta_coreg (per piccoli angoli)
  # Somma in quadratura su n_voxel contributi indipendenti:
  #   sigma_coreg = sqrt(n_voxel) * da_medio * delta_theta_coreg
  #              = area / sqrt(n_voxel) * (err_coreg / r_locale)
  delta_theta_coreg <- err_coreg_m / r_locale
  sigma_coreg <- (area_m2 / sqrt(n_voxel)) * delta_theta_coreg

  # --- B) Voxelizzazione ---
  # Il centro del voxel puo' essere spostato di voxel/2 rispetto alla
  # superficie reale. In direzione casuale rispetto alla normale,
  # lo spostamento efficace e' voxel / (2*sqrt(3)) (componente RMS
  # su 3 assi). Genera delta_theta_vox = (voxel/(2*sqrt(3))) / r_locale
  delta_theta_vox <- (voxel_size / (2 * sqrt(3))) / r_locale
  sigma_voxel <- (area_m2 / sqrt(n_voxel)) * delta_theta_vox

  # --- C) Stima normali (PCA su k vicini) ---
  # Il rumore di posizione del singolo punto (BLK360: ~0.006 m a 10 m)
  # genera un errore angolare sulla normale stimata con PCA su k punti:
  #   delta_theta_norm ~ sigma_scanner / r_locale
  # (approssimazione di primo ordine per superficie planare locale)
  delta_theta_norm <- sigma_rumore_scanner / r_locale
  delta_theta_norm <- min(delta_theta_norm, pi/6)  # cap a 30 gradi
  sigma_normali <- (area_m2 / sqrt(n_voxel)) * sin(delta_theta_norm)

  # --- Propagazione in quadratura ---
  sigma_geo <- sqrt(sigma_coreg^2 + sigma_voxel^2 + sigma_normali^2)

  list(
    sigma_coreg_m2   = round(sigma_coreg,  4),
    sigma_voxel_m2   = round(sigma_voxel,  4),
    sigma_normali_m2 = round(sigma_normali, 4),
    sigma_geo_m2     = round(sigma_geo,    4),
    ic95_geo_m2      = round(1.96 * sigma_geo, 4)
  )
}

# =============================================================================
# 4. FUNZIONI MUNC
# =============================================================================
calcola_normali <- function(xyz, k = 20) {
  n   <- nrow(xyz)
  N   <- matrix(0, n, 3)
  knn <- get.knn(xyz, k = k)
  pb  <- progress_bar$new(total=n, format="  [:bar] :percent ETA: :eta",
                           clear=FALSE)
  for (i in seq_len(n)) {
    pb$tick()
    vicini <- xyz[knn$nn.index[i,], , drop=FALSE]
    centro <- colMeans(vicini)
    D      <- sweep(vicini, 2, centro)
    eig    <- tryCatch(eigen(crossprod(D)/k, symmetric=TRUE),
                       error=function(e) NULL)
    if (!is.null(eig)) N[i,] <- eig$vectors[, 3]
  }
  N
}

voxelizza <- function(xyz, normali, pesi, voxel) {
  # Riduce la nuvola a un punto per voxel (media delle normali e dei pesi)
  # Rende il calcolo MUNC indipendente dalla densita locale di punti.
  # Fan et al. 2025 eq. 7: voxelizzazione prima del calcolo dell'area.
  ix <- floor(xyz[,1] / voxel)
  iy <- floor(xyz[,2] / voxel)
  iz <- floor(xyz[,3] / voxel)
  key <- paste(ix, iy, iz, sep="_")

  # Un rappresentante per voxel: media delle normali e dei pesi
  idx_rep <- tapply(seq_len(nrow(xyz)), key, function(i) i[1])
  idx_rep <- unlist(idx_rep)

  # Per le normali: media vettoriale per voxel
  N_out <- matrix(0, length(idx_rep), 3)
  p_out <- if (!is.null(pesi)) numeric(length(idx_rep)) else NULL

  voxel_keys <- unique(key)
  for (k in seq_along(voxel_keys)) {
    ii <- which(key == voxel_keys[k])
    N_media <- colMeans(normali[ii, , drop=FALSE])
    # Rinormalizza
    len_n <- sqrt(sum(N_media^2)) + 1e-10
    N_out[k,] <- N_media / len_n
    if (!is.null(pesi)) p_out[k] <- mean(pesi[ii], na.rm=TRUE)
  }

  list(normali=N_out, pesi=p_out, n_voxel=length(idx_rep))
}

area_MUNC_con_errore <- function(xyz, normali, voxel, pesi = NULL) {
  # Passo 1: voxelizza (rende risultato indipendente dalla densita)
  vox <- voxelizza(xyz, normali, pesi, voxel)
  N   <- vox$normali
  p   <- vox$pesi

  # Passo 2: MUNC su punti voxelizzati
  len   <- sqrt(rowSums(N^2)) + 1e-10
  N_hat <- N / len
  munc  <- apply(abs(N_hat), 1, max)
  munc[munc < 1e-6] <- 1e-6
  da    <- (voxel^2) / munc

  # Applica pesi (P(foglia) da RF.prob)
  if (!is.null(p)) da_pesata <- da * pmin(pmax(p, 0), 1) else da_pesata <- da

  n   <- length(da_pesata)
  tot <- sum(da_pesata)

  # Standard error: sd(contributi individuali) * sqrt(n)
  # Misura la variabilita geometrica interna della superficie stimata.
  # NON include errore di classificazione (quello e nei range FP/FN).
  se  <- sd(da_pesata) * sqrt(n)

  list(
    area_m2   = tot,
    se_m2     = se,
    ci95_low  = max(0, tot - 1.96*se),
    ci95_high = tot + 1.96*se,
    cv_pct    = round(100 * se / (tot + 1e-10), 2),
    n_voxel   = vox$n_voxel,
    n_punti   = length(munc)
  )
}

calcola_range_fp_fn <- function(xyz_leaf, normali_leaf, rf_prob_leaf,
                                 xyz_wood, normali_wood, rf_prob_wood,
                                 voxel) {
  # Calcola range asimmetrico dell'area:
  #   FP (falsi positivi): punti foglia che potrebbero essere legno
  #      -> area sovrastimata = MUNC pesato per (1 - P(foglia))
  #   FN (falsi negativi): punti legno che potrebbero essere foglia
  #      -> area mancante = MUNC pesato per P(foglia) sui punti legno
  #
  # Range credibile: [area_stimata - FP, area_stimata + FN]

  # FP: per i punti foglia, peso = (1 - P(foglia)) = P(legno)
  pesi_fp <- 1 - pmin(pmax(rf_prob_leaf, 0), 1)
  vox_fp  <- voxelizza(xyz_leaf, normali_leaf, pesi_fp, voxel)
  N_fp    <- vox_fp$normali
  munc_fp <- apply(abs(N_fp / (sqrt(rowSums(N_fp^2)) + 1e-10)), 1, max)
  munc_fp[munc_fp < 1e-6] <- 1e-6
  fp_m2   <- sum((voxel^2) / munc_fp * pmin(pmax(vox_fp$pesi, 0), 1))

  # FN: per i punti legno, P(foglia) = 1 - RF.prob (RF.prob e P(legno) qui)
  # Nota: RF.prob sui punti legno e P(foglia) dal modello RF
  # quindi usiamo direttamente rf_prob_wood come peso
  pesi_fn <- pmin(pmax(rf_prob_wood, 0), 1)
  fn_m2   <- 0
  if (!is.null(xyz_wood) && nrow(xyz_wood) >= 10) {
    vox_fn  <- voxelizza(xyz_wood, normali_wood, pesi_fn, voxel)
    N_fn    <- vox_fn$normali
    munc_fn <- apply(abs(N_fn / (sqrt(rowSums(N_fn^2)) + 1e-10)), 1, max)
    munc_fn[munc_fn < 1e-6] <- 1e-6
    fn_m2   <- sum((voxel^2) / munc_fn * pmin(pmax(vox_fn$pesi, 0), 1))
  }

  list(fp_m2 = fp_m2, fn_m2 = fn_m2)
}

# =============================================================================
# 3. PLOT DIAGNOSTICI
# =============================================================================
crea_plot <- function(xyz, classe, tree_name, outdir) {
  df <- data.frame(X=xyz[,1], Y=xyz[,2], Z=xyz[,3],
                   Classe=factor(classe,
                                 levels=c(params$class_wood, params$class_leaf),
                                 labels=c("Legno","Foglia")))
  colori <- c("Legno"="#8B4513", "Foglia"="#2d8a4e")

  p1 <- ggplot(df, aes(X, Y, color=Classe)) +
    geom_point(size=0.2, alpha=0.4) +
    scale_color_manual(values=colori) + coord_equal() +
    labs(title=paste0("Vista dall'alto - ", tree_name)) +
    theme_minimal(base_size=11)
  ggsave(file.path(outdir, "01_topview.png"), p1, width=8, height=8, dpi=150)

  p2 <- ggplot(df, aes(X, Z, color=Classe)) +
    geom_point(size=0.2, alpha=0.4) +
    scale_color_manual(values=colori) +
    labs(title=paste0("Vista laterale - ", tree_name),
         x="X (m)", y="Altezza (m)") +
    theme_minimal(base_size=11)
  ggsave(file.path(outdir, "02_sideview.png"), p2, width=10, height=6, dpi=150)

  df$fascia <- cut(df$Z, breaks=20)
  prof <- df %>%
    group_by(fascia, Classe) %>%
    summarise(n=n(), .groups="drop") %>%
    group_by(fascia) %>%
    mutate(pct=100*n/sum(n))

  p3 <- ggplot(prof, aes(x=fascia, y=pct, fill=Classe)) +
    geom_col() +
    scale_fill_manual(values=colori) + coord_flip() +
    labs(title=paste0("Profilo verticale - ", tree_name),
         x="Fascia altezza", y="% punti") +
    theme_minimal(base_size=10)
  ggsave(file.path(outdir, "03_profilo.png"), p3, width=8, height=6, dpi=150)
}

# =============================================================================
# 4. REPORT TESTUALE
# =============================================================================
scrivi_report <- function(tree_name, n_pts, n_leaf, n_wood,
                           munc, fp_fn, errori, err_coreg_m,
                           range_tot_min, range_tot_max,
                           prob_media, concordanza_canupo,
                           params, outdir) {

  righe_classif <- if (!is.null(fp_fn)) {
    c("",
      "--- ERRORE DI CLASSIFICAZIONE RF ---",
      "  (componente asimmetrica: FP e FN pesati per RF.prob)",
      sprintf("  FP area sovrastimata: %.4f m2  (punti foglia mal classificati)",
              fp_fn$fp_m2),
      sprintf("  FN area mancante:     %.4f m2  (punti legno mal classificati)",
              fp_fn$fn_m2)
    )
  } else character(0)

  righe_geo <- c(
    "",
    "--- ERRORI GEOMETRICI (propagazione in quadratura) ---",
    sprintf("  Co-registrazione (Cyclone): %.4f m2  [err_gruppo=%.3f m]",
            errori$sigma_coreg_m2, err_coreg_m),
    sprintf("  Voxelizzazione:             %.4f m2  [voxel=%.3f m]",
            errori$sigma_voxel_m2, params$voxel_size),
    sprintf("  Stima normali (PCA k=%d):   %.4f m2",
            params$k_normali, errori$sigma_normali_m2),
    sprintf("  SIGMA_GEO totale:           %.4f m2",
            errori$sigma_geo_m2),
    sprintf("  IC 95%% geometrico:       +/- %.4f m2",
            errori$ic95_geo_m2)
  )

  righe_range <- c(
    "",
    "--- RANGE FINALE AREA FOGLIARE ---",
    "  (errore classif. asimmetrico + IC95% geometrico)",
    sprintf("  Minimo:   %.4f m2  (area - FP - 1.96*sigma_geo)",
            range_tot_min),
    sprintf("  Centrale: %.4f m2  (stima puntuale)", munc$area_m2),
    sprintf("  Massimo:  %.4f m2  (area + FN + 1.96*sigma_geo)",
            range_tot_max),
    sprintf("  Ampiezza range: %.4f m2", range_tot_max - range_tot_min)
  )

  righe <- c(
    "=======================================================",
    "  RF CLASSIFY - REPORT",
    "=======================================================",
    paste("Albero:      ", tree_name),
    paste("Data analisi:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "--- CLASSIFICAZIONE RANDOM FOREST ---",
    sprintf("  Punti totali:       %d", n_pts),
    sprintf("  Foglia:             %d (%.1f%%)", n_leaf, 100*n_leaf/n_pts),
    sprintf("  Legno:              %d (%.1f%%)", n_wood, 100*n_wood/n_pts),
    sprintf("  P(foglia) media:    %.3f",
            if (!is.na(prob_media)) prob_media else NA),
    if (!is.na(concordanza_canupo))
      sprintf("  Concordanza CANUPO: %.1f%%", concordanza_canupo)
    else
      "  Concordanza CANUPO: non disponibile",
    "",
    "--- AREA FOGLIARE (MUNC + voxelizzazione) ---",
    "  Metodo: pesata per RF.prob | voxelizzata (indip. da densita punti)",
    sprintf("  Voxel usati: %d / %d punti foglia",
            munc$n_voxel, n_leaf),
    "",
    sprintf("  Area stimata:  %.4f m2", munc$area_m2),
    sprintf("  SE geometrico: %.4f m2  (variab. geometrica interna)",
            munc$se_m2),
    sprintf("  CV:            %.2f%%", munc$cv_pct),
    righe_classif,
    righe_geo,
    righe_range,
    "",
    "--- PARAMETRI ---",
    paste("  k normali:   ", params$k_normali),
    paste("  Voxel size:  ", params$voxel_size, "m"),
    paste("  Err. co-reg.:", err_coreg_m, "m (Cyclone REGISTER 360 PLUS)"),
    "======================================================="
  )
  writeLines(righe, file.path(outdir, "report.txt"))
  cat("  Report salvato.\n")
}

# =============================================================================
# 5. PROCESSA SINGOLO LAS CLASSIFICATO
# =============================================================================
processa_rf_las <- function(percorso, params) {
  # Nome albero: rimuove suffisso _RF
  tree_name <- gsub("_RF$", "",
                    file_path_sans_ext(basename(percorso)))
  cat(sprintf("\n--- %s ---\n", tree_name))
  t0 <- proc.time()

  outdir <- file.path(params$output_dir, tree_name)
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)

  # Carica LAS
  las <- readLAS(percorso, select="*")
  if (is.null(las) || npoints(las)==0) {
    cat("  ERRORE: file vuoto\n"); return(NULL)
  }
  cat(sprintf("  Punti totali: %d\n", npoints(las)))

  df    <- las@data
  xyz   <- cbind(df$X, df$Y, df$Z)
  n_pts <- npoints(las)

  # Leggi classificazione RF
  rf_field <- intersect(c("RF.class","RF_class"), names(df))
  if (length(rf_field)==0)
    stop("Campo RF.class non trovato. Campi presenti: ",
         paste(names(df), collapse=", "))

  classe <- as.integer(df[[rf_field[1]]])
  cat(sprintf("  Foglia: %d (%.1f%%) | Legno: %d (%.1f%%)\n",
              sum(classe==params$class_leaf),
              100*sum(classe==params$class_leaf)/n_pts,
              sum(classe==params$class_wood),
              100*sum(classe==params$class_wood)/n_pts))

  # Leggi RF.prob per usarla come peso MUNC
  prob_field <- intersect(c("RF.prob","RF_prob"), names(df))
  rf_prob    <- if (length(prob_field)>0) as.numeric(df[[prob_field[1]]]) else NULL

  # Leggi CANUPO.class per concordanza
  concordanza <- NA_real_
  if ("CANUPO.class" %in% names(df)) {
    canupo      <- as.integer(df[["CANUPO.class"]])
    concordanza <- round(100*mean(classe==canupo, na.rm=TRUE), 1)
    cat(sprintf("  Concordanza RF vs CANUPO: %.1f%%\n", concordanza))
  }

  # Indici foglia e legno
  idx_leaf <- which(classe == params$class_leaf)
  idx_wood <- which(classe == params$class_wood)
  if (length(idx_leaf) < 10) {
    cat("  WARN: meno di 10 punti foglia\n"); return(NULL)
  }

  # RF.prob sui punti foglia (peso MUNC) e legno (stima FN)
  pesi_leaf  <- if (!is.null(rf_prob)) rf_prob[idx_leaf] else NULL
  pesi_wood  <- if (!is.null(rf_prob)) rf_prob[idx_wood] else NULL
  prob_media <- if (!is.null(pesi_leaf)) round(mean(pesi_leaf, na.rm=TRUE), 3) else NA

  # Normali punti foglia
  cat("  Calcolo normali (foglia)...\n")
  normali_leaf <- calcola_normali(xyz[idx_leaf,], k=params$k_normali)

  # Area MUNC con voxelizzazione (indipendente dalla densita)
  cat("  Calcolo area MUNC (con voxelizzazione)...\n")
  munc <- area_MUNC_con_errore(xyz[idx_leaf,], normali_leaf,
                                params$voxel_size, pesi=pesi_leaf)
  cat(sprintf("  Voxel usati: %d / %d punti foglia\n",
              munc$n_voxel, length(idx_leaf)))
  cat(sprintf("  Area stimata: %.4f m2\n", munc$area_m2))

  # Range FP/FN (normali sui punti legno necessarie per FN)
  fp_fn <- NULL
  if (!is.null(rf_prob) && length(idx_wood) >= 10) {
    cat("  Calcolo normali (legno per stima FN)...\n")
    normali_wood <- calcola_normali(xyz[idx_wood,], k=params$k_normali)
    fp_fn <- calcola_range_fp_fn(
      xyz[idx_leaf,], normali_leaf, pesi_leaf,
      xyz[idx_wood,], normali_wood, pesi_wood,
      params$voxel_size
    )
    cat(sprintf("  FP (area sovrastimata): %.4f m2\n", fp_fn$fp_m2))
    cat(sprintf("  FN (area mancante):     %.4f m2\n", fp_fn$fn_m2))
    cat(sprintf("  Range credibile: [%.4f, %.4f] m2\n",
                munc$area_m2 - fp_fn$fp_m2,
                munc$area_m2 + fp_fn$fn_m2))
  }

  # Propagazione errori geometrici
  err_coreg <- get_err_coreg(tree_name, coreg)
  errori    <- propaga_errori(
    area_m2    = munc$area_m2,
    munc_result = munc,
    voxel_size  = params$voxel_size,
    k_normali   = params$k_normali,
    err_coreg_m = err_coreg
  )

  # Range finale combinato:
  # - componente geometrica (simmetrica): +/- 1.96 * sigma_geo
  # - componente classificazione (asimmetrica): -FP / +FN
  # Il range totale e' la somma delle due componenti
  if (!is.null(fp_fn)) {
    range_tot_min <- munc$area_m2 - fp_fn$fp_m2 - errori$ic95_geo_m2
    range_tot_max <- munc$area_m2 + fp_fn$fn_m2 + errori$ic95_geo_m2
  } else {
    range_tot_min <- munc$area_m2 - errori$ic95_geo_m2
    range_tot_max <- munc$area_m2 + errori$ic95_geo_m2
  }
  range_tot_min <- max(0, range_tot_min)

  cat(sprintf("  Errori: coreg=%.4f | voxel=%.4f | normali=%.4f | sigma_geo=%.4f m2\n",
              errori$sigma_coreg_m2, errori$sigma_voxel_m2,
              errori$sigma_normali_m2, errori$sigma_geo_m2))
  cat(sprintf("  Range finale: [%.4f, %.4f] m2\n",
              range_tot_min, range_tot_max))

  # Plot
  crea_plot(xyz, classe, tree_name, outdir)

  # Report
  scrivi_report(tree_name,
                n_pts              = n_pts,
                n_leaf             = length(idx_leaf),
                n_wood             = length(idx_wood),
                munc               = munc,
                fp_fn              = fp_fn,
                errori             = errori,
                err_coreg_m        = err_coreg,
                range_tot_min      = range_tot_min,
                range_tot_max      = range_tot_max,
                prob_media         = prob_media,
                concordanza_canupo = concordanza,
                params             = params,
                outdir             = outdir)

  elapsed <- (proc.time()-t0)["elapsed"]
  cat(sprintf("  Completato in %.1f s\n", elapsed))

  data.frame(
    tree_name              = tree_name,
    n_pts                  = n_pts,
    n_leaf                 = length(idx_leaf),
    n_wood                 = length(idx_wood),
    pct_leaf               = round(100*length(idx_leaf)/n_pts, 1),
    prob_media_leaf        = prob_media,
    # Area stimata (m2)
    area_m2                = round(munc$area_m2, 4),
    se_geo_m2              = round(munc$se_m2, 4),
    cv_pct                 = munc$cv_pct,
    n_voxel                = munc$n_voxel,
    # Errore classificazione RF (asimmetrico)
    fp_m2                  = if (!is.null(fp_fn)) round(fp_fn$fp_m2, 4) else NA_real_,
    fn_m2                  = if (!is.null(fp_fn)) round(fp_fn$fn_m2, 4) else NA_real_,
    # Errori geometrici (propagazione in quadratura)
    sigma_coreg_m2         = errori$sigma_coreg_m2,
    sigma_voxel_m2         = errori$sigma_voxel_m2,
    sigma_normali_m2       = errori$sigma_normali_m2,
    sigma_geo_m2           = errori$sigma_geo_m2,
    err_coreg_input_m      = err_coreg,
    # Range finale (m2): classif. asimmetrico + IC95% geometrico
    range_min_m2           = round(range_tot_min, 4),
    range_max_m2           = round(range_tot_max, 4),
    range_ampiezza_m2      = round(range_tot_max - range_tot_min, 4),
    concordanza_canupo_pct = concordanza,
    elapsed_sec            = round(elapsed, 1),
    status                 = "OK",
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 6. BATCH
# =============================================================================
files <- list.files(params$input_dir, pattern="_RF\\.las$",
                    full.names=TRUE, ignore.case=TRUE)
if (length(files)==0)
  stop("Nessun file *_RF.las trovato in: ", params$input_dir)

cat(sprintf("=== RF REPORT BATCH (%d file) ===\n", length(files)))
for (f in files) cat("  -", basename(f), "\n")
cat(rep("=", 60), "\n", sep="")

summary_df <- data.frame()
t_batch    <- proc.time()

for (i in seq_along(files)) {
  cat(sprintf("\n[%d/%d] ", i, length(files)))
  riga <- tryCatch(
    processa_rf_las(files[i], params),
    error=function(e) {
      cat("  ERRORE:", conditionMessage(e), "\n")
      data.frame(
        tree_name = gsub("_RF$", "", file_path_sans_ext(basename(files[i]))),
        status    = paste("ERR:", conditionMessage(e)),
        stringsAsFactors = FALSE
      )
    }
  )
  if (!is.null(riga)) summary_df <- dplyr::bind_rows(summary_df, riga)
  cat(rep("-", 60), "\n", sep="")
}

# Salva CSV
write.csv2(summary_df,
           file.path(params$output_dir, "00_RF_summary.csv"),
           row.names=FALSE)

# Plot comparativi
ok <- summary_df[!is.na(summary_df$area_cm2) & summary_df$status=="OK", ]

if (nrow(ok) > 1) {

  # Area fogliare con IC 95%
  p_area <- ggplot(ok, aes(x=reorder(tree_name, area_m2),
                            y=area_m2, fill=area_m2)) +
    geom_col(show.legend=FALSE) +
    geom_errorbar(aes(ymin=ci95_low_m2, ymax=ci95_high_m2),
                  width=0.3, color="gray30") +
    scale_fill_viridis() + coord_flip() +
    labs(title="Area fogliare RF per albero (IC 95%)",
         x="Albero", y="Area (m2)") +
    theme_minimal(base_size=12)
  ggsave(file.path(params$output_dir, "00_compare_area_RF.png"),
         p_area, width=9, height=6, dpi=150)

  # Range credibile FP/FN
  ok_range <- ok[!is.na(ok$range_min_m2) & !is.na(ok$range_max_m2), ]
  if (nrow(ok_range) > 1) {
    p_range <- ggplot(ok_range, aes(x=reorder(tree_name, area_m2))) +
      geom_segment(aes(xend=reorder(tree_name, area_m2),
                       y=range_min_m2, yend=range_max_m2),
                   color="gray60", linewidth=3, alpha=0.4) +
      geom_point(aes(y=area_m2), color="#2d8a4e", size=3.5) +
      geom_point(aes(y=range_min_m2), shape="|", size=5, color="#E87722") +
      geom_point(aes(y=range_max_m2), shape="|", size=5, color="#3B6FA0") +
      coord_flip() +
      labs(title="Range area fogliare RF - propagazione errori completa",
           subtitle="verde = stima | arancio = min (area - FP - 1.96*sigma_geo) | blu = max (area + FN + 1.96*sigma_geo)",
           x="Albero", y="Area (m2)") +
      theme_minimal(base_size=12)
    ggsave(file.path(params$output_dir, "00_compare_range_RF.png"),
           p_range, width=10, height=6, dpi=150)

    # Plot decomposizione errori per albero
    err_long <- ok_range %>%
      select(tree_name, area_m2, sigma_coreg_m2, sigma_voxel_m2,
             sigma_normali_m2, fp_m2, fn_m2) %>%
      pivot_longer(cols=c(sigma_coreg_m2, sigma_voxel_m2, sigma_normali_m2,
                          fp_m2, fn_m2),
                   names_to="componente", values_to="errore_m2") %>%
      mutate(componente = recode(componente,
        "sigma_coreg_m2"   = "Co-registrazione",
        "sigma_voxel_m2"   = "Voxelizzazione",
        "sigma_normali_m2" = "Stima normali",
        "fp_m2"            = "FP (classif.)",
        "fn_m2"            = "FN (classif.)"
      ))

    p_err_decomp <- ggplot(err_long,
                           aes(x=reorder(tree_name, area_m2),
                               y=errore_m2, fill=componente)) +
      geom_col(position="stack") +
      scale_fill_manual(values=c(
        "Co-registrazione" = "#E87722",
        "Voxelizzazione"   = "#F2C14E",
        "Stima normali"    = "#A8D5BA",
        "FP (classif.)"    = "#3B6FA0",
        "FN (classif.)"    = "#8B4513"
      )) +
      coord_flip() +
      labs(title="Decomposizione errori per albero",
           subtitle="Barre impilate: contributo di ogni componente all'incertezza totale",
           x="Albero", y="Errore (m2)", fill="Componente") +
      theme_minimal(base_size=12) +
      theme(legend.position="bottom")
    ggsave(file.path(params$output_dir, "00_errori_decomposizione.png"),
           p_err_decomp, width=10, height=6, dpi=150)
  }

  # Concordanza RF vs CANUPO
  ok_conc <- ok[!is.na(ok$concordanza_canupo_pct), ]
  if (nrow(ok_conc) > 1) {
    p_conc <- ggplot(ok_conc,
                     aes(x=reorder(tree_name, concordanza_canupo_pct),
                         y=concordanza_canupo_pct,
                         fill=concordanza_canupo_pct)) +
      geom_col(show.legend=FALSE) +
      scale_fill_viridis(option="magma") + coord_flip() +
      labs(title="Concordanza RF vs CANUPO per albero",
           subtitle="100% = classificazione identica",
           x="Albero", y="Concordanza (%)") +
      theme_minimal(base_size=12)
    ggsave(file.path(params$output_dir, "00_RF_vs_CANUPO.png"),
           p_conc, width=9, height=6, dpi=150)
  }

  # Heatmap
  heat_cols <- intersect(c("pct_leaf","area_m2","prob_media_leaf",
                            "concordanza_canupo_pct"), names(ok))
  heat_data <- ok %>%
    select(tree_name, all_of(heat_cols)) %>%
    pivot_longer(-tree_name, names_to="metrica", values_to="valore") %>%
    group_by(metrica) %>%
    mutate(valore_norm=(valore-min(valore,na.rm=TRUE))/
             (max(valore,na.rm=TRUE)-min(valore,na.rm=TRUE)+1e-10)) %>%
    ungroup()

  p_heat <- ggplot(heat_data, aes(x=metrica, y=tree_name, fill=valore_norm)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(valore,1)), size=2.8, color="white") +
    scale_fill_viridis(name="Norm.") +
    labs(title="Heatmap RF", x="", y="Albero") +
    theme_minimal(base_size=10) +
    theme(axis.text.x=element_text(angle=30, hjust=1))
  ggsave(file.path(params$output_dir, "00_heatmap_RF.png"),
         p_heat, width=10,
         height=max(4, nrow(ok)*0.5+2), dpi=150)
}

# Feature importance (se modello disponibile)
if (file.exists(params$modello_path)) {
  rf_obj <- readRDS(params$modello_path)
  imp    <- sort(rf_obj$modello$variable.importance, decreasing=TRUE)
  df_imp <- data.frame(feature=names(imp), importanza=imp,
                       stringsAsFactors=FALSE)
  df_imp$feature <- factor(df_imp$feature, levels=rev(df_imp$feature))
  df_imp$tipo <- dplyr::case_when(
    grepl("^lin|^plan|^sfer|^curv|^dens", df_imp$feature) ~ "Geometrica",
    grepl("^R$|^G$|^B$|^ExG",            df_imp$feature) ~ "RGB",
    grepl("intensity",                    df_imp$feature) ~ "Intensita",
    grepl("return",                       df_imp$feature) ~ "Ritorni",
    grepl("Z_norm",                       df_imp$feature) ~ "Altezza",
    TRUE                                                   ~ "Altro"
  )
  p_imp <- ggplot(df_imp, aes(x=feature, y=importanza, fill=tipo)) +
    geom_col() +
    scale_fill_viridis(discrete=TRUE) + coord_flip() +
    labs(title="Importanza feature - Random Forest",
         subtitle="Permutation importance",
         x="", y="Importanza", fill="Tipo") +
    theme_minimal(base_size=11)
  ggsave(file.path(params$output_dir, "00_RF_feature_importance.png"),
         p_imp, width=10,
         height=max(5, nrow(df_imp)*0.25+2), dpi=150)
  cat("  Feature importance salvata.\n")
}

elapsed_tot <- (proc.time()-t_batch)["elapsed"]
cat(sprintf("\n=== COMPLETATO in %.1f min ===\n", elapsed_tot/60))
cat(sprintf("Piante OK: %d / %d\n",
            sum(summary_df$status=="OK", na.rm=TRUE), nrow(summary_df)))
cat("Output:", params$output_dir, "\n\n")

cols_print <- intersect(c("tree_name","pct_leaf","area_m2",
                           "sigma_geo_m2","fp_m2","fn_m2",
                           "range_min_m2","range_max_m2","range_ampiezza_m2",
                           "concordanza_canupo_pct","status"), names(summary_df))
print(summary_df[, cols_print, drop=FALSE], row.names=FALSE)
