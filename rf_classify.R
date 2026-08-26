# =============================================================================
# RF_CLASSIFY - Random Forest per classificazione foglia/legno
# Alternativa a CANUPO per nuvole di punti TLS BLK360
#
# WORKFLOW:
#   1. Legge i punti di training (legno.las + foglia.las)
#   2. Calcola feature geometriche multiscala + radiometriche
#   3. Addestra Random Forest con cross-validazione 5-fold
#   4. Classifica le 10 piante e salva LAS con RF.class
#   5. Report comparativo RF vs CANUPO (se disponibile)
#
# INPUT:
#   - training_legno_path: .las con punti legno annotati
#   - training_foglia_path: .las con punti foglia annotati
#   - input_dir: cartella con i .las delle 10 piante
#
# OUTPUT PER PIANTA:
#   - NomeAlbero_RF.las con campo RF.class (1=foglia, 2=legno)
#     e RF.prob (probabilita P(foglia))
#
# OUTPUT BATCH:
#   - 00_RF_summary.csv
#   - 00_RF_feature_importance.png
#   - 00_RF_confusion_matrix.png
#   - 00_RF_vs_CANUPO.png (se i LAS hanno anche CANUPO.class)
# =============================================================================

# =============================================================================
# 0. PACCHETTI
# =============================================================================
pacchetti <- c("lidR", "ranger", "caret", "FNN", "ggplot2",
               "dplyr", "tidyr", "viridis", "tools", "progress")

mancanti <- pacchetti[!pacchetti %in% installed.packages()[,"Package"]]
if (length(mancanti) > 0) {
  cat("Installazione:", paste(mancanti, collapse=", "), "\n")
  install.packages(mancanti, dependencies=TRUE)
}
invisible(lapply(pacchetti, library, character.only=TRUE))

# =============================================================================
# 1. PARAMETRI
# =============================================================================
params <- list(

  # --- Percorsi training ---
  training_legno_path  = "C:/Users/emagi/OneDrive/Desktop/Petiole_XLA/TLS/Isolated_trees/SOR_per_classification/wood_training.las",
  training_foglia_path = "C:/Users/emagi/OneDrive/Desktop/Petiole_XLA/TLS/Isolated_trees/SOR_per_classification/leaf_training.las",

  # --- Piante da classificare ---
  input_dir  = "C:/Users/emagi/OneDrive/Desktop/Petiole_XLA/TLS/Isolated_trees/SOR_per_classification/classified",
  output_dir = "C:/Users/emagi/OneDrive/Desktop/Petiole_XLA/TLS/Isolated_trees/SOR_per_classification/output_RF",

  # --- Feature geometriche multiscala ---
  # Scale in metri: da piccola (dettaglio foglia) a grande (ramo)
  # Min ~3x risoluzione BLK360 (~0.02m), Max ~diametro ramo secondario
  scale_geometria = c(0.05, 0.10, 0.20, 0.40),

  # --- Feature da includere ---
  usa_rgb        = TRUE,   # R, G, B normalizzati
  usa_intensity  = TRUE,   # Intensita ritorno laser
  usa_return_num = TRUE,   # Return Number + Number of Returns + ratio
                           # (BLK360: "Return Number" e "Number of Returns")
  usa_altezza    = TRUE,   # Z normalizzato (gia presente nel LAS)
  usa_densita    = TRUE,   # Densita locale (punti in sfera r=0.1m)

  # --- Random Forest ---
  n_trees        = 500,    # alberi nel forest (500 e stabile, 1000 piu preciso)
  n_fold_cv      = 5,      # fold per cross-validazione
  max_train_pts  = 50000,  # max punti per classe nel training
                           # (sottocampiona se necessario per velocita)
  seed           = 42,

  # --- Soglia probabilita per classificazione ---
  # P(foglia) >= soglia -> foglia, altrimenti legno
  # 0.5 = massima accuratezza | >0.5 = piu conservativo sulle foglie
  soglia_foglia  = 0.50,

  # Codici classe output
  class_leaf = 1L,
  class_wood = 2L
)

if (!dir.exists(params$output_dir))
  dir.create(params$output_dir, recursive=TRUE)

cat("=== RF CLASSIFY - Random Forest per TLS ===\n")
cat("Scale geometriche:", paste(params$scale_geometria, collapse=", "), "m\n\n")

# =============================================================================
# 2. ESTRAZIONE FEATURE
# Calcola feature geometriche + radiometriche per ogni punto
# =============================================================================
estrai_feature <- function(las, scale_geom, params, label="") {
  cat(sprintf("  Estrazione feature%s (%d punti)...\n",
              label, npoints(las)))

  df   <- las@data
  xyz  <- cbind(df$X, df$Y, df$Z)
  n    <- nrow(df)

  # ---- Feature radiometriche (dirette dal LAS) ----
  feat <- data.frame(row.names=seq_len(n))

  if (params$usa_altezza)
    feat$Z_norm <- df$Z  # gia normalizzato

  if (params$usa_intensity && "Intensity" %in% names(df))
    feat$intensity <- as.numeric(df$Intensity)

  # lidR legge i nomi dei campi sostituendo spazi con punti
  # "Return Number"    -> ReturnNumber o Return.Number
  # "Number of Returns"-> NumberOfReturns o Number.of.Returns
  # Proviamo entrambe le varianti
  rn_field  <- intersect(c("ReturnNumber", "Return.Number",
                            "Return_Number"), names(df))
  nor_field <- intersect(c("NumberOfReturns", "Number.of.Returns",
                            "Number_of_Returns", "NumberReturns"), names(df))

  if (params$usa_return_num && length(rn_field) > 0)
    feat$return_num <- as.numeric(df[[rn_field[1]]])

  # Number of Returns: utile perche rami danno spesso ritorni multipli
  if (params$usa_return_num && length(nor_field) > 0)
    feat$n_returns <- as.numeric(df[[nor_field[1]]])

  # Ratio ritorno (Return Number / Number of Returns)
  # 1.0 = unico ritorno | <1.0 = primo di piu ritorni
  if (!is.null(feat$return_num) && !is.null(feat$n_returns))
    feat$return_ratio <- feat$return_num / (feat$n_returns + 1e-10)

  if (params$usa_rgb) {
    if (all(c("R","G","B") %in% names(df))) {
      # Normalizza in 0-1 (BLK360 usa 16 bit)
      mx <- if (max(df$R, na.rm=TRUE) > 255) 65535 else 255
      feat$R <- as.numeric(df$R) / mx
      feat$G <- as.numeric(df$G) / mx
      feat$B <- as.numeric(df$B) / mx
      # ExG (Excess Green) - utile per foglie verdi
      feat$ExG <- 2*feat$G - feat$R - feat$B
    } else {
      cat("  WARN: RGB non trovato nel LAS\n")
    }
  }

  # ---- Feature geometriche multiscala ----
  # Per ogni scala r: calcola PCA sui k vicini in sfera di raggio r
  # -> eigenvalues normalizzati -> linearita, planarità, sfericità
  cat("  Feature geometriche multiscala")
  for (r in scale_geom) {
    cat(sprintf(" [%.2fm]", r))

    # Numero vicini nella sfera di raggio r (stima adattiva)
    # Usiamo knn con k fisso = 20, poi filtriamo per distanza
    k <- min(30, n - 1)
    knn <- get.knn(xyz, k=k)

    lin_v  <- numeric(n)
    plan_v <- numeric(n)
    sfer_v <- numeric(n)
    curv_v <- numeric(n)
    dens_v <- numeric(n)

    pb <- progress_bar$new(total=n, format="  [:bar] :percent",
                           clear=TRUE, width=40)
    for (i in seq_len(n)) {
      pb$tick()
      # Prendi vicini entro raggio r
      dists <- knn$nn.dist[i,]
      idx   <- knn$nn.index[i, dists <= r]

      if (length(idx) < 3) {
        # Troppo pochi vicini: valori di default (punto isolato)
        lin_v[i] <- 1; next
      }

      pts    <- xyz[idx, , drop=FALSE]
      centro <- colMeans(pts)
      D      <- sweep(pts, 2, centro)
      eig    <- tryCatch(
        eigen(crossprod(D)/nrow(D), symmetric=TRUE),
        error=function(e) NULL
      )
      if (is.null(eig)) { lin_v[i] <- 1; next }

      lam <- sort(eig$values, decreasing=TRUE)
      lam <- pmax(lam, 0)
      tot <- sum(lam) + 1e-10

      # Feature dimensionali (Demantke et al. 2011)
      lin_v[i]  <- (lam[1] - lam[2]) / tot          # linearita
      plan_v[i] <- (lam[2] - lam[3]) / tot          # planarità
      sfer_v[i] <- lam[3] / tot                      # sfericità
      curv_v[i] <- lam[3] / (lam[1] + 1e-10)        # curvatura
      dens_v[i] <- length(idx) / (4/3 * pi * r^3)   # densita locale
    }

    sufx <- gsub("\\.", "", sprintf("_r%03d", round(r*100)))
    feat[[paste0("lin",  sufx)]]  <- lin_v
    feat[[paste0("plan", sufx)]]  <- plan_v
    feat[[paste0("sfer", sufx)]]  <- sfer_v
    feat[[paste0("curv", sufx)]]  <- curv_v
    if (params$usa_densita)
      feat[[paste0("dens", sufx)]] <- dens_v
  }
  cat("\n")

  feat
}

# =============================================================================
# 3. CARICA TRAINING SET
# =============================================================================
carica_training <- function(params) {
  cat("=== CARICAMENTO TRAINING SET ===\n")

  las_legno  <- readLAS(params$training_legno_path,  select="*")
  las_foglia <- readLAS(params$training_foglia_path, select="*")

  cat(sprintf("  Legno:  %d punti\n", npoints(las_legno)))
  cat(sprintf("  Foglia: %d punti\n", npoints(las_foglia)))

  # Sottocampiona se necessario (velocita + bilanciamento)
  set.seed(params$seed)
  n_max <- params$max_train_pts

  if (npoints(las_legno) > n_max) {
    cat(sprintf("  Sottocampionamento legno: %d -> %d punti\n",
                npoints(las_legno), n_max))
    idx <- sample(npoints(las_legno), n_max)
    las_legno <- las_legno[idx]
  }
  if (npoints(las_foglia) > n_max) {
    cat(sprintf("  Sottocampionamento foglia: %d -> %d punti\n",
                npoints(las_foglia), n_max))
    idx <- sample(npoints(las_foglia), n_max)
    las_foglia <- las_foglia[idx]
  }

  # Estrai feature
  feat_legno  <- estrai_feature(las_legno,  params$scale_geometria, params, " (legno)")
  feat_foglia <- estrai_feature(las_foglia, params$scale_geometria, params, " (foglia)")

  # Assembla dataset
  df_legno  <- feat_legno;  df_legno$classe  <- "legno"
  df_foglia <- feat_foglia; df_foglia$classe <- "foglia"
  df        <- bind_rows(df_legno, df_foglia)
  df$classe <- factor(df$classe, levels=c("foglia","legno"))

  # Rimuovi colonne con NA o varianza zero
  df <- df[, colSums(is.na(df)) == 0, drop=FALSE]
  df <- df[, apply(df[, -ncol(df), drop=FALSE], 2,
                   function(x) is.character(x) || var(as.numeric(x), na.rm=TRUE) > 1e-10),
           drop=FALSE]

  cat(sprintf("  Training set: %d punti x %d feature\n",
              nrow(df), ncol(df)-1))
  cat("  Feature:", paste(names(df)[-ncol(df)], collapse=", "), "\n\n")

  df
}

# =============================================================================
# 4. TRAINING RANDOM FOREST CON CROSS-VALIDAZIONE
# =============================================================================
addestra_RF <- function(df, params) {
  cat("=== TRAINING RANDOM FOREST ===\n")

  set.seed(params$seed)
  X <- df[, names(df) != "classe", drop=FALSE]
  y <- df$classe

  # Cross-validazione k-fold per stimare accuratezza
  cat(sprintf("  Cross-validazione %d-fold...\n", params$n_fold_cv))
  folds  <- createFolds(y, k=params$n_fold_cv, list=TRUE)
  acc_cv <- numeric(params$n_fold_cv)
  ba_cv  <- numeric(params$n_fold_cv)

  for (f in seq_along(folds)) {
    idx_val <- folds[[f]]
    idx_trn <- setdiff(seq_len(nrow(df)), idx_val)

    mod_f <- ranger(classe ~ ., data=df[idx_trn,],
                    num.trees=200,         # meno alberi per CV (velocita)
                    probability=TRUE,
                    seed=params$seed)

    prob_val <- predict(mod_f, data=df[idx_val,])$predictions[,"foglia"]
    pred_val <- factor(ifelse(prob_val >= params$soglia_foglia, "foglia", "legno"),
                       levels=c("foglia","legno"))
    cm       <- confusionMatrix(pred_val, y[idx_val], positive="foglia")

    acc_cv[f] <- cm$overall["Accuracy"]
    ba_cv[f]  <- cm$byClass["Balanced Accuracy"]
    cat(sprintf("    Fold %d: Accuracy=%.3f | Balanced Accuracy=%.3f\n",
                f, acc_cv[f], ba_cv[f]))
  }

  cat(sprintf("  CV media: Accuracy=%.3f±%.3f | Balanced Accuracy=%.3f±%.3f\n\n",
              mean(acc_cv), sd(acc_cv), mean(ba_cv), sd(ba_cv)))

  # Modello finale su tutto il training set
  cat(sprintf("  Training modello finale (%d alberi)...\n", params$n_trees))
  modello <- ranger(classe ~ ., data=df,
                    num.trees      = params$n_trees,
                    probability    = TRUE,
                    importance     = "permutation",
                    seed           = params$seed)

  cat("  Training completato.\n\n")

  # Importanza feature
  imp <- sort(modello$variable.importance, decreasing=TRUE)
  cat("  Feature importanza (top 10):\n")
  for (i in seq_len(min(10, length(imp))))
    cat(sprintf("    %-20s %.4f\n", names(imp)[i], imp[i]))
  cat("\n")

  list(modello  = modello,
       acc_cv   = acc_cv,
       ba_cv    = ba_cv,
       feature  = names(X),
       n_train  = nrow(df))
}

# =============================================================================
# 5. PLOT FEATURE IMPORTANCE
# =============================================================================
plot_importance <- function(rf_obj, outdir) {
  imp <- sort(rf_obj$modello$variable.importance, decreasing=TRUE)
  df  <- data.frame(feature=names(imp), importanza=imp,
                    stringsAsFactors=FALSE)
  df$feature <- factor(df$feature, levels=rev(df$feature))

  # Raggruppa per tipo (geometrica vs radiometrica)
  df$tipo <- case_when(
    grepl("^lin|^plan|^sfer|^curv|^dens", df$feature) ~ "Geometrica",
    grepl("^R$|^G$|^B$|^ExG", df$feature)            ~ "RGB",
    grepl("intensity", df$feature)                    ~ "Intensita",
    grepl("return", df$feature)                       ~ "Ritorni",
    grepl("Z_norm", df$feature)                       ~ "Altezza",
    TRUE                                               ~ "Altro"
  )

  p <- ggplot(df, aes(x=feature, y=importanza, fill=tipo)) +
    geom_col() +
    scale_fill_viridis(discrete=TRUE) +
    coord_flip() +
    labs(title="Importanza feature - Random Forest",
         subtitle="Permutation importance",
         x="", y="Importanza", fill="Tipo") +
    theme_minimal(base_size=11)
  ggsave(file.path(outdir, "00_RF_feature_importance.png"),
         p, width=10, height=max(5, nrow(df)*0.25+2), dpi=150)
  cat("  Plot feature importance salvato.\n")
}

# =============================================================================
# 6. CLASSIFICA SINGOLO LAS
# =============================================================================
classifica_las <- function(percorso, rf_obj, params, outdir) {
  tree_name <- file_path_sans_ext(basename(percorso))
  cat(sprintf("\n--- %s ---\n", tree_name))
  t0 <- proc.time()

  las <- readLAS(percorso, select="*")
  if (is.null(las) || npoints(las)==0) {
    cat("  ERRORE: file vuoto\n")
    return(NULL)
  }

  # Estrai feature (stesse del training)
  feat <- estrai_feature(las, params$scale_geometria, params)

  # Mantieni solo le feature usate nel training
  feat_comuni <- intersect(rf_obj$feature, names(feat))
  feat_mancanti <- setdiff(rf_obj$feature, names(feat))
  if (length(feat_mancanti) > 0)
    cat("  WARN: feature mancanti:", paste(feat_mancanti, collapse=", "), "\n")

  # Predizione
  cat("  Classificazione RF...\n")
  feat_pred <- feat[, feat_comuni, drop=FALSE]
  # Gestisci NA con imputazione mediana
  for (col in names(feat_pred))
    feat_pred[[col]][is.na(feat_pred[[col]])] <- median(feat_pred[[col]], na.rm=TRUE)

  pred      <- predict(rf_obj$modello, data=feat_pred)
  prob_leaf <- pred$predictions[, "foglia"]
  classe_rf <- ifelse(prob_leaf >= params$soglia_foglia,
                      params$class_leaf, params$class_wood)

  n_foglia <- sum(classe_rf == params$class_leaf)
  n_legno  <- sum(classe_rf == params$class_wood)
  cat(sprintf("  Foglia: %d (%.1f%%) | Legno: %d (%.1f%%)\n",
              n_foglia, 100*n_foglia/npoints(las),
              n_legno,  100*n_legno /npoints(las)))

  # Aggiunge campi RF al LAS
  las <- add_lasattribute(las, as.integer(classe_rf), "RF.class",
                          "RF classification: 1=foglia 2=legno")
  las <- add_lasattribute(las, as.numeric(round(prob_leaf, 4)), "RF.prob",
                          "RF P(foglia)")

  # Confronto con CANUPO se disponibile
  confronto <- NULL
  df_data   <- las@data
  if ("CANUPO.class" %in% names(df_data)) {
    canupo   <- as.integer(df_data[["CANUPO.class"]])
    # Mappa CANUPO ai nostri codici (1=foglia, 2=legno)
    concord  <- mean(classe_rf == canupo, na.rm=TRUE)
    cat(sprintf("  Concordanza RF vs CANUPO: %.1f%%\n", 100*concord))
    confronto <- list(concordanza_pct = round(100*concord, 1),
                      n_rf_fog_can_leg = sum(classe_rf==1L & canupo==2L),
                      n_rf_leg_can_fog = sum(classe_rf==2L & canupo==1L))
  }

  # Salva LAS classificato
  out_path <- file.path(outdir, paste0(tree_name, "_RF.las"))
  writeLAS(las, out_path)
  cat(sprintf("  Salvato: %s\n", basename(out_path)))

  elapsed <- (proc.time()-t0)["elapsed"]
  cat(sprintf("  Completato in %.1f s\n", elapsed))

  list(tree_name      = tree_name,
       n_pts          = npoints(las),
       n_foglia       = n_foglia,
       n_legno        = n_legno,
       pct_foglia     = round(100*n_foglia/npoints(las), 1),
       prob_media     = round(mean(prob_leaf[classe_rf==1L]), 3),
       concordanza_pct = if (!is.null(confronto)) confronto$concordanza_pct else NA,
       n_rf_fog_can_leg = if (!is.null(confronto)) confronto$n_rf_fog_can_leg else NA,
       n_rf_leg_can_fog = if (!is.null(confronto)) confronto$n_rf_leg_can_fog else NA,
       elapsed_sec    = round(elapsed, 1),
       status         = "OK")
}

# =============================================================================
# 7. BATCH
# =============================================================================
run_batch_RF <- function(params) {

  # Training
  df_train <- carica_training(params)
  rf_obj   <- addestra_RF(df_train, params)

  # Plot feature importance
  plot_importance(rf_obj, params$output_dir)

  # Salva modello per riuso futuro
  modello_path <- file.path(params$output_dir, "modello_RF.rds")
  saveRDS(rf_obj, modello_path)
  cat(sprintf("Modello salvato: %s\n\n", modello_path))

  # Classifica le 10 piante
  files <- list.files(params$input_dir, pattern="\\.las$",
                      full.names=TRUE, ignore.case=TRUE)
  if (length(files)==0) stop("Nessun .las in: ", params$input_dir)

  cat(sprintf("=== CLASSIFICAZIONE BATCH (%d file) ===\n", length(files)))

  risultati <- list()
  for (i in seq_along(files)) {
    cat(sprintf("\n[%d/%d] ", i, length(files)))
    res <- tryCatch(
      classifica_las(files[i], rf_obj, params, params$output_dir),
      error=function(e) {
        cat("  ERRORE:", conditionMessage(e), "\n")
        list(tree_name=file_path_sans_ext(basename(files[i])),
             status=paste("ERR:", conditionMessage(e)))
      }
    )
    if (!is.null(res)) risultati[[i]] <- res
    cat(rep("-", 60), "\n", sep="")
  }

  # Summary
  summary_df <- bind_rows(lapply(risultati, as.data.frame))
  write.csv2(summary_df,
             file.path(params$output_dir, "00_RF_summary.csv"),
             row.names=FALSE)

  # Plot comparativo RF vs CANUPO (se concordanza disponibile)
  ok <- summary_df[summary_df$status=="OK" & !is.na(summary_df$concordanza_pct),]
  if (nrow(ok) > 1) {

    # Concordanza per albero
    p_conc <- ggplot(ok, aes(x=reorder(tree_name, concordanza_pct),
                              y=concordanza_pct, fill=concordanza_pct)) +
      geom_col(show.legend=FALSE) +
      scale_fill_viridis(option="magma") + coord_flip() +
      labs(title="Concordanza RF vs CANUPO per albero",
           subtitle="100% = classificazione identica",
           x="Albero", y="Concordanza (%)") +
      theme_minimal(base_size=12)
    ggsave(file.path(params$output_dir, "00_RF_vs_CANUPO.png"),
           p_conc, width=9, height=6, dpi=150)

    # Punti discordanti per albero
    disc <- ok %>%
      select(tree_name, n_rf_fog_can_leg, n_rf_leg_can_fog) %>%
      pivot_longer(-tree_name, names_to="tipo", values_to="n") %>%
      mutate(tipo = ifelse(tipo=="n_rf_fog_can_leg",
                           "RF=foglia, CANUPO=legno",
                           "RF=legno, CANUPO=foglia"))

    p_disc <- ggplot(disc, aes(x=reorder(tree_name, n), y=n, fill=tipo)) +
      geom_col(position="dodge") +
      scale_fill_manual(values=c("RF=foglia, CANUPO=legno"  = "#2d8a4e",
                                  "RF=legno, CANUPO=foglia"  = "#8B4513")) +
      coord_flip() +
      labs(title="Punti classificati diversamente da RF e CANUPO",
           x="Albero", y="N punti discordanti", fill="") +
      theme_minimal(base_size=12) +
      theme(legend.position="bottom")
    ggsave(file.path(params$output_dir, "00_RF_discordanza.png"),
           p_disc, width=10, height=6, dpi=150)
  }

  # Riepilogo CV
  cat(sprintf("\n=== RIEPILOGO ===\n"))
  cat(sprintf("Cross-validazione %d-fold:\n", params$n_fold_cv))
  cat(sprintf("  Balanced Accuracy: %.3f ± %.3f\n",
              mean(rf_obj$ba_cv), sd(rf_obj$ba_cv)))
  cat(sprintf("  Accuracy:          %.3f ± %.3f\n",
              mean(rf_obj$acc_cv), sd(rf_obj$acc_cv)))
  cat(sprintf("Piante classificate: %d / %d\n",
              sum(summary_df$status=="OK", na.rm=TRUE), nrow(summary_df)))
  cat("Output:", params$output_dir, "\n")

  invisible(list(rf=rf_obj, summary=summary_df))
}

# =============================================================================
# 8. ESECUZIONE
# =============================================================================
if (!dir.exists(params$input_dir))
  stop("Cartella input non trovata: ", params$input_dir)
if (!file.exists(params$training_legno_path))
  stop("File training legno non trovato: ", params$training_legno_path)
if (!file.exists(params$training_foglia_path))
  stop("File training foglia non trovato: ", params$training_foglia_path)

risultati <- run_batch_RF(params)
