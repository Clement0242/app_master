library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)

# ===== CONSTANTES =====
GRAVITY <- 9.81
DEFAULT_TAKEOFF_LANDING_THRESHOLD <- 0.10
DEFAULT_MOVEMENT_WINDOW <- 10
DEFAULT_MOVEMENT_CONSECUTIVE <- 50
MIN_FLIGHT_TIME <- 0.1
MAX_FLIGHT_TIME <- 2.0
MIN_PUSH_TIME <- 0.05
MAX_PUSH_TIME <- 5.0

# ===== CLASSE RÉSULTAT =====
JumpAnalysisResult <- function(body_mass = NULL, body_weight = NULL,
                               movement_start_time = NA, takeoff_time = NA,
                               landing_time = NA, push_time = NA,
                               flight_time = NA, relative_threshold = NA,
                               mean_force_push = NA, mean_force_push_rel = NA,
                               mean_net_force_push = NA, mean_net_force_push_rel = NA,
                               max_force_push = NA, max_force_push_rel = NA,
                               max_net_force_push = NA, max_net_force_push_rel = NA,
                               net_impulse = NA, takeoff_velocity = NA,
                               jump_height = NA,
                               mean_velocity_push = NA,
                               max_velocity_push = NA,
                               velocity_data = NULL,
                               hpo = NA,
                               Pmax = NA,
                               Pmax_rel = NA,
                               Vdec_opt = NA,
                               SFVopt = NA,
                               is_valid = FALSE, validation_messages = character()) {
  structure(
    list(
      body_mass = body_mass,
      body_weight = body_weight,
      movement_start_time = movement_start_time,
      takeoff_time = takeoff_time,
      landing_time = landing_time,
      push_time = push_time,
      flight_time = flight_time,
      relative_threshold = relative_threshold,
      mean_force_push = mean_force_push,
      mean_force_push_rel = mean_force_push_rel,
      mean_net_force_push = mean_net_force_push,
      mean_net_force_push_rel = mean_net_force_push_rel,
      max_force_push = max_force_push,
      max_force_push_rel = max_force_push_rel,
      max_net_force_push = max_net_force_push,
      max_net_force_push_rel = max_net_force_push_rel,
      net_impulse = net_impulse,
      takeoff_velocity = takeoff_velocity,
      jump_height = jump_height,
      mean_velocity_push = mean_velocity_push,
      max_velocity_push = max_velocity_push,
      velocity_data = velocity_data,
      hpo = hpo,
      Pmax = Pmax,
      Pmax_rel = Pmax_rel,
      Vdec_opt = Vdec_opt,
      SFVopt = SFVopt,
      is_valid = is_valid,
      validation_messages = validation_messages
    ),
    class = "JumpAnalysisResult"
  )
}

# ===== FONCTIONS D'ANALYSE =====
import_force_data <- function(file_path) {
  if (!file.exists(file_path)) stop("Fichier non trouvé: ", file_path)
  
  tryCatch({
    raw_data <- read_csv(file_path, show_col_types = FALSE, col_types = cols(.default = col_double()))
    
    if (ncol(raw_data) < 2) stop("Au moins 2 colonnes requises")
    
    colnames(raw_data)[1:2] <- c("abs_time_s", "Fz_N")
    
    data <- raw_data %>%
      select(abs_time_s, Fz_N) %>%
      mutate(
        abs_time_s = as.numeric(abs_time_s),
        Fz_N = as.numeric(Fz_N)
      ) %>%
      filter(!is.na(abs_time_s), !is.na(Fz_N))
    
    if (nrow(data) == 0) stop("Aucune donnée valide")
    return(data)
  }, error = function(e) stop("Erreur import: ", e$message))
}

calculate_body_weight <- function(body_mass) {
  if (!is.numeric(body_mass) || body_mass <= 0) 
    stop("Masse invalide")
  return(body_mass * GRAVITY)
}

select_data_range <- function(data, start_time, end_time) {
  if (start_time >= end_time) stop("Plage invalide")
  
  selected_data <- data %>%
    filter(abs_time_s >= start_time & abs_time_s <= end_time)
  
  if (nrow(selected_data) == 0) stop("Aucune donnée dans la plage")
  return(selected_data)
}

calculate_relative_threshold <- function(data, n_samples = 250) {
  n_available <- min(n_samples, nrow(data))
  mean(data$Fz_N[1:n_available]) * 0.1
}

detect_movement_start <- function(data, body_weight, relative_threshold,
                                  window_size = 10, consecutive_points = 50) {
  force_deviation <- abs(data$Fz_N - body_weight)
  
  if (nrow(data) >= window_size) {
    smoothed_deviation <- rep(NA, length(force_deviation))
    half_window <- floor(window_size / 2)
    
    for (i in (half_window + 1):(length(force_deviation) - half_window)) {
      start_idx <- i - half_window
      end_idx <- i + half_window
      smoothed_deviation[i] <- mean(force_deviation[start_idx:end_idx])
    }
    deviation_to_use <- ifelse(is.na(smoothed_deviation), force_deviation, smoothed_deviation)
  } else {
    deviation_to_use <- force_deviation
  }
  
  above_threshold <- deviation_to_use > relative_threshold
  rle_result <- rle(above_threshold)
  cumsum_lengths <- cumsum(rle_result$lengths)
  
  for (i in seq_along(rle_result$values)) {
    if (rle_result$values[i] && rle_result$lengths[i] >= consecutive_points) {
      start_index <- if (i == 1) 1 else cumsum_lengths[i - 1] + 1
      return(data$abs_time_s[start_index])
    }
  }
  
  start_index <- which(deviation_to_use > relative_threshold)[1]
  if (is.na(start_index)) return(NA)
  return(data$abs_time_s[start_index])
}

detect_takeoff <- function(data, movement_start_time, body_weight, threshold_pct = 0.10) {
  if (is.na(movement_start_time)) return(NA)
  
  takeoff_threshold <- body_weight * threshold_pct
  post_movement_data <- data %>% filter(abs_time_s >= movement_start_time)
  
  if (nrow(post_movement_data) == 0) return(NA)
  
  takeoff_index <- which(post_movement_data$Fz_N < takeoff_threshold)[1]
  if (is.na(takeoff_index)) return(NA)
  
  return(post_movement_data$abs_time_s[takeoff_index])
}

detect_landing <- function(data, takeoff_time, body_weight, threshold_pct = 0.10) {
  if (is.na(takeoff_time)) return(NA)
  
  landing_threshold <- body_weight * threshold_pct
  post_takeoff_data <- data %>% filter(abs_time_s > takeoff_time)
  
  if (nrow(post_takeoff_data) == 0) return(NA)
  
  landing_index <- which(post_takeoff_data$Fz_N > landing_threshold)[1]
  if (is.na(landing_index)) return(NA)
  
  return(post_takeoff_data$abs_time_s[landing_index])
}

validate_results <- function(movement_start_time, takeoff_time, landing_time,
                             push_time, flight_time) {
  messages <- character()
  is_valid <- TRUE
  
  if (any(is.na(c(movement_start_time, takeoff_time, landing_time)))) {
    messages <- c(messages, "Événements non détectés")
    is_valid <- FALSE
  } else {
    if (movement_start_time >= takeoff_time) {
      messages <- c(messages, "Chronologie incorrecte")
      is_valid <- FALSE
    }
    if (takeoff_time >= landing_time) {
      messages <- c(messages, "Chronologie incorrecte")
      is_valid <- FALSE
    }
    
    if (!is.na(flight_time)) {
      if (flight_time < MIN_FLIGHT_TIME || flight_time > MAX_FLIGHT_TIME) {
        messages <- c(messages, sprintf("Temps de vol hors limites: %.0f ms", flight_time * 1000))
        is_valid <- FALSE
      }
    }
    
    if (!is.na(push_time)) {
      if (push_time < MIN_PUSH_TIME || push_time > MAX_PUSH_TIME) {
        messages <- c(messages, sprintf("Temps poussée hors limites: %.0f ms", push_time * 1000))
        is_valid <- FALSE
      }
    }
  }
  
  if (is_valid) messages <- c(messages, "✓ Valide")
  
  return(list(is_valid = is_valid, messages = messages))
}

# ===== CALCULS DU PROFIL OPTIMAL (Samozino et al. 2008, 2012) =====

# Fonction 1 : Calcul de la distance de poussée (hpo) par intégration de la vitesse
calculate_push_distance <- function(velocity_data) {
  if (is.null(velocity_data) || nrow(velocity_data) < 2) {
    return(NA)
  }
  
  tryCatch({
    # Intégration numérique : distance = ∫ v(t) dt
    # Méthode des trapèzes : somme des aires (v1 + v2)/2 × dt
    dt <- diff(velocity_data$abs_time_s)
    v_mean_intervals <- (velocity_data$velocity_ms[-1] + velocity_data$velocity_ms[-nrow(velocity_data)]) / 2
    
    # Distance totale de poussée
    hpo <- sum(v_mean_intervals * dt)
    
    # Vérification cohérence (doit être > 0 et < 1m pour un squat jump)
    if (hpo < 0 || hpo > 1.5) {
      warning(paste("Distance de poussée incohérente:", round(hpo, 3), "m"))
      return(NA)
    }
    
    return(hpo)
  }, error = function(e) {
    warning(paste("Erreur calcul hpo:", e$message))
    return(NA)
  })
}

# Fonction 2 : Calcul du profil force-vitesse optimal (Samozino et al. 2012)
calculate_optimal_FV_profile <- function(Pmax, body_mass, hpo) {
  # Vérification des inputs
  if (any(is.na(c(Pmax, body_mass, hpo))) || 
      Pmax <= 0 || body_mass <= 0 || hpo <= 0) {
    return(list(Vdec_opt = NA, SFVopt = NA))
  }
  
  tryCatch({
    # Formules de Samozino et Morin (2008, 2012)

    # Calcul de la vitesse de décollage optimale
    # Vdec_opt^2 = 2 * hpo * ((Pmax / (m * hpo)) - (g / 2))
    term_inside_sqrt <- 2 * hpo * ((Pmax / (body_mass * hpo)) - (GRAVITY / 2))

    if (term_inside_sqrt <= 0) {
      warning("Impossible de calculer Vdec_opt : argument négatif sous la racine")
      return(list(Vdec_opt = NA, SFVopt = NA))
    }

    Vdec_opt <- sqrt(term_inside_sqrt)

    # Calcul de la pente optimale force-vitesse
    # SFVopt = - (g * (hpo + Vdec_opt^2 / (2 * g))) / (2 * (Pmax / m))
    denominator <- 2 * (Pmax / body_mass)

    if (denominator == 0) {
      warning("Impossible de calculer SFVopt : dénominateur nul")
      return(list(Vdec_opt = NA, SFVopt = NA))
    }

    SFVopt <- - (GRAVITY * (hpo + (Vdec_opt^2 / (2 * GRAVITY)))) / denominator

    return(list(Vdec_opt = Vdec_opt, SFVopt = SFVopt))

  }, error = function(e) {
    warning(paste("Erreur calcul profil optimal:", e$message))
    return(list(Vdec_opt = NA, SFVopt = NA))
  })
}

# ===== FONCTION PRINCIPALE AVEC CALCUL DE VITESSE =====
analyze_jump_data_with_selection <- function(selected_data, body_mass,
                                             takeoff_landing_threshold = 0.10,
                                             movement_window_size = 10,
                                             movement_consecutive_points = 50) {
  tryCatch({
    body_weight <- calculate_body_weight(body_mass)
    relative_threshold <- calculate_relative_threshold(selected_data)
    
    movement_start_time <- detect_movement_start(selected_data, body_weight, 
                                                 relative_threshold,
                                                 movement_window_size, 
                                                 movement_consecutive_points)
    takeoff_time <- detect_takeoff(selected_data, movement_start_time, 
                                   body_weight, takeoff_landing_threshold)
    landing_time <- detect_landing(selected_data, takeoff_time, 
                                   body_weight, takeoff_landing_threshold)
    
    push_time <- if (!is.na(movement_start_time) && !is.na(takeoff_time)) 
      takeoff_time - movement_start_time else NA
    
    flight_time <- if (!is.na(takeoff_time) && !is.na(landing_time)) 
      landing_time - takeoff_time else NA
    
    # Initialisation
    mean_force_push <- mean_force_push_rel <- NA
    mean_net_force_push <- mean_net_force_push_rel <- NA
    max_force_push <- max_force_push_rel <- NA
    max_net_force_push <- max_net_force_push_rel <- NA
    net_impulse <- takeoff_velocity <- jump_height <- NA
    mean_velocity_push <- max_velocity_push <- NA
    velocity_data <- NULL
    
    # Initialisation des nouvelles métriques
    hpo <- NA
    Pmax <- NA
    Pmax_rel <- NA
    Vdec_opt <- NA
    SFVopt <- NA

    # ===== CALCULS DYNAMIQUES AVEC VITESSE =====
    if (!is.na(movement_start_time) && !is.na(takeoff_time)) {
      push_data <- selected_data %>%
        filter(abs_time_s >= movement_start_time & abs_time_s <= takeoff_time) %>%
        arrange(abs_time_s) %>%
        mutate(
          # 1. Force nette (sans le poids)
          force_nette_N = Fz_N - body_weight,
          
          # 2. Accélération instantanée (a = F_nette / m)
          acceleration_ms2 = force_nette_N / body_mass,
          
          # 3. Intervalle de temps entre mesures
          dt = c(0, diff(abs_time_s)),
          
          # 4. Vitesse instantanée (intégration de l'accélération)
          velocity_ms = cumsum(acceleration_ms2 * dt)
        )
      
      # Stocker les données de vitesse pour le graphique
      velocity_data <- push_data %>% select(abs_time_s, velocity_ms, acceleration_ms2)
      
      # Intervalle de temps moyen
      dt_mean <- mean(diff(push_data$abs_time_s))
      
      # Forces moyennes et max
      mean_force_push <- mean(push_data$Fz_N)
      mean_force_push_rel <- mean_force_push / body_mass
      
      mean_net_force_push <- mean(push_data$force_nette_N)
      mean_net_force_push_rel <- mean_net_force_push / body_mass
      
      max_force_push <- max(push_data$Fz_N)
      max_force_push_rel <- max_force_push / body_mass
      
      max_net_force_push <- max(push_data$force_nette_N)
      max_net_force_push_rel <- max_net_force_push / body_mass
      
      # Impulsion et vitesse de décollage
      net_impulse <- sum(push_data$force_nette_N) * dt_mean
      takeoff_velocity <- net_impulse / body_mass
      jump_height <- (takeoff_velocity^2) / (2 * GRAVITY)

      # ===== NOUVELLES MÉTRIQUES DE VITESSE =====
      mean_velocity_push <- mean(push_data$velocity_ms, na.rm = TRUE)
      max_velocity_push <- max(push_data$velocity_ms, na.rm = TRUE)

      # ===== MÉTRIQUES SUPPLÉMENTAIRES =====
      hpo <- calculate_push_distance(velocity_data)

      if (!is.na(mean_net_force_push) && !is.na(mean_velocity_push)) {
        Pmax <- mean_net_force_push * mean_velocity_push
        Pmax_rel <- Pmax / body_mass
      }

      optimal_profile <- calculate_optimal_FV_profile(Pmax, body_mass, hpo)
      Vdec_opt <- optimal_profile$Vdec_opt
      SFVopt <- optimal_profile$SFVopt
    }

    validation <- validate_results(movement_start_time, takeoff_time, landing_time,
                                   push_time, flight_time)

    result <- JumpAnalysisResult(
      body_mass = body_mass,
      body_weight = body_weight,
      movement_start_time = movement_start_time,
      takeoff_time = takeoff_time,
      landing_time = landing_time,
      push_time = push_time,
      flight_time = flight_time,
      relative_threshold = relative_threshold,
      mean_force_push = mean_force_push,
      mean_force_push_rel = mean_force_push_rel,
      mean_net_force_push = mean_net_force_push,
      mean_net_force_push_rel = mean_net_force_push_rel,
      max_force_push = max_force_push,
      max_force_push_rel = max_force_push_rel,
      max_net_force_push = max_net_force_push,
      max_net_force_push_rel = max_net_force_push_rel,
      net_impulse = net_impulse,
      takeoff_velocity = takeoff_velocity,
      jump_height = jump_height,
      mean_velocity_push = mean_velocity_push,
      max_velocity_push = max_velocity_push,
      velocity_data = velocity_data,
      hpo = hpo,
      Pmax = Pmax,
      Pmax_rel = Pmax_rel,
      Vdec_opt = Vdec_opt,
      SFVopt = SFVopt,
      is_valid = validation$is_valid,
      validation_messages = validation$messages
    )

    return(result)
    
  }, error = function(e) {
    return(JumpAnalysisResult(
      body_mass = body_mass,
      is_valid = FALSE,
      validation_messages = paste("Erreur:", e$message)
    ))
  })
}

# ---- FIGURE COMPOSITE A-E : rapport graphique complet ----
make_jump_report <- function(selected_data, results,
                             athlete_name = "Athlète",
                             session_date = Sys.Date(),
                             sampling_hz = NA_real_) {
  stopifnot(!is.null(selected_data), !is.null(results))

  # Récup infos/événements
  bw   <- results$body_weight
  m    <- results$body_mass
  t0   <- results$movement_start_time
  t_to <- results$takeoff_time
  t_ld <- results$landing_time
  push_time <- results$push_time
  vdec <- results$takeoff_velocity
  h    <- results$jump_height
  vel_df <- results$velocity_data

  # 1) Reconstruire les données "push" et puissances instantanées
  # On ne garde que la poussée pour B/C/E et calculs avancés
  push_df <- selected_data |>
    dplyr::filter(abs_time_s >= t0, abs_time_s <= t_to) |>
    dplyr::arrange(abs_time_s) |>
    dplyr::mutate(
      force_nette = Fz_N - bw,
      dt = c(0, diff(abs_time_s))
    )

  # Joindre la vitesse/accélération calculées dans analyze_jump_data_with_selection()
  if (!is.null(vel_df) && nrow(vel_df) > 1) {
    push_df <- dplyr::left_join(push_df, vel_df, by = "abs_time_s")
  } else {
    # Si pas de vitesse (ne devrait pas arriver vu ton code), on met des NA
    push_df$velocity_ms <- NA_real_
    push_df$acceleration_ms2 <- NA_real_
  }

  # Puissance instantanée (W) = force_nette (N) * vitesse (m/s)
  push_df <- push_df |>
    dplyr::mutate(power_w = force_nette * velocity_ms)

  # Pmax instantanée (vraie) et relative
  Pmax_inst <- suppressWarnings(max(push_df$power_w, na.rm = TRUE))
  if (!is.finite(Pmax_inst)) Pmax_inst <- NA_real_
  Pmax_rel  <- if (is.finite(Pmax_inst) && !is.na(m) && m > 0) Pmax_inst / m else NA_real_

  # 2) RFD (méthode 1/3 – 2/3 de la poussée)
  # t1 = t0 + 1/3 * push_time ; t2 = t0 + 2/3 * push_time
  t1 <- if (!is.na(push_time)) t0 + push_time/3 else NA_real_
  t2 <- if (!is.na(push_time)) t0 + 2*push_time/3 else NA_real_

  # Interpolation force nette aux instants t1 et t2
  # On protège l’interpolation si t1/t2 en dehors du domaine
  interp_force <- function(tq) {
    if (is.na(tq)) return(NA_real_)
    if (tq < min(push_df$abs_time_s, na.rm = TRUE) ||
        tq > max(push_df$abs_time_s, na.rm = TRUE)) return(NA_real_)
    approx(x = push_df$abs_time_s, y = push_df$force_nette, xout = tq)$y
  }
  F1 <- interp_force(t1)
  F2 <- interp_force(t2)

  RFD_13_23 <- if (is.finite(F1) && is.finite(F2) && is.finite(t1) && is.finite(t2) && (t2 > t1)) {
    (F2 - F1) / (t2 - t1)  # N/s
  } else NA_real_

  # Temps au pic de force (en % de la poussée)
  idx_Fmax <- suppressWarnings(which.max(push_df$force_nette))
  t_Fmax   <- if (length(idx_Fmax) && is.finite(push_df$abs_time_s[idx_Fmax])) push_df$abs_time_s[idx_Fmax] else NA_real_
  time_to_Fmax_pct <- if (is.finite(t_Fmax) && !is.na(push_time) && push_time > 0) {
    100 * (t_Fmax - t0) / push_time
  } else NA_real_

  # 3) PANELS

  # ----- Panel A : Force-Temps complet (avec zones/événements) -----
  # On utilise la fenêtre sélectionnée par l’utilisateur (inclut classiquement poussée + vol + réception)
  full_df <- selected_data |> dplyr::arrange(abs_time_s)

  # Définition des zones (statique / poussée / vol / réception) pour shading
  # NB : si un des temps est NA, on omet la zone correspondante.
  zones <- list(
    static    = if (is.finite(t0)) c(xmin = min(full_df$abs_time_s), xmax = t0) else NULL,
    push      = if (is.finite(t0) && is.finite(t_to)) c(xmin = t0, xmax = t_to) else NULL,
    flight    = if (is.finite(t_to) && is.finite(t_ld)) c(xmin = t_to, xmax = t_ld) else NULL,
    reception = if (is.finite(t_ld)) c(xmin = t_ld, xmax = max(full_df$abs_time_s)) else NULL
  )

  pA <- ggplot(full_df, aes(abs_time_s, Fz_N)) +
    { if (!is.null(zones$static))    annotate("rect", xmin = zones$static["xmin"], xmax = zones$static["xmax"], ymin = -Inf, ymax = Inf, fill = "grey90") } +
    { if (!is.null(zones$push))      annotate("rect", xmin = zones$push["xmin"], xmax = zones$push["xmax"], ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.08) } +
    { if (!is.null(zones$flight))    annotate("rect", xmin = zones$flight["xmin"], xmax = zones$flight["xmax"], ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5) } +
    { if (!is.null(zones$reception)) annotate("rect", xmin = zones$reception["xmin"], xmax = zones$reception["xmax"], ymin = -Inf, ymax = Inf, fill = "grey90") } +
    geom_hline(yintercept = bw, linetype = "dashed", size = 0.7, color = "blue") +
    geom_line(size = 0.7, color = "blue") +
    { if (is.finite(t0)) geom_vline(xintercept = t0, color = "darkgreen", size = 0.9) } +
    { if (is.finite(t_to)) geom_vline(xintercept = t_to, color = "orange", size = 0.9) } +
    { if (is.finite(t_ld)) geom_vline(xintercept = t_ld, color = "purple", size = 0.9) } +
    labs(title = "A — Force verticale (N) vs Temps (s)",
         subtitle = sprintf("%s — %s — Masse: %.1f kg", athlete_name, as.character(session_date), m),
         x = "Temps (s)", y = "Force (N)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # Annotations clés Panel A
  ann_A <- c(
    if (is.finite(push_time)) sprintf("Temps poussée: %d ms", round(push_time*1000)) else NULL,
    if (is.finite(floor(results$flight_time*1000))) sprintf("Temps vol: %d ms", round(results$flight_time*1000)) else NULL,
    if (is.finite(h)) sprintf("Hauteur: %.1f cm", h*100) else NULL
  )
  if (length(ann_A)) {
    pA <- pA + annotate("text", x = Inf, y = Inf, label = paste(ann_A, collapse = "  •  "),
                        hjust = 1.02, vjust = 1.5, size = 3.5)
  }

  # ----- Panel B : Vitesse (poussée) -----
  pB <- ggplot(push_df, aes(abs_time_s, velocity_ms)) +
    geom_line(size = 1.0, color = "orange") +
    { if (is.finite(results$mean_velocity_push)) geom_hline(yintercept = results$mean_velocity_push, linetype = "dashed", color = "orange") } +
    { if (is.finite(vdec)) geom_hline(yintercept = vdec, linetype = "dotted", color = "red") } +
    labs(title = "B — Vitesse (m/s) pendant la poussée",
         x = "Temps (s)", y = "Vitesse (m/s)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # ----- Panel C : Accélération (poussée) -----
  pC <- ggplot(push_df, aes(abs_time_s, acceleration_ms2)) +
    geom_line(size = 1.0, color = "darkgreen") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "C — Accélération (m/s²) pendant la poussée",
         x = "Temps (s)", y = "Accélération (m/s²)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # ----- Panel D : Métriques principales (tableau texte) -----
  # On construit un "dashboard" textuel
  met_lines <- c(
    sprintf("Hauteur saut: %s", if (is.finite(h)) sprintf("%.1f cm", 100*h) else "N/A"),
    sprintf("Temps vol: %s", if (is.finite(results$flight_time)) sprintf("%d ms", round(results$flight_time*1000)) else "N/A"),
    sprintf("Temps poussée: %s", if (is.finite(push_time)) sprintf("%d ms", round(push_time*1000)) else "N/A"),
    sprintf("Force max: %s", if (isTRUE(is.finite(max(selected_data$Fz_N)))) sprintf("%.0f N (%.1f N/kg)", max(selected_data$Fz_N, na.rm = TRUE), max(selected_data$Fz_N, na.rm = TRUE)/m) else "N/A"),
    sprintf("Force moy (poussée): %s", if (is.finite(results$mean_force_push)) sprintf("%.0f N (%.1f N/kg)", results$mean_force_push, results$mean_force_push/m) else "N/A"),
    sprintf("Impulsion nette: %s", if (is.finite(results$net_impulse)) sprintf("%.1f N·s", results$net_impulse) else "N/A"),
    sprintf("Vitesse décollage: %s", if (is.finite(vdec)) sprintf("%.2f m/s", vdec) else "N/A"),
    sprintf("Vitesse moy (poussée): %s", if (is.finite(results$mean_velocity_push)) sprintf("%.2f m/s", results$mean_velocity_push) else "N/A"),
    sprintf("Puissance max (inst.): %s", if (is.finite(Pmax_inst)) sprintf("%.0f W (%.1f W/kg)", Pmax_inst, Pmax_rel) else "N/A"),
    sprintf("RFD (1/3–2/3): %s", if (is.finite(RFD_13_23)) sprintf("%.0f N/s", RFD_13_23) else "N/A"),
    sprintf("Ratio Fmax/Poids: %s", if (isTRUE(is.finite(max(selected_data$Fz_N)))) sprintf("%.2f", max(selected_data$Fz_N, na.rm = TRUE)/bw) else "N/A"),
    sprintf("Temps au pic de force: %s", if (is.finite(time_to_Fmax_pct)) sprintf("%.0f %% poussée", time_to_Fmax_pct) else "N/A")
  )
  d_text <- data.frame(x = 0, y = rev(seq_along(met_lines)), lab = rev(met_lines))
  pD <- ggplot(d_text, aes(x, y, label = lab)) +
    geom_text(hjust = 0, size = 3.6) +
    xlim(0, 1) + ylim(0, length(met_lines)+1) +
    labs(title = "D — Métriques principales") +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # ----- Panel E : Profil technique (Force nette vs Vitesse) -----
  pE <- ggplot(push_df, aes(x = velocity_ms, y = force_nette)) +
    geom_point(size = 1.6, alpha = 0.6, color = "blue") +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
    labs(title = "E — Force nette (N) vs Vitesse (m/s)",
         x = "Vitesse (m/s)", y = "Force nette (N)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  # 4) Assemblage (portrait) : A pleine largeur ; B|C ; D|E
  top <- pA
  mid <- pB + pC
  bot <- pD + pE
  fig <- top / mid / bot +
    plot_annotation(
      caption = sprintf("Méthodo : seuil décollage 10%% du poids | %s | %s",
                        if (is.finite(sampling_hz)) sprintf("Fréq. échantillonnage %.0f Hz", sampling_hz) else "Fréq. échantillonnage N/A",
                        "Modèle force-vitesse : intégration vitesse, Pmax instantanée"),
      theme = theme(plot.caption = element_text(size = 9))
    )

  # On renvoie l’objet ggplot patchwork + une liste de métriques utiles
  list(
    figure = fig,
    metrics = list(
      Pmax_inst = Pmax_inst,
      Pmax_rel = Pmax_rel,
      RFD_13_23 = RFD_13_23,
      time_to_Fmax_pct = time_to_Fmax_pct
    )
  )
}

# ---- Helper pratique pour sauvegarder en PDF ----
save_jump_report_pdf <- function(report_obj, file = "rapport_saut.pdf", width = 8.27, height = 11.69) {
  # dimensions par défaut = A4 portrait en pouces (approx 8.27 x 11.69)
  stopifnot(!is.null(report_obj$figure))
  ggplot2::ggsave(filename = file, plot = report_obj$figure, width = width, height = height, units = "in")
  invisible(file)
}

# ===== UI =====
ui <- dashboardPage(
  dashboardHeader(title = "Analyse Saut + Vitesse"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Analyse", tabName = "analysis", icon = icon("chart-line")),
      menuItem("Vitesse", tabName = "velocity", icon = icon("tachometer-alt")),
      menuItemOutput("report_menu"),
      menuItem("Paramètres", tabName = "settings", icon = icon("cogs"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "analysis",
              fluidRow(
                box(
                  title = "1. Chargement", status = "primary", solidHeader = TRUE, width = 4,
                  fileInput("file", "Fichier .txt", accept = c(".txt")),
                  numericInput("body_mass", "Masse (kg)", 70, min = 30, max = 150, step = 0.1),
                  textInput("athlete_name", "Nom de l’athlète", value = ""),
                  hr(),

                  h4("2. Sélection de plage"),
                  conditionalPanel(
                    condition = "output.file_loaded",
                    fluidRow(
                      column(6, numericInput("start_time", "Début (s)", 0, step = 0.1)),
                      column(6, numericInput("end_time", "Fin (s)", 5, step = 0.1))
                    ),
                    actionButton("select_range", "Valider", class = "btn-info", width = "100%"),
                    br(), br(),
                    conditionalPanel(
                      condition = "output.range_selected",
                      actionButton("analyze", "Analyser", class = "btn-success", width = "100%")
                    )
                  )
                ),
                
                box(
                  title = "Force verticale", status = "success", solidHeader = TRUE, width = 8,
                  plotOutput("main_plot", height = "500px")
                )
              ),
              
              fluidRow(
                conditionalPanel(
                  condition = "output.analysis_done",
                  valueBoxOutput("push_time_box"),
                  valueBoxOutput("mean_velocity_box"),
                  valueBoxOutput("jump_height_box")
                )
              ),
              
              fluidRow(
                conditionalPanel(
                  condition = "output.analysis_done",
                  box(
                    title = "Résultats", status = "info", solidHeader = TRUE, width = 12,
                    tableOutput("results_table")
                  )
                )
              )
      ),
      
      # === NOUVEL ONGLET VITESSE ===
      tabItem(tabName = "velocity",
              fluidRow(
                conditionalPanel(
                  condition = "output.analysis_done",
                  box(
                    title = "Vitesse instantanée pendant la poussée",
                    status = "warning", solidHeader = TRUE, width = 12,
                    plotOutput("velocity_plot", height = "400px")
                  )
                )
              ),

              fluidRow(
                conditionalPanel(
                  condition = "output.analysis_done",
                  box(
                    title = "Accélération instantanée pendant la poussée",
                    status = "info", solidHeader = TRUE, width = 12,
                    plotOutput("acceleration_plot", height = "400px")
                  )
                )
              )
      ),

      tabItem(tabName = "report",
              fluidRow(
                conditionalPanel(
                  condition = "output.analysis_done",
                  box(
                    title = "Rapport complet du saut",
                    status = "primary", solidHeader = TRUE, width = 12,
                    plotOutput("report_plot", height = "900px"),
                    br(),
                    downloadButton("download_report", "Télécharger en PDF", icon = icon("file-download"))
                  )
                )
              )
      ),

      tabItem(tabName = "settings",
              box(
                title = "Paramètres avancés", status = "warning", solidHeader = TRUE, width = 6,
                numericInput("movement_window", "Fenêtre lissage", 10, min = 3, max = 20),
                numericInput("consecutive_points", "Points consécutifs", 50, min = 5, max = 100),
                numericInput("takeoff_threshold", "Seuil décollage (%)", 0.10, min = 0.05, max = 0.20, step = 0.01),
                actionButton("reset_params", "Réinitialiser", class = "btn-warning")
              )
      )
    )
  )
)

# ===== SERVER =====
server <- function(input, output, session) {
  
  original_data <- reactiveVal(NULL)
  selected_data <- reactiveVal(NULL)
  analysis_results <- reactiveVal(NULL)
  report_obj <- reactiveVal(NULL)
  
  observeEvent(input$reset_params, {
    updateNumericInput(session, "movement_window", value = 10)
    updateNumericInput(session, "consecutive_points", value = 50)
    updateNumericInput(session, "takeoff_threshold", value = 0.10)
  })
  
  observeEvent(input$file, {
    req(input$file)
    selected_data(NULL)
    analysis_results(NULL)
    report_obj(NULL)
    tryCatch({
      data <- import_force_data(input$file$datapath)
      original_data(data)

      updateNumericInput(session, "start_time",
                         value = round(min(data$abs_time_s), 2),
                         min = min(data$abs_time_s),
                         max = max(data$abs_time_s))
      
      updateNumericInput(session, "end_time", 
                         value = round(max(data$abs_time_s), 2),
                         min = min(data$abs_time_s),
                         max = max(data$abs_time_s))
      
      showNotification("Fichier chargé!", type = "success", duration = 3)
    }, error = function(e) {
      showNotification(paste("Erreur:", e$message), type = "error")
    })
  })
  
  observeEvent(input$select_range, {
    req(original_data(), input$start_time, input$end_time)
    tryCatch({
      data <- original_data()
      selected <- select_data_range(data, input$start_time, input$end_time)
      selected_data(selected)
      analysis_results(NULL)
      report_obj(NULL)
      showNotification(paste(nrow(selected), "échantillons sélectionnés"),
                       type = "success", duration = 3)
    }, error = function(e) {
      showNotification(paste("Erreur:", e$message), type = "error")
    })
  })
  
  observeEvent(input$analyze, {
    req(selected_data())
    showNotification("Analyse...", type = "message", duration = NULL, id = "analyzing")
    
    tryCatch({
      results <- analyze_jump_data_with_selection(
        selected_data = selected_data(),
        body_mass = input$body_mass,
        takeoff_landing_threshold = input$takeoff_threshold,
        movement_window_size = input$movement_window,
        movement_consecutive_points = input$consecutive_points
      )

      analysis_results(results)
      athlete_name <- input$athlete_name
      if (is.null(athlete_name) || !nzchar(trimws(athlete_name))) {
        athlete_name <- "Athlète"
      }

      selected <- selected_data()
      sampling_hz <- NA_real_
      if (!is.null(selected) && nrow(selected) > 1) {
        dt <- diff(selected$abs_time_s)
        dt <- dt[is.finite(dt) & dt > 0]
        if (length(dt) > 0) {
          sampling_hz <- 1 / mean(dt)
        }
      }

      report <- tryCatch({
        make_jump_report(
          selected_data = selected,
          results = results,
          athlete_name = athlete_name,
          session_date = Sys.Date(),
          sampling_hz = sampling_hz
        )
      }, error = function(err) {
        showNotification(paste("Erreur rapport:", err$message), type = "error")
        NULL
      })
      report_obj(report)
      removeNotification("analyzing")
      showNotification("Analyse terminée!", type = "success", duration = 3)
    }, error = function(e) {
      removeNotification("analyzing")
      report_obj(NULL)
      showNotification(paste("Erreur:", e$message), type = "error")
    })
  })
  
  output$file_loaded <- reactive(!is.null(original_data()))
  outputOptions(output, "file_loaded", suspendWhenHidden = FALSE)
  
  output$range_selected <- reactive(!is.null(selected_data()))
  outputOptions(output, "range_selected", suspendWhenHidden = FALSE)
  
  output$analysis_done <- reactive(!is.null(analysis_results()))
  outputOptions(output, "analysis_done", suspendWhenHidden = FALSE)

  output$report_menu <- renderMenu({
    req(report_obj())
    menuItem("Rapport", tabName = "report", icon = icon("file-pdf"))
  })

  output$report_plot <- renderPlot({
    report <- report_obj()
    req(report)
    report$figure
  })

  output$download_report <- downloadHandler(
    filename = function() {
      name <- input$athlete_name
      if (is.null(name) || !nzchar(trimws(name))) {
        name <- "rapport_saut"
      } else {
        name <- gsub("[^[:alnum:]]+", "_", tolower(trimws(name)))
        name <- gsub("_+", "_", name)
        name <- trimws(name, which = "both", whitespace = "_")
        if (!nzchar(name)) name <- "rapport_saut"
      }
      paste0(name, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      report <- report_obj()
      if (is.null(report)) {
        stop("Rapport indisponible")
      }
      save_jump_report_pdf(report, file = file)
    }
  )
  
  # === GRAPHIQUE PRINCIPAL ===
  output$main_plot <- renderPlot({
    data <- original_data()
    selected <- selected_data()
    results <- analysis_results()
    
    if (is.null(data)) {
      return(ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "Uploadez un fichier"), 
                                  size = 6, color = "gray") +
               xlim(0, 1) + ylim(0, 1) + theme_void())
    }
    
    body_weight <- input$body_mass * GRAVITY
    
    p <- ggplot(data, aes(x = abs_time_s, y = Fz_N)) +
      geom_line(color = if(is.null(selected)) "black" else "lightgray", size = 0.5) +
      geom_hline(yintercept = body_weight, color = "blue", linetype = "dashed", size = 1) +
      labs(title = "Force verticale", x = "Temps (s)", y = "Force (N)") +
      theme_minimal() + theme(plot.title = element_text(face = "bold"))
    
    if (!is.null(selected)) {
      p <- p +
        annotate("rect", 
                 xmin = min(selected$abs_time_s), 
                 xmax = max(selected$abs_time_s),
                 ymin = -Inf, ymax = Inf,
                 fill = "lightblue", alpha = 0.2) +
        geom_line(data = selected, aes(x = abs_time_s, y = Fz_N), 
                  color = "black", size = 1)
      
      if (!is.null(results)) {
        if (!is.na(results$movement_start_time)) 
          p <- p + geom_vline(xintercept = results$movement_start_time, color = "green", size = 1)
        if (!is.na(results$takeoff_time)) 
          p <- p + geom_vline(xintercept = results$takeoff_time, color = "orange", size = 1)
        if (!is.na(results$landing_time)) 
          p <- p + geom_vline(xintercept = results$landing_time, color = "purple", size = 1)
      }
    }
    
    return(p)
  })
  
  # === GRAPHIQUE VITESSE ===
  output$velocity_plot <- renderPlot({
    results <- analysis_results()
    if (is.null(results) || is.null(results$velocity_data)) return(NULL)
    
    ggplot(results$velocity_data, aes(x = abs_time_s, y = velocity_ms)) +
      geom_line(color = "#E69F00", size = 1.5) +
      geom_hline(yintercept = results$mean_velocity_push, 
                 color = "blue", linetype = "dashed", size = 1) +
      geom_hline(yintercept = results$takeoff_velocity, 
                 color = "red", linetype = "dotted", size = 1) +
      annotate("text", x = Inf, y = results$mean_velocity_push, 
               label = paste("Moy =", round(results$mean_velocity_push, 2), "m/s"), 
               hjust = 1.1, vjust = -0.5, color = "blue") +
      annotate("text", x = Inf, y = results$takeoff_velocity, 
               label = paste("Décollage =", round(results$takeoff_velocity, 2), "m/s"), 
               hjust = 1.1, vjust = -0.5, color = "red") +
      labs(
        title = "Vitesse instantanée calculée depuis la force",
        x = "Temps (s)", 
        y = "Vitesse (m/s)"
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
  })
  
  # === GRAPHIQUE ACCÉLÉRATION ===
  output$acceleration_plot <- renderPlot({
    results <- analysis_results()
    if (is.null(results) || is.null(results$velocity_data)) return(NULL)
    
    ggplot(results$velocity_data, aes(x = abs_time_s, y = acceleration_ms2)) +
      geom_line(color = "#009E73", size = 1.5) +
      geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      labs(
        title = "Accélération instantanée (a = F_nette / masse)",
        x = "Temps (s)", 
        y = "Accélération (m/s²)"
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"))
  })
  
  # === VALUE BOXES ===
  output$push_time_box <- renderValueBox({
    results <- analysis_results()
    if (is.null(results) || is.na(results$push_time)) {
      valueBox("N/A", "Temps poussée", icon = icon("arrow-up"), color = "gray")
    } else {
      valueBox(paste(round(results$push_time * 1000), "ms"), 
               "Temps poussée", icon = icon("arrow-up"), color = "blue")
    }
  })
  
  output$mean_velocity_box <- renderValueBox({
    results <- analysis_results()
    if (is.null(results) || is.na(results$mean_velocity_push)) {
      valueBox("N/A", "Vitesse moy poussée", icon = icon("tachometer-alt"), color = "gray")
    } else {
      valueBox(paste(round(results$mean_velocity_push, 2), "m/s"), 
               "Vitesse moy poussée", icon = icon("tachometer-alt"), color = "orange")
    }
  })
  
  output$jump_height_box <- renderValueBox({
    results <- analysis_results()
    if (is.null(results) || is.na(results$jump_height)) {
      valueBox("N/A", "Hauteur saut", icon = icon("chart-bar"), color = "gray")
    } else {
      valueBox(paste(round(results$jump_height * 100), "cm"), 
               "Hauteur saut", icon = icon("chart-bar"), color = "green")
    }
  })
  
  # === TABLEAU RÉSULTATS ===
  output$results_table <- renderTable({
    results <- analysis_results()
    if (is.null(results)) return(NULL)
    
    data.frame(
      Métrique = c(
        "Masse", "Poids", "Seuil",
        "Début mouvement", "Décollage", "Atterrissage",
        "Temps poussée", "Temps vol",
        "Force moy poussée", "Force nette moy",
        "Impulsion nette",
        "Vitesse MOY poussée", "Vitesse MAX poussée", "Vitesse décollage",
        "Hauteur saut",
        "Distance poussée (hpo)", "Puissance maximale (Pmax)", "Pente optimale (SFVopt)",
        "Validité"
      ),
      Valeur = c(
        paste(round(results$body_mass, 1), "kg"),
        paste(round(results$body_weight, 1), "N"),
        paste(round(results$relative_threshold, 2), "N"),
        paste(round(results$movement_start_time, 3), "s"),
        paste(round(results$takeoff_time, 3), "s"),
        paste(round(results$landing_time, 3), "s"),
        if(is.na(results$push_time)) "N/A" else paste(round(results$push_time * 1000), "ms"),
        if(is.na(results$flight_time)) "N/A" else paste(round(results$flight_time * 1000), "ms"),
        if(is.na(results$mean_force_push)) "N/A" else paste(round(results$mean_force_push, 1), "N"),
        if(is.na(results$mean_net_force_push)) "N/A" else paste(round(results$mean_net_force_push, 1), "N"),
        if(is.na(results$net_impulse)) "N/A" else paste(round(results$net_impulse, 2), "N·s"),
        if(is.na(results$mean_velocity_push)) "N/A" else paste(round(results$mean_velocity_push, 2), "m/s"),
        if(is.na(results$max_velocity_push)) "N/A" else paste(round(results$max_velocity_push, 2), "m/s"),
        if(is.na(results$takeoff_velocity)) "N/A" else paste(round(results$takeoff_velocity, 2), "m/s"),
        if(is.na(results$jump_height)) "N/A" else paste(round(results$jump_height * 100, 1), "cm"),
        if(is.na(results$hpo)) "N/A" else paste(round(results$hpo * 100, 1), "cm"),
        if(is.na(results$Pmax)) "N/A" else paste(round(results$Pmax, 1), "W"),
        if(is.na(results$SFVopt)) "N/A" else paste(round(results$SFVopt, 1), "N/(m/s)"),
        ifelse(results$is_valid, "✓ Valide", "✗ Non valide")
      )
    )
  }, bordered = TRUE, striped = TRUE, spacing = 'l')
}

shinyApp(ui = ui, server = server)