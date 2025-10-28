library(ggplot2)
library(dplyr)
library(readr)

# Constantes globales
GRAVITY <- 9.81
LINES_TO_SKIP <- 13
DEFAULT_BASELINE_SAMPLES <- 500
DEFAULT_TAKEOFF_LANDING_THRESHOLD <- 0.10
DEFAULT_MOVEMENT_WINDOW <- 10
DEFAULT_MOVEMENT_CONSECUTIVE <- 50
MIN_FLIGHT_TIME <- 0.1
MAX_FLIGHT_TIME <- 2.0
MIN_PUSH_TIME <- 0.05
MAX_PUSH_TIME <- 5.0

# Classe pour les résultats
JumpAnalysisResult <- function(body_mass = NULL, body_weight = NULL,
                               movement_start_time = NA, takeoff_time = NA,
                               landing_time = NA, push_time = NA,
                               flight_time = NA, relative_threshold = NA,
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
      is_valid = is_valid,
      validation_messages = validation_messages
    ),
    class = "JumpAnalysisResult"
  )
}

print.JumpAnalysisResult <- function(x, ...) {
  cat("=== Résultats d'Analyse de Saut ===\n")
  cat(sprintf("Masse corporelle: %.1f kg\n", x$body_mass))
  cat(sprintf("Poids corporel: %.1f N\n", x$body_weight))
  cat(sprintf("Seuil relatif: %.2f N\n", x$relative_threshold))
  cat("\n--- Événements temporels ---\n")
  cat(sprintf("Début du mouvement: %.3f s\n", x$movement_start_time))
  cat(sprintf("Décollage: %.3f s\n", x$takeoff_time))
  cat(sprintf("Atterrissage: %.3f s\n", x$landing_time))
  cat("\n--- Métriques ---\n")
  cat(sprintf("Temps de poussée: %.0f ms\n", x$push_time * 1000))
  cat(sprintf("Temps de vol: %.0f ms\n", x$flight_time * 1000))
  cat(sprintf("\nValidité: %s\n", ifelse(x$is_valid, "✓ Valide", "✗ Non valide")))
  if (length(x$validation_messages) > 0) {
    cat("Messages de validation:\n")
    for (msg in x$validation_messages) {
      cat(sprintf("  - %s\n", msg))
    }
  }
}

# Import des données
import_force_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Le fichier spécifié n'existe pas: ", file_path)
  }
  
  tryCatch({
    raw_data <- read_csv(file_path)
    
    if (ncol(raw_data) < 2) {
      stop("Le fichier doit contenir au moins 2 colonnes (temps et force)")
    }
    
    colnames(raw_data)[1:2] <- c("abs_time_s", "Fz_N")
    
    data <- raw_data %>%
      select(abs_time_s, Fz_N) %>%
      mutate(
        abs_time_s = as.numeric(abs_time_s),
        Fz_N = as.numeric(Fz_N)
      ) %>%
      filter(!is.na(abs_time_s), !is.na(Fz_N))
    
    if (nrow(data) == 0) {
      stop("Aucune donnée valide trouvée après nettoyage")
    }
    
    return(data)
    
  }, error = function(e) {
    stop("Erreur lors de l'importation: ", e$message)
  })
}

# Calcul du poids corporel
calculate_body_weight <- function(body_mass) {
  if (!is.numeric(body_mass) || body_mass <= 0) {
    stop("La masse corporelle doit être un nombre positif")
  }
  return(body_mass * GRAVITY)
}

# Sélection de plage de données
select_data_range <- function(data, start_time, end_time) {
  if (start_time >= end_time) {
    stop("Le temps de début doit être inférieur au temps de fin")
  }
  
  selected_data <- data %>%
    filter(abs_time_s >= start_time & abs_time_s <= end_time)
  
  if (nrow(selected_data) == 0) {
    stop("Aucune donnée dans la plage sélectionnée")
  }
  
  return(selected_data)
}

# Calcul du seuil relatif
calculate_relative_threshold <- function(data, n_samples = 250) {
  n_available <- min(n_samples, nrow(data))
  baseline_values <- mean(data$Fz_N[1:n_available])
  relative_threshold <- mean(data$Fz_N[1:n_available]) * 0.1 # facteur 2
  return(relative_threshold)
}


# Détection du début du mouvement
detect_movement_start <- function(data, body_weight, relative_threshold,
                                  window_size = DEFAULT_MOVEMENT_WINDOW,
                                  consecutive_points = DEFAULT_MOVEMENT_CONSECUTIVE) {
  
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
      if (i == 1) {
        start_index <- 1
      } else {
        start_index <- cumsum_lengths[i - 1] + 1
      }
      
      movement_start_time <- data$abs_time_s[start_index]
      return(movement_start_time)
    }
  }
  
  start_index <- which(deviation_to_use > relative_threshold)[1]
  
  if (is.na(start_index)) {
    return(NA)
  }
  
  movement_start_time <- data$abs_time_s[start_index]
  return(movement_start_time)
}

# Détection du décollage
detect_takeoff <- function(data, movement_start_time, body_weight, 
                           threshold_pct = DEFAULT_TAKEOFF_LANDING_THRESHOLD) {
  
  if (is.na(movement_start_time)) {
    return(NA)
  }
  
  takeoff_threshold <- body_weight * threshold_pct
  
  post_movement_data <- data %>%
    filter(abs_time_s >= movement_start_time)
  
  if (nrow(post_movement_data) == 0) {
    return(NA)
  }
  
  takeoff_index <- which(post_movement_data$Fz_N < takeoff_threshold)[1]
  
  if (is.na(takeoff_index)) {
    return(NA)
  }
  
  takeoff_time <- post_movement_data$abs_time_s[takeoff_index]
  return(takeoff_time)
}

# Détection de l'atterrissage
detect_landing <- function(data, takeoff_time, body_weight,
                           threshold_pct = DEFAULT_TAKEOFF_LANDING_THRESHOLD) {
  
  if (is.na(takeoff_time)) {
    return(NA)
  }
  
  landing_threshold <- body_weight * threshold_pct
  
  post_takeoff_data <- data %>%
    filter(abs_time_s > takeoff_time)
  
  if (nrow(post_takeoff_data) == 0) {
    return(NA)
  }
  
  landing_index <- which(post_takeoff_data$Fz_N > landing_threshold)[1]
  
  if (is.na(landing_index)) {
    return(NA)
  }
  
  landing_time <- post_takeoff_data$abs_time_s[landing_index]
  return(landing_time)
}

# Validation des résultats
validate_results <- function(movement_start_time, takeoff_time, landing_time,
                             push_time, flight_time) {
  messages <- character()
  is_valid <- TRUE
  
  if (any(is.na(c(movement_start_time, takeoff_time, landing_time)))) {
    messages <- c(messages, "Un ou plusieurs événements non détectés")
    is_valid <- FALSE
  } else {
    if (movement_start_time >= takeoff_time) {
      messages <- c(messages, "Ordre chronologique incorrect: début >= décollage")
      is_valid <- FALSE
    }
    
    if (takeoff_time >= landing_time) {
      messages <- c(messages, "Ordre chronologique incorrect: décollage >= atterrissage")
      is_valid <- FALSE
    }
    
    if (!is.na(flight_time)) {
      if (flight_time < MIN_FLIGHT_TIME) {
        messages <- c(messages, sprintf("Temps de vol trop court: %.0f ms", flight_time * 1000))
        is_valid <- FALSE
      } else if (flight_time > MAX_FLIGHT_TIME) {
        messages <- c(messages, sprintf("Temps de vol trop long: %.0f ms", flight_time * 1000))
        is_valid <- FALSE
      }
    }
    
    if (!is.na(push_time)) {
      if (push_time < MIN_PUSH_TIME) {
        messages <- c(messages, sprintf("Temps de poussée trop court: %.0f ms", push_time * 1000))
        is_valid <- FALSE
      } else if (push_time > MAX_PUSH_TIME) {
        messages <- c(messages, sprintf("Temps de poussée trop long: %.0f ms", push_time * 1000))
        is_valid <- FALSE
      }
    }
  }
  
  if (is_valid) {
    messages <- c(messages, "Tous les critères de validation sont respectés")
  }
  
  return(list(is_valid = is_valid, messages = messages))
}

# Fonction principale d'analyse
analyze_jump_data_with_selection <- function(selected_data, body_mass,
                                             takeoff_landing_threshold = DEFAULT_TAKEOFF_LANDING_THRESHOLD,
                                             movement_window_size = DEFAULT_MOVEMENT_WINDOW,
                                             movement_consecutive_points = DEFAULT_MOVEMENT_CONSECUTIVE) {
  
  tryCatch({
    body_weight <- calculate_body_weight(body_mass)
    relative_threshold <- calculate_relative_threshold(selected_data)
    
    movement_start_time <- detect_movement_start(selected_data, body_weight, relative_threshold,
                                                 movement_window_size, movement_consecutive_points)
    
    takeoff_time <- detect_takeoff(selected_data, movement_start_time, body_weight, takeoff_landing_threshold)
    landing_time <- detect_landing(selected_data, takeoff_time, body_weight, takeoff_landing_threshold)
    
    push_time <- if (!is.na(movement_start_time) && !is.na(takeoff_time)) {
      takeoff_time - movement_start_time
    } else NA
    
    flight_time <- if (!is.na(takeoff_time) && !is.na(landing_time)) {
      landing_time - takeoff_time
    } else NA
    
    validation <- validate_results(movement_start_time, takeoff_time, landing_time, push_time, flight_time)
    
    result <- JumpAnalysisResult(
      body_mass = body_mass,
      body_weight = body_weight,
      movement_start_time = movement_start_time,
      takeoff_time = takeoff_time,
      landing_time = landing_time,
      push_time = push_time,
      flight_time = flight_time,
      relative_threshold = relative_threshold,
      is_valid = validation$is_valid,
      validation_messages = validation$messages
    )
    
    return(result)za
    
  }, error = function(e) {
    return(JumpAnalysisResult(
      body_mass = body_mass,
      is_valid = FALSE,
      validation_messages = paste("Erreur d'analyse:", e$message)
    ))
  })
}


analyze_jump_data_with_selection <- function(selected_data, body_mass,
                                             takeoff_landing_threshold = DEFAULT_TAKEOFF_LANDING_THRESHOLD,
                                             movement_window_size = DEFAULT_MOVEMENT_WINDOW,
                                             movement_consecutive_points = DEFAULT_MOVEMENT_CONSECUTIVE) {
  tryCatch({
    body_weight <- calculate_body_weight(body_mass)
    relative_threshold <- calculate_relative_threshold(selected_data)
    
    # Détection des événements
    movement_start_time <- detect_movement_start(selected_data, body_weight, relative_threshold,
                                                 movement_window_size, movement_consecutive_points)
    takeoff_time <- detect_takeoff(selected_data, movement_start_time, body_weight, takeoff_landing_threshold)
    landing_time <- detect_landing(selected_data, takeoff_time, body_weight, takeoff_landing_threshold)
    
    # Temps caractéristiques
    push_time <- if (!is.na(movement_start_time) && !is.na(takeoff_time)) takeoff_time - movement_start_time else NA
    flight_time <- if (!is.na(takeoff_time) && !is.na(landing_time)) landing_time - takeoff_time else NA
    
    # === Calculs dynamiques ===
    if (!is.na(movement_start_time) && !is.na(takeoff_time)) {
      push_data <- selected_data %>%
        filter(abs_time_s >= movement_start_time & abs_time_s <= takeoff_time)
      
      dt <- mean(diff(push_data$abs_time_s))
      
      # Forces moyennes
      mean_force_push <- mean(push_data$Fz_N)
      mean_force_push_rel <- mean_force_push / body_mass
      
      # Forces moyennes nettes (sans le poids)
      mean_net_force_push <- mean(push_data$Fz_N - body_weight)
      mean_net_force_push_rel <- mean_net_force_push / body_mass
      
      # Impulsion nette
      net_impulse <- sum(push_data$Fz_N - body_weight) * dt
      
      # Vitesse de décollage
      takeoff_velocity <- net_impulse / body_mass
      
      # Hauteur de saut à partir de la vitesse
      jump_height <- (takeoff_velocity^2) / (2 * GRAVITY)
    } else {
      mean_force_push <- mean_force_push_rel <- NA
      mean_net_force_push <- mean_net_force_push_rel <- NA
      net_impulse <- takeoff_velocity <- jump_height <- NA
    }
    
    # Validation standard
    validation <- validate_results(movement_start_time, takeoff_time, landing_time, push_time, flight_time)
    
    # Création de l’objet résultat complet
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
      net_impulse = net_impulse,
      takeoff_velocity = takeoff_velocity,
      jump_height = jump_height,
      is_valid = validation$is_valid,
      validation_messages = validation$messages
    )
    
    return(result)
    
  }, error = function(e) {
    return(JumpAnalysisResult(
      body_mass = body_mass,
      is_valid = FALSE,
      validation_messages = paste("Erreur d'analyse:", e$message)
    ))
  })
}