library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
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
    # Formule de Samozino et al. (2012)
    # Vitesse de décollage théorique optimale
    # Vdec_opt = sqrt(Pmax × hpo / body_mass)
    
    Vdec_opt <- sqrt((Pmax * hpo) / body_mass)
    
    # Pente optimale de la relation force-vitesse (Slope Force-Velocity optimal)
    # SFVopt = -(Pmax / Vdec_opt²)
    # ou de façon équivalente : SFVopt = -body_mass / hpo
    
    SFVopt <- -(body_mass * GRAVITY) / hpo
    
    # Vérifications de cohérence
    if (Vdec_opt < 0 || Vdec_opt > 5) {
      warning(paste("Vdec_opt incohérent:", round(Vdec_opt, 2), "m/s"))
      return(list(Vdec_opt = NA, SFVopt = NA))
    }
    
    if (abs(SFVopt) < 100 || abs(SFVopt) > 10000) {
      warning(paste("SFVopt incohérent:", round(SFVopt, 1), "N/(m/s)"))
      return(list(Vdec_opt = NA, SFVopt = NA))
    }
    
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

# ===== UI =====
ui <- dashboardPage(
  dashboardHeader(title = "Analyse Saut + Vitesse"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Analyse", tabName = "analysis", icon = icon("chart-line")),
      menuItem("Vitesse", tabName = "velocity", icon = icon("tachometer-alt")),
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
  
  observeEvent(input$reset_params, {
    updateNumericInput(session, "movement_window", value = 10)
    updateNumericInput(session, "consecutive_points", value = 50)
    updateNumericInput(session, "takeoff_threshold", value = 0.10)
  })
  
  observeEvent(input$file, {
    req(input$file)
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
      removeNotification("analyzing")
      showNotification("Analyse terminée!", type = "success", duration = 3)
    }, error = function(e) {
      removeNotification("analyzing")
      showNotification(paste("Erreur:", e$message), type = "error")
    })
  })
  
  output$file_loaded <- reactive(!is.null(original_data()))
  outputOptions(output, "file_loaded", suspendWhenHidden = FALSE)
  
  output$range_selected <- reactive(!is.null(selected_data()))
  outputOptions(output, "range_selected", suspendWhenHidden = FALSE)
  
  output$analysis_done <- reactive(!is.null(analysis_results()))
  outputOptions(output, "analysis_done", suspendWhenHidden = FALSE)
  
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