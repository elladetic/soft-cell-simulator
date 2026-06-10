library(shiny)
library(bslib)
library(ggplot2)

# --- Algorithm ---

# One iteration step on a single chain (matrix with columns x, y).
# Chain starts at P+ and ends at P-.
cut_chain_step <- function(chain, rho) {
  n <- nrow(chain)
  if (n < 3) return(chain)

  P_plus  <- chain[1, ]
  P_minus <- chain[n, ]
  new_chain <- matrix(P_plus, nrow = 1)

  for (i in 1:(n - 2)) {
    a_prev <- chain[i, ]
    a_curr <- chain[i + 1, ]
    a_next <- chain[i + 2, ]

    l_i <- rho * a_prev + (1 - rho) * a_curr
    r_i <- (1 - rho) * a_curr + rho * a_next

    if (i == 1) {
      proj_l <- c(P_plus[1], l_i[2])
      l_i <- (1 - rho) * l_i + rho * proj_l
    }
    if (i == n - 2) {
      proj_r <- c(P_minus[1], r_i[2])
      r_i <- (1 - rho) * r_i + rho * proj_r
    }

    new_chain <- rbind(new_chain, l_i, r_i)
  }
  rbind(new_chain, matrix(P_minus, nrow = 1))
}

# Returns a data frame of cut points for a chain (for visualisation).
# Columns: x, y, type ("l_star" / "r_star" for boundary-corrected, "l" / "r" for interior)
get_cut_points <- function(chain, rho) {
  n <- nrow(chain)
  if (n < 3) return(NULL)

  P_plus  <- chain[1, ]
  P_minus <- chain[n, ]
  rows <- list()

  for (i in 1:(n - 2)) {
    a_prev <- chain[i, ]
    a_curr <- chain[i + 1, ]
    a_next <- chain[i + 2, ]

    l_i <- rho * a_prev + (1 - rho) * a_curr
    r_i <- (1 - rho) * a_curr + rho * a_next

    l_type <- "l"
    r_type <- "r"

    if (i == 1) {
      proj_l <- c(P_plus[1], l_i[2])
      l_i    <- (1 - rho) * l_i + rho * proj_l
      l_type <- "l_star"
    }
    if (i == n - 2) {
      proj_r <- c(P_minus[1], r_i[2])
      r_i    <- (1 - rho) * r_i + rho * proj_r
      r_type <- "r_star"
    }

    rows[[length(rows) + 1]] <- c(l_i, l_type)
    rows[[length(rows) + 1]] <- c(r_i, r_type)
  }

  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(df) <- c("x", "y", "type")
  df$x <- as.numeric(df$x)
  df$y <- as.numeric(df$y)
  df
}

# Run n_iter steps and return all intermediate states
run_iterations <- function(A0, B0, rho, n_iter) {
  states <- vector("list", n_iter + 1)
  states[[1]] <- list(A = A0, B = B0)
  A <- A0; B <- B0
  for (k in seq_len(n_iter)) {
    A <- cut_chain_step(A, rho)
    B <- cut_chain_step(B, rho)
    states[[k + 1]] <- list(A = A, B = B)
  }
  states
}

# Build a closed polygon data frame from chains A and B
chains_to_polygon <- function(A, B) {
  pts <- rbind(A, B[nrow(B):1, ])
  data.frame(x = pts[, 1], y = pts[, 2])
}

# --- Initial polygon: diamond (rotated square) ---
# P+ = (0,1), P- = (0,-1); right chain A, left chain B
A0 <- matrix(c(0, 1,   1, 0,   0, -1), ncol = 2, byrow = TRUE)
B0 <- matrix(c(0, 1,  -1, 0,   0, -1), ncol = 2, byrow = TRUE)

P_plus_coords  <- data.frame(x = 0, y =  1, label = "P⁺")
P_minus_coords <- data.frame(x = 0, y = -1, label = "P⁻")

# --- UI ---
ui <- page_sidebar(
  title = "Iterativna konstrukcija mekih ćelija",
  theme = bs_theme(version = 5, preset = "shiny"),
  sidebar = sidebar(
    sliderInput("rho", "Parametar ρ",
                min = 0.05, max = 0.45, value = 0.25, step = 0.01),
    sliderInput("n_iter", "Broj iteracija k",
                min = 0, max = 12, value = 0, step = 1,
                animate = animationOptions(interval = 600, loop = FALSE)),
    hr(),
    input_switch("show_initial",    "Prikaži početni poligon",   value = TRUE),
    input_switch("show_vertices",   "Prikaži vrhove lanca",       value = FALSE),
    input_switch("show_cut_points", "Prikaži rezne točke",        value = FALSE)
  ),
  card(
    full_screen = TRUE,
    card_header("Simulacija"),
    plotOutput("plot", height = "520px")
  )
)

# --- Server ---
server <- function(input, output, session) {

  states <- reactive({
    run_iterations(A0, B0, input$rho, 12)
  })

  output$plot <- renderPlot({
    k     <- input$n_iter + 1
    state <- states()[[k]]
    A     <- state$A
    B     <- state$B

    df_cell <- chains_to_polygon(A, B)

    p <- ggplot(df_cell, aes(x, y)) +
      geom_polygon(fill = "steelblue", alpha = 0.2, color = "black", linewidth = 0.9) +
      coord_equal(xlim = c(-1.5, 1.5), ylim = c(-1.4, 1.4)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            legend.title    = element_blank()) +
      labs(x = NULL, y = NULL,
           title = paste0("k = ", input$n_iter, "   |   ρ = ", input$rho))

    # Početni poligon (isprekidano)
    if (input$show_initial) {
      df_init <- chains_to_polygon(A0, B0)
      p <- p + geom_polygon(
        data = df_init, aes(x, y),
        fill = NA, color = "gray60", linetype = "dashed", linewidth = 0.6
      )
    }

    # Unutarnji vrhovi tekućih lanaca
    if (input$show_vertices) {
      # inner vertices only (exclude P+ and P-)
      inner_A <- as.data.frame(A[2:(nrow(A) - 1), , drop = FALSE])
      inner_B <- as.data.frame(B[2:(nrow(B) - 1), , drop = FALSE])
      names(inner_A) <- names(inner_B) <- c("x", "y")
      df_inner <- rbind(inner_A, inner_B)

      if (nrow(df_inner) > 0) {
        p <- p + geom_point(
          data = df_inner, aes(x, y),
          shape = 21, fill = "white", color = "black", size = 2.5
        )
      }
    }

    # Rezne točke za sljedeći korak
    if (input$show_cut_points) {
      cp_A <- get_cut_points(A, input$rho)
      cp_B <- get_cut_points(B, input$rho)

      if (!is.null(cp_A) && !is.null(cp_B)) {
        df_cp <- rbind(cp_A, cp_B)
        # Reclassify for legend: boundary-corrected vs. interior
        df_cp$vrsta <- ifelse(
          df_cp$type %in% c("l_star", "r_star"),
          "Rubna korekcija",
          "Rezna točka"
        )

        p <- p + geom_point(
          data = df_cp, aes(x, y, fill = vrsta),
          shape = 21, color = "white", size = 3.5
        ) +
          scale_fill_manual(values = c(
            "Rezna točka"       = "#E07B39",
            "Rubna korekcija" = "#9B2335"
          ))
      }
    }

    # P+ i P- — uvijek vidljivi
    p <- p +
      geom_point(data = P_plus_coords,  aes(x, y),
                 shape = 23, fill = "#2ecc71", color = "black", size = 5) +
      geom_point(data = P_minus_coords, aes(x, y),
                 shape = 23, fill = "#e74c3c", color = "black", size = 5) +
      geom_text(data = P_plus_coords,  aes(x, y, label = label),
                nudge_x = 0.12, nudge_y = 0.08, size = 4.5, fontface = "bold") +
      geom_text(data = P_minus_coords, aes(x, y, label = label),
                nudge_x = 0.12, nudge_y = -0.08, size = 4.5, fontface = "bold")

    p
  })
}

shinyApp(ui, server)
