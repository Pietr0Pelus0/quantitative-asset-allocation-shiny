# ============================================================
# APPLICATIVO: Asset Allocation Dashboard
# ============================================================
# Descrizione: Dashboard interattiva per il confronto di diverse 
# strategie di allocazione:
# Minima Varianza, Tangency, ERC, Concentrazione, Target Vol,
# Equal Weight. Integra analisi CAPM per ciascuna strategia.
#
# ------------------------------------------------------------
# AUTORE: Pietro Peluso
# Data: 13/04/2026
# Note: All'esecuzione dell'app è affiancata una relazione tecnica che approfondisce
# e contestualizza i risultati prodotti dalla dashboard, analizzandone le
# implicazioni finanziarie e valutandone la coerenza rispetto al profilo
# di un investitore con elevata tolleranza al rischio e orizzonte di
# investimento di lungo periodo.
# ============================================================


# ------------------------------------------------------------
# SEZIONE 1: LIBRERIE
# ------------------------------------------------------------
# Caricamento dei pacchetti necessari.
# suppressPackageStartupMessages() evita output verbose
# durante l'avvio (utile in ambienti Shiny).
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)                
  library(quantmod)             
  library(xts)                  
  library(PerformanceAnalytics) 
  library(quadprog)             
  library(DT)                   
  library(ggplot2)              
  library(corrplot)             
  library(tidyr)                
})


# ------------------------------------------------------------
# SEZIONE 2: FUNZIONI ANALITICHE (Helper Functions)
# ------------------------------------------------------------
# Contiene tutte le funzioni core utilizzate dal server.
# Struttura:
# 2a -> Frontiera efficiente
# 2b -> Strategie portafoglio
# 2c -> Grafici
# 2d -> Analisi rischio
# ------------------------------------------------------------


# ------------------------------------------------------------
# 2a. FRONTIERA EFFICIENTE (Long-Only)
# ------------------------------------------------------------
# Calcola la frontiera efficiente imponendo:
#   - somma pesi = 1
#   - pesi >= 0 (no short selling)
#
# INPUT:
#   mu     -> vettore rendimenti attesi
#   Sigma  -> matrice covarianza
#   n_points -> numero punti frontiera
#
# OUTPUT:
#   list con:
#     df      -> (ret, vol)
#     weights -> matrice pesi
#     w_mv    -> pesi min var
#     ret_mv  -> rendimento min var
#     vol_mv  -> volatilità min var
# ------------------------------------------------------------
efficient_frontier_long_only <- function(mu, Sigma, n_points = 120) {
  n_assets <- length(mu)
  
  # Vincoli base:
  # - somma pesi = 1
  # - pesi >= 0
  Amat <- cbind(rep(1, n_assets), diag(n_assets))
  bvec <- c(1, rep(0, n_assets))
  
  # --- Portafoglio a Minima Varianza ---
  w_min_var_sol <- solve.QP(2 * Sigma, rep(0, n_assets), Amat, bvec, meq = 1)
  
  # Pulizia numerica + normalizzazione
  w_min_var <- pmax(w_min_var_sol$solution, 0)
  w_min_var <- w_min_var / sum(w_min_var)
  names(w_min_var) <- names(mu)
  
  mu_min_var  <- sum(w_min_var * mu)
  vol_min_var <- sqrt(drop(t(w_min_var) %*% Sigma %*% w_min_var))
  
  # Range rendimenti target
  mu_max_ret     <- max(mu)
  target_returns <- seq(mu_min_var, mu_max_ret, length.out = n_points)
  
  results        <- data.frame(ret = numeric(n_points), vol = numeric(n_points))
  weights_matrix <- matrix(0, nrow = n_points, ncol = n_assets)
  
  # --- Costruzione frontiera ---
  for (i in 1:n_points) {
    Amat_i <- cbind(mu, rep(1, n_assets), diag(n_assets))
    bvec_i <- c(target_returns[i], 1, rep(0, n_assets))
    
    sol <- tryCatch(
      solve.QP(2 * Sigma, rep(0, n_assets), Amat_i, bvec_i, meq = 2),
      error = function(e) NULL
    )
    
    if (!is.null(sol)) {
      results$ret[i] <- target_returns[i]
      results$vol[i] <- sqrt(drop(t(sol$solution) %*% Sigma %*% sol$solution))
      
      # Clamp numerico pesi negativi
      weights_matrix[i, ] <- pmax(0, sol$solution)
    }
  }
  
  colnames(weights_matrix) <- names(mu)
  
  # Filtra eventuali soluzioni non valide
  valid <- results$vol > 0
  
  list(
    df      = results[valid, ],
    weights = weights_matrix[valid, ],
    w_mv    = w_min_var,
    ret_mv  = mu_min_var,
    vol_mv  = vol_min_var
  )
}


# ------------------------------------------------------------
# 2b. STRATEGIE DI PORTAFOGLIO
# ------------------------------------------------------------

# --- 1. Minimum Variance ---
# Riutilizza risultato già calcolato (se disponibile)
# per evitare costi computazionali duplicati.
get_min_var <- function(mu, Sigma, ef = NULL) {
  
  if (!is.null(ef) && !is.null(ef$w_mv)) {
    return(list(w = ef$w_mv, ret = ef$ret_mv, vol = ef$vol_mv))
  }
  
  # Fallback: risoluzione diretta
  n    <- length(mu)
  Amat <- cbind(rep(1, n), diag(n))
  bvec <- c(1, rep(0, n))
  
  sol <- solve.QP(2 * Sigma, rep(0, n), Amat, bvec, meq = 1)
  
  w <- pmax(sol$solution, 0)
  w <- w / sum(w)
  names(w) <- names(mu)
  
  list(
    w   = w,
    ret = sum(w * mu),
    vol = sqrt(drop(t(w) %*% Sigma %*% w))
  )
}


# --- 2. Tangency Portfolio (Max Sharpe) ---
# Metodo QP con rinormalizzazione finale.
get_tangency <- function(mu, Sigma, rf = 0) {
  k <- length(mu)
  if (k < 2) return(NULL)
  
  # Assicura naming coerente
  if (is.null(names(mu)) || any(names(mu) == "")) {
    names(mu) <- paste0("Asset_", seq_len(k))
  }
  
  ex <- mu - rf
  
  # Controllo fattibilità
  if (all(!is.finite(ex)) || all(ex <= 0)) return(NULL)
  
  Dmat <- 2 * Sigma
  dvec <- rep(0, k)
  
  Amat <- cbind(ex, diag(k))
  bvec <- c(1, rep(0, k))
  
  sol <- tryCatch(
    solve.QP(Dmat, dvec, Amat, bvec, meq = 1),
    error = function(e) NULL
  )
  
  if (is.null(sol)) return(NULL)
  
  w <- pmax(sol$solution, 0)
  if (sum(w) <= 0) return(NULL)
  
  w <- w / sum(w)
  names(w) <- names(mu)
  
  ret <- sum(w * mu)
  vol <- sqrt(drop(t(w) %*% Sigma %*% w))
  
  sharpe <- if (is.finite(vol) && vol > 0) (ret - rf) / vol else NA_real_
  
  list(w = w, ret = ret, vol = vol, sharpe = sharpe)
}


# --- 3. Equal Risk Contribution (ERC) ---
# Risk parity tramite ottimizzazione convessa.
get_erc <- function(mu, Sigma) {
  N <- ncol(Sigma)
  
  # Funzione obiettivo
  eval_f <- function(x) { 
    0.5 * as.numeric(t(x) %*% Sigma %*% x) - sum(log(x)) 
  }
  
  # Gradiente analitico (fondamentale per stabilità numerica)
  eval_grad <- function(x) { 
    as.vector(Sigma %*% x) - 1/x 
  }
  
  x0 <- rep(1/N, N)
  
  res <- optim(
    par    = x0, 
    fn     = eval_f, 
    gr     = eval_grad, 
    method = "L-BFGS-B", 
    lower  = rep(1e-6, N)
  )
  
  # Normalizzazione finale
  w_final <- res$par / sum(res$par)
  names(w_final) <- colnames(Sigma)
  
  list(
    w   = w_final,
    ret = sum(w_final * mu),
    vol = sqrt(as.numeric(t(w_final) %*% Sigma %*% w_final))
  )
}


# --- 4. Portafoglio con vincolo di concentrazione ---
# Limita peso massimo per asset.
get_concentrated <- function(mu, Sigma, user_limit) {
  n <- length(mu)
  
  user_limit_dec <- user_limit / 100
  
  # Evita infeasibility matematica
  limit <- max(user_limit_dec, (1/n) + 0.001)
  
  Amat <- cbind(rep(1, n), diag(n), -diag(n))
  bvec <- c(1, rep(0, n), rep(-limit, n))
  
  sol <- try(solve.QP(2 * Sigma, rep(0, n), Amat, bvec, meq = 1), silent = TRUE)
  
  if (inherits(sol, "try-error")) {
    w <- rep(1/n, n)
  } else {
    w <- sol$solution
  }
  
  names(w) <- names(mu)
  
  # Pulizia numerica
  w[w < 1e-7] <- 0
  w <- w / sum(w)
  
  list(
    w          = w,
    ret        = sum(w * mu),
    vol        = sqrt(drop(t(w) %*% Sigma %*% w)),
    limit_used = limit
  )
}


# --- 5. Equal Weight ---
# Benchmark base (naive diversification)
get_equal_weight <- function(mu, Sigma) {
  n <- length(mu)
  w <- rep(1/n, n)
  names(w) <- names(mu)
  
  list(
    w   = w,
    ret = sum(w * mu),
    vol = sqrt(drop(t(w) %*% Sigma %*% w))
  )
}


# ------------------------------------------------------------
# 2c. FUNZIONI GRAFICHE
# ------------------------------------------------------------

# --- Barplot pesi ---
plot_weights_bar <- function(weights) {
  df <- data.frame(Asset = names(weights), Peso = as.numeric(weights))
  
  ggplot(df, aes(x = reorder(Asset, -Peso), y = Peso, fill = Asset)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    labs(x = "Asset", y = "Peso (%)") +
    theme(legend.position = "none")
}


# --- Frontiera + punto portafoglio ---
plot_frontier_point <- function(ef_df, p_ret, p_vol) {
  ggplot(ef_df, aes(x = vol, y = ret)) +
    geom_line(color = "darkgrey", size = 1) +
    annotate("point", x = p_vol, y = p_ret, color = "red", size = 4) +
    theme_minimal() +
    labs(x = "Volatilità", y = "Rendimento Atteso")
}


# ------------------------------------------------------------
# 2d. ANALISI DEL RISCHIO
# ------------------------------------------------------------

# --- Risk Decomposition ---
get_risk_decomposition <- function(w, Sigma) {
  w <- as.matrix(w)
  
  port_vol <- as.numeric(sqrt(t(w) %*% Sigma %*% w))
  
  # Marginal Risk Contribution
  mrc <- (Sigma %*% w) / port_vol
  
  # Percentage Risk Contribution (somma = 1)
  prc <- (w * mrc) / port_vol
  
  # Dominanza rischio vs peso
  dom <- mrc / port_vol
  
  data.frame(
    MRC       = as.numeric(mrc),
    PRC       = as.numeric(prc),
    Dominance = as.numeric(dom)
  )
}


# --- Scatter Peso vs Rischio ---
plot_risk_scatter <- function(w, prc, title) {
  df_plot <- data.frame(
    Ticker  = names(w),
    Peso    = as.numeric(w),
    Rischio = as.numeric(prc)
  )
  
  ax_max <- max(df_plot$Peso, df_plot$Rischio) * 1.1
  
  ggplot(df_plot, aes(x = Peso, y = Rischio, label = Ticker)) +
    geom_point(color = "royalblue", size = 3) +
    geom_text(vjust = -1, size = 3.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_x_continuous(labels = scales::percent, limits = c(0, ax_max)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, ax_max)) +
    labs(
      title    = title,
      subtitle = "Equilibrio Peso-Rischio (ERC ideale sulla diagonale)",
      x        = "Peso (%)",
      y        = "Rischio (%)"
    ) +
    theme_minimal()
}


# ============================================================
# SEZIONE 3: INTERFACCIA UTENTE (UI)
# ============================================================
# Definisce il layout visivo dell'applicazione Shiny.
#
# Struttura generale:
#   - Sidebar: input utente (configurazione analisi)
#   - Main Panel: output organizzati in tab tematiche
#
# Logica UX:
#   - L'utente definisce universo e parametri
#   - Clicca "Esegui Analisi"
#   - Naviga tra le strategie per confronto
# ============================================================


ui <- fluidPage(
  
  # ----------------------------------------------------------
  # HEADER APPLICAZIONE
  # ----------------------------------------------------------
  # Nota: titolo statico, eventualmente rendibile dinamico
  # con versione/app info in futuro
  # ----------------------------------------------------------
  titlePanel("Asset Allocation Dashboard - Pietro Peluso"),
  
  
  sidebarLayout(
    
    # ==========================================================
    # SIDEBAR PANEL
    # ==========================================================
    # Contiene tutti gli input controllabili dall'utente.
    # Le modifiche NON triggerano automaticamente i calcoli:
    # è necessario premere il bottone "run".
    # ==========================================================
    sidebarPanel(
      
      helpText("Configurazione Universo Investibile"),
      
      # --------------------------------------------------------
      # INPUT: TICKERS
      # --------------------------------------------------------
      # Lista asset separati da virgola.
      # Parsing gestito lato server (split + trim).
      # Attenzione: ticker non validi generano errori quantmod.
      # --------------------------------------------------------
      textInput(
        "tickers",
        "Tickers (separati da virgola):",
        value = "AAPL,MSFT,ASML,JPM,HSBC,JNJ,NVO,XOM,TTE,PG,BABA,HD,RIO,NEE,CAT"
      ),
      
      
      # --------------------------------------------------------
      # INPUT: RANGE TEMPORALE
      # --------------------------------------------------------
      # Definisce finestra per download dati storici.
      # Nota: Yahoo Finance può avere dati mancanti su alcuni asset.
      # --------------------------------------------------------
      dateInput("from", "Data Inizio:", value = "2023-03-21"),
      dateInput("to",   "Data Fine:",   value = "2026-03-20"),
      
      hr(),
      
      
      # --------------------------------------------------------
      # INPUT: VINCOLO CONCENTRAZIONE
      # --------------------------------------------------------
      # Limite massimo peso per asset (%).
      # Usato nella strategia "Concentrated".
      # --------------------------------------------------------
      sliderInput(
        "max_w_slider",
        "Vincolo Max Concentrazione (%):", 
        min = 5, max = 100, value = 16, step = 1,
        post = "%"
      ),
      
      # Info dinamica su limite minimo teorico (1/N)
      tags$small(uiOutput("min_limit_info")),
      
      
      # --------------------------------------------------------
      # INPUT: TARGET VOLATILITY
      # --------------------------------------------------------
      # Volatilità target annualizzata (%).
      # Il portafoglio selezionato sarà il punto della frontiera
      # più vicino a questo valore.
      # --------------------------------------------------------
      sliderInput(
        "target_vol",
        "Scegli Volatilità Target (Annua %):",
        min = 5, max = 50, value = 17, step = 1, 
        post = "%"
      ),
      
      helpText("La strategia seleziona il portafoglio più vicino alla volatilità target."),
      
      
      # --------------------------------------------------------
      # INPUT: RISK-FREE RATE
      # --------------------------------------------------------
      # Ticker del proxy risk-free.
      # Esempio: ^IRX (13-week T-Bill).
      # Conversione a rendimento gestita lato server.
      # --------------------------------------------------------
      textInput(
        "rf_ticker",
        "Ticker Risk-Free (es. ^IRX)",
        value = "^IRX"
      ),
      
      
      # --------------------------------------------------------
      # INPUT: BENCHMARK DI MERCATO
      # --------------------------------------------------------
      # Usato per:
      #   - regressione CAPM
      #   - calcolo beta e R^2
      # --------------------------------------------------------
      textInput(
        "mkt_ticker",
        "Ticker Benchmark Mercato (es. ^GSPC)",
        value = "^GSPC"
      ),
      
      
      # --------------------------------------------------------
      # ACTION BUTTON
      # --------------------------------------------------------
      # Trigger esplicito dell'analisi.
      # Best practice Shiny: evita ricalcoli ad ogni input change.
      # --------------------------------------------------------
      actionButton(
        "run",
        "Esegui Analisi",
        class = "btn-primary",
        width = "100%"
      )
    ),
    
    
    # ==========================================================
    # MAIN PANEL
    # ==========================================================
    # Contiene gli output suddivisi per tab.
    # Ogni tab rappresenta:
    #   - una vista dati
    #   - oppure una strategia di portafoglio
    # ==========================================================
    mainPanel(
      
      tabsetPanel(
        
        # ----------------------------------------------------------
        # TAB 1: UNIVERSO INVESTIBILE
        # ----------------------------------------------------------
        # Scopo: esplorazione dati input
        # ----------------------------------------------------------
        tabPanel(
          "Universo investibile",
          
          fluidRow(
            column(6,
                   plotOutput("pricePlot")  # Prezzi adjusted
            ),
            column(6,
                   plotOutput("retPlot")    # Returns log o semplici
            )
          ),
          
          hr(),
          
          h4("Matrice di Correlazione"),
          plotOutput("corrPlot", height = "500px"),
          
          hr(),
          
          # Statistiche riassuntive:
          # - rendimento medio annualizzato
          # - volatilità annualizzata
          DTOutput("statsTable")
          
          # --- BLOCCO SPERIMENTALE (disattivato) ---
          # plot confronto strategie (utile per futura estensione UX)
          # hr(),
          # h3("Confronto Strategie"),
          # plotOutput("comparisonPlot", height = "600px")
        ),
        
        
        # ==========================================================
        # TEMPLATE STRATEGIE
        # ==========================================================
        # Struttura standard per ogni strategia:
        #   1. Frontiera + Barplot pesi
        #   2. Scatter Peso vs Rischio
        #   3. Tabella dettagli asset
        #   4. Performance + CAPM
        # ==========================================================
        
        
        # --- Minima Varianza ---
        tabPanel(
          "Minima Varianza",
          h3("Portafoglio a Minima Varianza"),
          fluidRow(
            column(6, plotOutput("mvFrontier")),
            column(6, plotOutput("mvBar"))
          ),
          hr(),
          plotOutput("mvScatter"),
          hr(),
          DTOutput("mvTable"),
          hr(),
          fluidRow(
            column(4,
                   h4("Performance Strategia", style = "font-weight: bold;"),
                   tableOutput("mvSummary")),
            column(4,
                   h4("CAPM di Portafoglio", style = "font-weight: bold;"),
                   tableOutput("mvCapm"))
          )
        ),
        
        
        # --- Tangency ---
        tabPanel(
          "Tangency",
          h3("Portafoglio di Tangenza"),
          fluidRow(
            column(6, plotOutput("tanFrontier")),
            column(6, plotOutput("tanBar"))
          ),
          hr(),
          plotOutput("tanScatter"),
          hr(),
          DTOutput("tanTable"),
          hr(),
          fluidRow(
            column(4,
                   h4("Performance Strategia", style = "font-weight: bold;"),
                   tableOutput("tanSummary")),
            column(4,
                   h4("CAPM di Portafoglio", style = "font-weight: bold;"),
                   tableOutput("tanCapm"))
          )
        ),
        
        
        # --- Equal Risk Contribution ---
        tabPanel(
          "Equal Risk Contribution",
          h3("Portafoglio con Equal Risk Contribution"),
          fluidRow(
            column(6, plotOutput("ercFrontier")),
            column(6, plotOutput("ercBar"))
          ),
          hr(),
          plotOutput("ercScatter"),
          hr(),
          DTOutput("ercTable"),
          hr(),
          fluidRow(
            column(4,
                   h4("Performance Strategia", style = "font-weight: bold;"),
                   tableOutput("ercSummary")),
            column(4,
                   h4("CAPM di Portafoglio", style = "font-weight: bold;"),
                   tableOutput("ercCapm"))
          )
        ),
        
        
        # --- Vincolo Concentrazione ---
        tabPanel(
          "Vincolo Concentrazione",
          h3("Portafoglio con Vincolo di Concentrazione"),
          fluidRow(
            column(6, plotOutput("concFrontier")),
            column(6, plotOutput("concBar"))
          ),
          hr(),
          plotOutput("concScatter"),
          hr(),
          DTOutput("concTable"),
          hr(),
          fluidRow(
            column(4,
                   h4("Performance Strategia", style = "font-weight: bold;"),
                   tableOutput("concSummary")),
            column(4,
                   h4("CAPM di Portafoglio", style = "font-weight: bold;"),
                   tableOutput("concCapm"))
          )
        ),
        
        
        # --- Target Volatility ---
        tabPanel(
          "Target Volatility",
          h3("Portafoglio con Target di Volatilità"),
          fluidRow(
            column(6, plotOutput("tvFrontier")),
            column(6, plotOutput("tvBar"))
          ),
          hr(),
          plotOutput("tvScatter"),
          hr(),
          DTOutput("tvTable"),
          hr(),
          fluidRow(
            column(4,
                   h4("Performance Strategia", style = "font-weight: bold;"),
                   tableOutput("tvSummary")),
            column(4,
                   h4("CAPM di Portafoglio", style = "font-weight: bold;"),
                   tableOutput("tvCapm"))
          )
        ),
        
        
        # --- Equal Weight ---
        tabPanel(
          "Equal Weight",
          h3("Portafoglio Equipesato"),
          fluidRow(
            column(6, plotOutput("ewFrontier")),
            column(6, plotOutput("ewBar"))
          ),
          hr(),
          plotOutput("ewScatter"),
          hr(),
          DTOutput("ewTable"),
          hr(),
          fluidRow(
            column(4,
                   h4("Performance Strategia", style = "font-weight: bold;"),
                   tableOutput("ewSummary")),
            column(4,
                   h4("CAPM di Portafoglio", style = "font-weight: bold;"),
                   tableOutput("ewCapm"))
          )
        )
      )
    )
  )
)



# ============================================================
# SEZIONE 4: SERVER
# Logica reattiva dell'applicazione Shiny.
#
# Struttura:
#   4a. Download e preprocessing dati (eventReactive)
#   4b. Calcolo frontiera efficiente (reactive)
#   4c. Preparazione dati per grafici Tab 1 (reactive)
#   4d. Output Tab 1: grafici esplorativi
#   4e. Output Tab 2-7: strategie (non incluse qui)
# ============================================================
server <- function(input, output, session) {
  
  # ----------------------------------------------------------
  # 4a. DOWNLOAD E PREPROCESSING DATI
  #
  # Trigger: click su input$run (eventReactive)
  #
  # Workflow:
  #   1. Parsing ticker inseriti dall'utente
  #   2. Download prezzi adjusted (asset + benchmark + RF)
  #   3. Calcolo rendimenti semplici
  #   4. Allineamento temporale delle serie
  #   5. Stima CAPM (Beta e R²) per ciascun asset
  #
  # Output:
  #   Lista contenente tutte le variabili necessarie
  #   per le successive elaborazioni reattive
  # ----------------------------------------------------------
  data_all <- eventReactive(input$run, {
    
    # ---- 1. Parsing ticker ----
    # Rimuove spazi e separa per virgola
    tickers <- unlist(strsplit(gsub(" ", "", input$tickers), ","))
    
    # ---- 2. Download prezzi asset ----
    # getSymbols salva gli oggetti nell'ambiente globale
    getSymbols(tickers, from = input$from, to = input$to, src = "yahoo", auto.assign = TRUE)
    
    # Estrazione prezzi adjusted e merge in un unico oggetto xts
    prices <- do.call(merge, lapply(tickers, function(x) Ad(get(x))))
    colnames(prices) <- tickers
    
    # ---- 3. Download benchmark di mercato ----
    mkt_prices <- Ad(getSymbols(
      input$mkt_ticker,
      from = input$from,
      to = input$to,
      auto.assign = FALSE
    ))
    colnames(mkt_prices) <- "Market"
    
    # ---- 4. Download risk-free ----
    # suppressWarnings evita warning su dati mancanti/interpolazioni
    rf_data <- suppressWarnings(
      getSymbols(input$rf_ticker, from = input$from, to = input$to, auto.assign = FALSE)
    )
    
    # Pulizia serie RF:
    # - na.locf: forward fill (tipico per tassi)
    # - na.omit: rimozione residui NA
    rf_data <- na.locf(rf_data, na.rm = FALSE)
    rf_data <- na.omit(rf_data)
    
    # Conversione:
    # da percentuale annua -> rendimento giornaliero decimale
    rf_daily_series <- Ad(rf_data) / 100 / 252
    colnames(rf_daily_series) <- "RF_Daily"
    
    
    # ---- 5. Calcolo rendimenti ----
    # NOTA: usiamo rendimenti semplici (non log)
    # perché compatibili con combinazioni lineari (Markowitz)
    returns        <- na.omit(Return.calculate(prices, method = "simple"))
    market_returns <- na.omit(Return.calculate(mkt_prices, method = "simple"))
    
    # ---- 6. Allineamento serie ----
    # Merge e rimozione NA per garantire coerenza temporale
    merged_all <- na.omit(merge(returns, market_returns, rf_daily_series))
    
    # Separazione componenti
    ret_assets <- merged_all[, 1:length(tickers)]
    ret_mkt    <- merged_all[, "Market"]
    rf_vec     <- merged_all[, "RF_Daily"]  # Serie RF giornaliera (per CAPM dinamico)
    
    # ---- 7. Risk-free scalare ----
    # Ultimo valore disponibile (annuo) -> utile per Sharpe ratio
    rf_rate  <- as.numeric(last(na.omit(Ad(rf_data)))) / 100
    rf_daily <- rf_rate / 252
    
    
    # ---- 8. Stima CAPM ----
    # Regressione: (R_i - Rf) ~ (Rm - Rf)
    capm_stats <- apply(ret_assets, 2, function(x) {
      
      # Excess returns
      excess_asset <- x - as.numeric(rf_vec)
      excess_mkt   <- ret_mkt - as.numeric(rf_vec)
      
      # Modello lineare
      model <- lm(excess_asset ~ excess_mkt)
      s <- summary(model)
      
      # Output: Beta e R²
      c(
        Beta = as.numeric(coef(model)[2]),
        R2   = s$r.squared
      )
    })
    
    
    # ---- 9. Output finale ----
    # Oggetto centralizzato per tutta la logica dell'app
    list(
      P          = prices,           # Prezzi adjusted (xts)
      R          = ret_assets,       # Rendimenti asset allineati
      mu         = colMeans(ret_assets),  # Rendimenti medi
      Sigma      = cov(ret_assets),       # Matrice covarianza
      rf         = rf_rate,               # RF annuo (scalare)
      Rm         = ret_mkt,               # Rendimenti mercato
      Rf_vec     = rf_vec,               # RF giornaliero (serie)
      stock_beta = capm_stats[1, ],      # Beta per asset
      stock_rsq  = capm_stats[2, ]       # R² per asset
    )
  })
  
  
  # ----------------------------------------------------------
  # 4b. CALCOLO FRONTIERA EFFICIENTE
  #
  # Reactive puro: si aggiorna automaticamente quando cambia
  # data_all()
  #
  # Usa 120 punti per ottenere una curva liscia
  # ----------------------------------------------------------
  markowitz_results <- reactive({
    d <- data_all()
    efficient_frontier_long_only(d$mu, d$Sigma, n_points = 120)
  })
  
  
  # ----------------------------------------------------------
  # 4c. TRASFORMAZIONE DATI PER GRAFICI
  #
  # Converte da formato xts (wide) a data.frame long
  # richiesto da ggplot2
  # ----------------------------------------------------------
  data_long <- reactive({
    d <- data_all()
    
    # ---- Prezzi ----
    df_p <- as.data.frame(d$P)
    df_p$Date <- as.Date(rownames(df_p))
    
    df_p_long <- tidyr::pivot_longer(
      df_p,
      cols = -Date,
      names_to = "Ticker",
      values_to = "Price"
    )
    
    # ---- Rendimenti ----
    df_r <- as.data.frame(d$R)
    df_r$Date <- as.Date(rownames(df_r))
    
    df_r_long <- tidyr::pivot_longer(
      df_r,
      cols = -Date,
      names_to = "Ticker",
      values_to = "Return"
    )
    
    list(
      prices = df_p_long,
      rets   = df_r_long
    )
  })
  
  
  # ----------------------------------------------------------
  # 4d. OUTPUT TAB 1: ANALISI ESPLORATIVA
  # ----------------------------------------------------------
  
  # ---- Prezzi ----
  output$pricePlot <- renderPlot({
    df <- data_long()$prices
    
    ggplot(df, aes(x = Date, y = Price, color = Ticker)) +
      geom_line(linewidth = 0.7) +
      labs(
        title    = "Prezzi Adjusted",
        subtitle = "Prezzi rettificati per dividendi e split (Yahoo Finance)",
        x = "Data",
        y = "Prezzo (USD)",
        color = "Titoli"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 15)
      )
  })
  
  # ---- Rendimenti ----
  output$retPlot <- renderPlot({
    df <- data_long()$rets
    
    ggplot(df, aes(x = Date, y = Return, color = Ticker)) +
      geom_line(alpha = 0.6, linewidth = 0.4) +
      scale_y_continuous(labels = scales::percent) +
      labs(
        title    = "Rendimenti Giornalieri (Semplici)",
        subtitle = "Variazioni percentuali giornaliere",
        x = "Data",
        y = "Rendimento (%)",
        color = "Titoli"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank()
      )
  })
  
  
  # ---- Correlazioni ----
  output$corrPlot <- renderPlot({
    cor_mat <- cor(data_all()$R)
    
    # Ordinamento gerarchico per evidenziare cluster
    corrplot(
      cor_mat,
      method = "circle",
      type = "upper",
      order = "hclust",
      tl.col = "black",
      diag = FALSE
    )
  })
  
  # ---- Statistiche riassuntive ----
  output$statsTable <- renderDT({
    stats <- data.frame(
      Media_Annua      = colMeans(data_all()$R) * 252,
      Volatilità_Annua = apply(data_all()$R, 2, sd) * sqrt(252)
    )
    
    datatable(round(stats, 4), options = list(pageLength = 15))
  })
  
  
  # ---- Vincoli sui pesi ----
  output$min_limit_info <- renderUI({
    req(input$tickers)
    
    # Parsing ticker
    tickers_list <- unlist(strsplit(gsub(" ", "", input$tickers), ","))
    n <- length(tickers_list)
    
    if (n == 0) return(NULL)
    
    # Minimo teorico per portafoglio long-only equiponderato
    min_theor <- ceiling((1 / n) * 100)
    
    if (input$max_w_slider < min_theor) {
      tags$div(
        style = "color: #d9534f; font-size: 0.85em; margin-top: 5px; font-weight: bold;",
        icon("exclamation-triangle"),
        paste0("Nota: con ", n, " titoli, il minimo reale è ", min_theor, "%")
      )
    } else {
      tags$div(
        style = "color: #666; font-size: 0.85em; margin-top: 5px;",
        paste0("Minimo teorico con ", n, " titoli: ", min_theor, "%")
      )
    }
  })
  

  # ----------------------------------------------------------
  # 4e. OUTPUT TAB 2-7: STRATEGIE DI PORTAFOGLIO
  #
  # observe() viene attivato DOPO il click su "Esegui Analisi".
  # Tutti gli output delle strategie vengono calcolati qui
  # per evitare ridondanze computazionali.
  #
  # Workflow per ogni strategia:
  #   1. Calcolo dei pesi ottimali
  #   2. Metriche portafoglio (ret, vol)
  #   3. Decomposizione del rischio:
  #        - MRC (Marginal Risk Contribution)
  #        - PRC (Percent Risk Contribution)
  #        - Dominance (impatto relativo)
  #   4. Rendering:
  #        - Frontiera + punto
  #        - Barplot pesi
  #        - Scatter rischio
  #        - Tabelle (performance, CAPM, dettaglio)
  # ----------------------------------------------------------
  observe({
    req(input$run)   # Garantisce che l'utente abbia avviato l'analisi
    
    d  <- data_all()
    ef <- markowitz_results()
    
    # ------------------------------------------------------
    # Risk-Free giornaliero
    # Coerente con i rendimenti giornalieri usati nel CAPM
    # ------------------------------------------------------
    rf_daily <- d$rf / 252
    
    
    # ======================================================
    # HELPER FUNCTIONS (riutilizzabili)
    # ======================================================
    
    # ------------------------------------------------------
    # CAPM TABLE
    #
    # Calcola Beta e R² del PORTAFOGLIO (non dei singoli asset)
    #
    # Modello:
    #   (Rp - Rf) = alpha + beta * (Rm - Rf) + errore
    #
    # Nota:
    # - usa RF dinamico (giornaliero)
    # - regressione su serie temporali allineate
    # ------------------------------------------------------
    compute_capm_table <- function(w) {
      
      # Rendimento portafoglio (serie temporale)
      p_ret <- d$R %*% w
      
      # Excess returns (giornalieri)
      p_excess   <- p_ret - as.numeric(d$Rf_vec)
      mkt_excess <- d$Rm  - as.numeric(d$Rf_vec)
      
      # Regressione CAPM
      m_capm <- lm(p_excess ~ mkt_excess)
      s_capm <- summary(m_capm)
      
      data.frame(
        Metrica = c("Beta Portafoglio", "R2"),
        Valore  = c(
          round(coef(m_capm)[2], 3),
          scales::percent(s_capm$r.squared, accuracy = 0.1)
        )
      )
    }
    
    
    # ------------------------------------------------------
    # SUMMARY TABLE
    #
    # Converte metriche giornaliere -> annuali:
    #   Ret_ann = ret * 252
    #   Vol_ann = vol * sqrt(252)
    #
    # Sharpe Ratio:
    #   (Ret_ann - rf_ann) / Vol_ann
    # ------------------------------------------------------
    compute_summary_table <- function(ret, vol) {
      
      # Annualizzazione
      ret_ann <- ret * 252
      vol_ann <- vol * sqrt(252)
      
      # Sharpe Ratio (RF già annuale)
      sharpe_ann <- (ret_ann - d$rf) / vol_ann
      
      data.frame(
        Metrica = c("Rendimento Atteso (Ann.)", "Volatilità (Ann.)", "Sharpe Ratio"),
        Valore  = c(
          scales::percent(ret_ann, accuracy = 0.01),
          scales::percent(vol_ann, accuracy = 0.01),
          round(sharpe_ann, 3)
        )
      )
    }
    
    
    # ------------------------------------------------------
    # DETAIL TABLE
    #
    # Tabella completa per asset:
    #   - Peso
    #   - Contributi al rischio (MRC, PRC)
    #   - Dominanza
    #   - Beta e R² (CAPM asset-level)
    #
    # Nota:
    # PRC somma a 1 → interpretazione percentuale del rischio totale
    # ------------------------------------------------------
    render_detail_table <- function(w, risk_decomp) {
      
      df <- data.frame(
        Ticker    = names(w),
        Peso      = as.numeric(w),
        MRC       = risk_decomp$MRC,
        PRC       = risk_decomp$PRC,
        Dominance = risk_decomp$Dominance,
        Beta      = d$stock_beta,
        R2        = d$stock_rsq
      )
      
      datatable(
        df,
        rownames = FALSE,
        options = list(
          dom = 't',
          pageLength = -1,
          columnDefs = list(list(className = 'dt-center', targets = "_all"))
        )
      ) %>%
        formatPercentage(c('Peso', 'PRC', 'R2'), 2) %>%
        formatRound(c('MRC', 'Dominance', 'Beta'), 4)
    }
    
    
    # ======================================================
    # STRATEGIA 1: MINIMA VARIANZA (MV)
    #
    # Obiettivo:
    #   minimizzare w' Σ w
    #
    # Caratteristiche:
    #   - nessun input su rendimenti attesi
    #   - fortemente guidata da covarianze
    # ======================================================
    res_mv   <- get_min_var(d$mu, d$Sigma, ef = ef)
    risk_mv  <- get_risk_decomposition(res_mv$w, d$Sigma)
    
    output$mvSummary  <- renderTable(compute_summary_table(res_mv$ret, res_mv$vol),
                                     colnames = FALSE, striped = TRUE)
    output$mvCapm     <- renderTable(compute_capm_table(res_mv$w),
                                     colnames = FALSE, striped = TRUE)
    output$mvFrontier <- renderPlot(
      plot_frontier_point(ef$df, res_mv$ret, res_mv$vol)
    )
    output$mvBar      <- renderPlot(plot_weights_bar(res_mv$w))
    output$mvScatter  <- renderPlot(
      plot_risk_scatter(res_mv$w, risk_mv$PRC,
                        "Rischio vs Peso: Minima Varianza")
    )
    output$mvTable    <- renderDT(render_detail_table(res_mv$w, risk_mv))
    
    
    # ======================================================
    # STRATEGIA 2: TANGENCY (MAX SHARPE)
    #
    # Obiettivo:
    #   massimizzare (μ - rf) / σ
    #
    # Nota:
    #   dipende esplicitamente dal risk-free
    # ======================================================
    res_tan  <- get_tangency(d$mu, d$Sigma, rf = rf_daily)
    risk_tan <- get_risk_decomposition(res_tan$w, d$Sigma)
    
    output$tanSummary  <- renderTable(compute_summary_table(res_tan$ret, res_tan$vol),
                                      colnames = FALSE, striped = TRUE)
    output$tanCapm     <- renderTable(compute_capm_table(res_tan$w),
                                      colnames = FALSE, striped = TRUE)
    output$tanFrontier <- renderPlot(
      plot_frontier_point(ef$df, res_tan$ret, res_tan$vol)
    )
    output$tanBar      <- renderPlot(plot_weights_bar(res_tan$w))
    output$tanScatter  <- renderPlot(
      plot_risk_scatter(res_tan$w, risk_tan$PRC,
                        "Rischio vs Peso: Tangency")
    )
    output$tanTable    <- renderDT(render_detail_table(res_tan$w, risk_tan))
    
    
    # ======================================================
    # STRATEGIA 3: EQUAL RISK CONTRIBUTION (ERC)
    #
    # Obiettivo:
    #   PRC_i = 1 / n
    #
    # Interpretazione:
    #   ogni asset contribuisce in modo uguale al rischio totale
    # ======================================================
    res_erc  <- get_erc(d$mu, d$Sigma)
    risk_erc <- get_risk_decomposition(res_erc$w, d$Sigma)
    
    output$ercSummary  <- renderTable(compute_summary_table(res_erc$ret, res_erc$vol),
                                      colnames = FALSE, striped = TRUE)
    output$ercCapm     <- renderTable(compute_capm_table(res_erc$w),
                                      colnames = FALSE, striped = TRUE)
    output$ercFrontier <- renderPlot(
      plot_frontier_point(ef$df, res_erc$ret, res_erc$vol)
    )
    output$ercBar      <- renderPlot(plot_weights_bar(res_erc$w))
    output$ercScatter  <- renderPlot(
      plot_risk_scatter(res_erc$w, risk_erc$PRC,
                        "Rischio vs Peso: ERC")
    )
    output$ercTable    <- renderDT(render_detail_table(res_erc$w, risk_erc))
    
    
    # ======================================================
    # STRATEGIA 4: VINCOLO DI CONCENTRAZIONE
    #
    # Vincolo:
    #   w_i ≤ max_w_slider
    #
    # Uso:
    #   simulazione vincoli regolamentari / UCITS / policy
    # ======================================================
    res_conc  <- get_concentrated(d$mu, d$Sigma, input$max_w_slider)
    risk_conc <- get_risk_decomposition(res_conc$w, d$Sigma)
    
    output$concSummary  <- renderTable(compute_summary_table(res_conc$ret, res_conc$vol),
                                       colnames = FALSE, striped = TRUE)
    output$concCapm     <- renderTable(compute_capm_table(res_conc$w),
                                       colnames = FALSE, striped = TRUE)
    output$concFrontier <- renderPlot(
      plot_frontier_point(ef$df, res_conc$ret, res_conc$vol)
    )
    output$concBar      <- renderPlot(plot_weights_bar(res_conc$w))
    output$concScatter  <- renderPlot(
      plot_risk_scatter(res_conc$w, risk_conc$PRC,
                        "Rischio vs Peso: Concentrazione")
    )
    output$concTable    <- renderDT(render_detail_table(res_conc$w, risk_conc))
    
    
    # ======================================================
    # STRATEGIA 5: TARGET VOLATILITY
    #
    # Logica:
    #   seleziona il punto della frontiera con σ più vicino
    #   al target impostato dall'utente
    #
    # Nota:
    #   conversione Annua -> Giornaliera necessaria
    # ======================================================
    target_val_daily <- (input$target_vol / 100) / sqrt(252)
    
    idx_tv <- which.min(abs(ef$df$vol - target_val_daily))
    
    tv_w   <- ef$weights[idx_tv, ]
    tv_ret <- ef$df$ret[idx_tv]
    tv_vol <- ef$df$vol[idx_tv]
    
    risk_tv <- get_risk_decomposition(tv_w, d$Sigma)
    
    output$tvSummary  <- renderTable(compute_summary_table(tv_ret, tv_vol),
                                     colnames = FALSE, striped = TRUE)
    output$tvCapm     <- renderTable(compute_capm_table(tv_w),
                                     colnames = FALSE, striped = TRUE)
    output$tvFrontier <- renderPlot(
      plot_frontier_point(ef$df, tv_ret, tv_vol)
    )
    output$tvBar      <- renderPlot(plot_weights_bar(tv_w))
    output$tvScatter  <- renderPlot(
      plot_risk_scatter(tv_w, risk_tv$PRC,
                        "Rischio vs Peso: Target Vol")
    )
    output$tvTable    <- renderDT(render_detail_table(tv_w, risk_tv))
    
    
    # ======================================================
    # STRATEGIA 6: EQUAL WEIGHT (Benchmark)
    #
    # w_i = 1/n
    #
    # Ruolo:
    #   benchmark naive per confronto
    # ======================================================
    res_ew  <- get_equal_weight(d$mu, d$Sigma)
    risk_ew <- get_risk_decomposition(res_ew$w, d$Sigma)
    
    output$ewSummary  <- renderTable(compute_summary_table(res_ew$ret, res_ew$vol),
                                     colnames = FALSE, striped = TRUE)
    output$ewCapm     <- renderTable(compute_capm_table(res_ew$w),
                                     colnames = FALSE, striped = TRUE)
    output$ewFrontier <- renderPlot(
      plot_frontier_point(ef$df, res_ew$ret, res_ew$vol)
    )
    output$ewBar      <- renderPlot(plot_weights_bar(res_ew$w))
    output$ewScatter  <- renderPlot(
      plot_risk_scatter(res_ew$w, risk_ew$PRC,
                        "Rischio vs Peso: Equal Weight")
    )
    output$ewTable    <- renderDT(render_detail_table(res_ew$w, risk_ew))
    
  })
  
}



# ============================================================
# SEZIONE 5: AVVIO APPLICAZIONE
#
# Punto di ingresso dell'applicazione Shiny.
#
# shinyApp():
#   - collega l'interfaccia utente (ui)
#   - con la logica server (server)
#
# Flusso:
#   1. L'utente interagisce con la UI
#   2. Gli input vengono inviati al server
#   3. Il server elabora (reactive / observe)
#   4. Gli output vengono renderizzati nella UI
#
# Nota:
# In ambiente di sviluppo (RStudio), questa chiamata:
#   - avvia l'app in una nuova finestra o browser
#   - mantiene attivo il loop reattivo finché l'app è aperta
#
# In produzione (es. shinyapps.io / Shiny Server):
#   - questa istruzione rappresenta il punto di deploy
# ============================================================
shinyApp(ui = ui, server = server)

