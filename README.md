# Advanced Quantitative Asset Allocation & Risk Budgeting Framework

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.2.0-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-Interactive%20Dashboard-brightgreen.svg)](https://shiny.posit.co/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Dashboard interattiva e framework di ottimizzazione quantitativa di portafoglio sviluppati in **R** e **Shiny**, affiancati da una relazione tecnica approfondita sui modelli di *Modern Portfolio Theory* (MPT), *Risk Parity* e *Capital Asset Pricing Model* (CAPM).

---

## Panoramica del Progetto

Il progetto implementa e confronta sei diverse strategie di allocazione di portafoglio su un universo azionario eterogeneo (Technology, Healthcare, Energy, Financials, Consumer Discretionary, Materials, Utilities, Industrials), integrando l'acquisizione automatizzata di serie storiche finanziarie, l'ottimizzazione vincolata e la scomposizione avanzata del rischio.

I risultati sono analizzati con particolare attenzione a:
- **Concentrazione dei pesi** e diversificazione settoriale.
- **Ripartizione del rischio** (*Marginal Risk Contribution* e *Percentage Risk Contribution*).
- **Sensibilità al rischio sistematico** (regressione CAPM: $\beta$ e $R^2$).
- Valutazione per profili di investimento ad **elevata tolleranza al rischio** e **lungo orizzonte temporale**.

---

## Strategie di Portafoglio Implementate

1. **Portafoglio a Varianza Minima (Global Minimum Variance):**  
   Minimizzazione del rischio totale del portafoglio tramite programmazione quadratica (`solve.QP`).
2. **Portafoglio di Tangenza (Maximum Sharpe Ratio):**  
   Massimizzazione dello Sharpe Ratio tramite trasformazione ausiliaria del vincolo di rendimento in formulazione quadratica.
3. **Equal Risk Contribution (ERC / Risk Parity):**  
   Allocazione bilanciata del rischio basata sulla formulazione log-penalizzata di Roncalli:
   $$\min_{w} \frac{1}{2} w^T \Sigma w - \lambda \sum_{i=1}^{n} \log(w_i)$$
   Risolta numericamente con algoritmo quasi-Newton `L-BFGS-B`.
4. **Portafoglio con Vincolo di Concentrazione:**  
   Controllo del rischio specifico tramite imposizione di un tetto massimo sul peso del singolo asset ($w_i \le w_{max}$).
5. **Portafoglio con Target di Volatilità:**  
   Selezione del punto ottimo lungo la frontiera efficiente corrispondente alla volatilità obiettivo.
6. **Portafoglio Equipesato (Equal Weight / 1/N):**  
   Benchmark euristico per la valutazione dell'efficacia delle ottimizzazioni parametriche.

---

## Stack Tecnologico & Librerie R

- **UI & Interattività:** `shiny`, `DT`
- **Data Ingestion & Time Series:** `quantmod`, `xts`, `PerformanceAnalytics`
- **Ottimizzazione Numerica:** `quadprog`, `stats::optim` (L-BFGS-B)
- **Data Wrangling & Visualizzazione:** `tidyr`, `ggplot2`, `corrplot`, `scales`

---

## Come Eseguire l'Applicazione

### Prerequisiti
Assicurati di avere installato R (versione $\ge 4.2$) e facoltativamente RStudio.

### 1. Clona il Repository
```bash
git clone https://github.com/<TUO-USERNAME>/quantitative-asset-allocation-shiny.git
cd quantitative-asset-allocation-shiny
```

### 2. Installa le Dipendenze
All'interno della console di R:
```R
install.packages(c(
  "shiny", "quantmod", "xts", "PerformanceAnalytics", 
  "quadprog", "DT", "ggplot2", "corrplot", "tidyr", "scales"
))
```

### 3. Avvia la Dashboard
```R
library(shiny)
runApp("app.R")
```

---

## Documentazione e Relazione Tecnica

All'interno della cartella [`docs/`](docs/) è disponibile il report completo in formato PDF:
- **`ProgettazioneStrategieAssetAllocation.pdf`**: Tratta la formulazione matematica dei problemi di ottimizzazione, le matrici di covarianza, l'interpretazione dei contributi marginali al rischio e le conclusioni empiriche sui portafogli analizzati.

---

## Autore
- **Pietro Peluso** 
