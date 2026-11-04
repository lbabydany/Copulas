# Modelación de Dependencia Oceanográfica Viento-Oleaje mediante Cópula de Gumbel

[![Julia](https://img.shields.io/badge/Language-Julia_1.9+-9558B2?style=flat&logo=julia)](https://julialang.org)
[!(https://img.shields.io/badge/Data-NOAA_NDBC-003366.svg)](https://www.ndbc.noaa.gov/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Proyecto final desarrollado para la asignatura **Temas Selectos de Probabilidad y Estadística** (Licenciatura en Actuaría, UNAM FES Acatlán).

---

## 📌 Contexto Físico y Relevancia Actuarial

En oceanografía física y gestión de riesgos hidrometeorológicos, modelar conjuntamente la **Velocidad del Viento ($m/s$)** y la **Altura Significativa del Oleaje ($m$)** es crítico para:
- Cuantificar riesgos de navegación e infraestructuras offshore/costeras.
- Evaluar probabilidades de falla estructural ante eventos extremos conjuntos (huracanes y marejadas).
- Modelar adecuadamente fenómenos asimétricos: el viento transfiere energía cinética al mar de forma rápida, pero el oleaje elevado persiste aun cuando el viento amaina (*mar de fondo*).

Asumir dependencia lineal o normalidad subestima severamente la cola superior. Por ello, se implementa el **Teorema de Sklar** para aislar la estructura de dependencia asimétrica mediante la **Cópula de Gumbel** tras modelar las distribuciones marginales.

---

## 🔬 Flujo de Trabajo Metodológico

```text
  Registros Crudos NOAA ──> Filtrado de Nulos (99.0) ──> Pseudo-observaciones (U, V)
                                                                    │
           ┌────────────────────────────────────────────────────────┴───────────────────────────┐
           ▼                                                                                    ▼
  Ajuste de Marginales                                                             Cópula Empírica vs Paramétrica
  ├─ Viento: Gamma(α, θ) vs Weibull(k, λ) [MLE]                                     ├─ Schweizer-Wolff y Rho de Spearman
  └─ Oleaje: Weibull(k, λ) [MLE]                                                    ├─ Diagonales Principal C(t,t) y Secundaria C(t, 1-t)
                                                                                    └─ Calibración ABC de Cópula de Gumbel: θ̂
---

## 🔬 Marco Metodológico

1. **Ajuste de Marginales ($F_X, F_Y$):**
   - Transformación de variables continuas al hipercubo unitario $[0, 1]^d$ mediante la función de distribución empírica (pseudo-observaciones) y distribuciones paramétricas seleccionadas.
   - Pruebas de bondad de ajuste (Kolmogorov-Smirnov / Anderson-Darling).

2. **Calibración de Cópulas:**
   - Estimación de parámetros por **Máxima Verosimilitud (MLE)** e inversión de medidas de asociación no paramétricas ($\tau$ de Kendall, $\rho_s$ de Spearman).
   - Familias evaluadas:
     - **Gaussiana:** Dependencia lineal simétrica sin dependencia de cola ($\lambda_U = \lambda_L = 0$).
     - **t-Student:** Dependencia simétrica en colas pesadas ($\lambda_U = \lambda_L > 0$).
     - **Clayton:** Dependencia fuerte en la cola inferior ($\lambda_L > 0, \lambda_U = 0$).
     - **Gumbel:** Dependencia fuerte en la cola superior ($\lambda_U > 0, \lambda_L = 0$).
     - **Frank:** Dependencia simétrica sin dependencia de colas extremas.

3. **Criterios de Selección y Diagnóstico:**
   - Comparación mediante **AIC** (Akaike Information Criterion) y **BIC** (Bayesian Information Criterion).
   - Inspección gráfica de contornos empíricos vs. teóricos y gráficos K-plot / Chi-plot.

---

## 📊 Datos Utilizados

- **Fuente:** *Registros históricos de boyas del National Data Buoy Center (NDBC - NOAA).
