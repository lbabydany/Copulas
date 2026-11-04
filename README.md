# Modelación de Dependencia Multivariada Mediante Cópulas

[![Julia](https://img.shields.io/badge/Language-Julia_1.10+-9558B2?style=flat&logo=julia)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Status](https://img.shields.io/badge/Status-Completed-success)](#)

> Proyecto final para la materia **Temas Selectos de Probabilidad y Estadística** (Licenciatura en Actuaría, UNAM FES Acatlán).

---

## 📌 Resumen Ejecutivo

En el análisis cuantitativo de riesgos y fenómenos multivariados, asumir normalidad multivariada o dependencia lineal (correlación de Pearson) subestima gravemente los eventos de cola y la dependencia asimétrica. 

Este proyecto implementa la teoría del **Teorema de Sklar** para desacoplar el comportamiento de las distribuciones marginales de la estructura de dependencia conjunta. A través de familias paramétricas de **cópulas elípticas** (Gaussiana, t-Student) y **arquimedianas** (Clayton, Gumbel, Frank), se modela y evalúa la estructura de riesgo conjunto en datos empíricos.

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

- **Fuente:** *(Describe brevemente la base de datos: e.g., series
