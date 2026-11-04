# Modelación de Dependencia Oceanográfica Viento-Oleaje mediante Cópula de Gumbel

[![Julia](https://img.shields.io/badge/Language-Julia_1.9+-9558B2?style=flat&logo=julia)](https://julialang.org)
[![Data](https://img.shields.io/badge/Data-NOAA_NDBC-003366.svg)](https://www.ndbc.noaa.gov/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Proyecto final desarrollado para la asignatura **Temas Selectos de Probabilidad y Estadística** (Licenciatura en Actuaría, UNAM FES Acatlán).

---

## 📌 Contexto Físico y Justificación Teórica

En la oceanografía física y el análisis actuarial de riesgos climáticos, la relación entre la **Velocidad del Viento ($m/s$)** y la **Altura Significativa del Oleaje ($m$)** representa un sistema acoplado complejo:
- **Mecanismo generador:** El viento transfiere energía cinética directamente a la superficie marina para originar el oleaje[cite: 1].
- **Asimetría temporal (*Mar de fondo*):** Vientos extremos generan incrementos inmediatos en el oleaje, pero alturas considerables persisten aun cuando el viento ha calmado[cite: 1].
- **Eventos extremos simultáneos:** Asumir dependencia lineal o normalidad subestima severamente las colas conjuntas[cite: 1]. Mediante el **Teorema de Sklar**, este proyecto separa el ajuste de las marginales de la estructura de dependencia asimétrica calibrada a través de la **Cópula de Gumbel**[cite: 1].

---

## 🔬 Flujo Metodológico

1. **Limpieza de Datos:** Tratamiento y filtrado de registros horarios de la boya oceanográfica de la NOAA (National Data Buoy Center), descartando valores centinela nulos (`99.0`)[cite: 1].
2. **Cópula Empírica y Pseudo-Observaciones:**
   - Transformación de observaciones continuas al hipercubo unitario $[0, 1]^2$ mediante rangos empíricos ($U_i, V_i$)[cite: 1].
   - Cuantificación de asociación no paramétrica: $\rho_s$ de Spearman y medida de Schweizer-Wolff[cite: 1].
   - Validación de fronteras respecto a las Cotas de Fréchet-Hoeffding en las diagonales $C_n(t, t)$ y $C_n(t, 1-t)$[cite: 1].
3. **Ajuste Paramétrico de Marginales (MLE):**
   - **Viento:** Comparación de bondad de ajuste entre $\text{Gamma}(\alpha, \theta)$ y $\text{Weibull}(k, \lambda)$ por Máxima Verosimilitud[cite: 1].
   - **Oleaje:** Calibración de distribución $\text{Weibull}(k, \lambda)$ por Máxima Verosimilitud[cite: 1].
4. **Calibración de la Cópula de Gumbel por ABC:**
   - Debido a la concentración asimétrica en el cuadrante $(1, 1)$ (dependencia de cola superior $\lambda_U > 0$)[cite: 1], se calibra el parámetro $\theta$ de Gumbel mediante un algoritmo de **Cómputo Bayesiano Aproximado (Approximate Bayesian Computation - ABC)** minimizando el error cuadrático respecto a la cópula empírica[cite: 1].

---
