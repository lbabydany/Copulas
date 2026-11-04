# Modelación de Dependencia Oceanográfica Viento-Oleaje mediante Cópula de Gumbel

[![Julia](https://img.shields.io/badge/Language-Julia_1.9+-9558B2?style=flat&logo=julia)](https://julialang.org)
[![Data](https://img.shields.io/badge/Data-NOAA_NDBC-003366.svg)](https://www.ndbc.noaa.gov/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Proyecto final desarrollado para la asignatura **Temas Selectos de Probabilidad y Estadística** (Licenciatura en Actuaría, UNAM FES Acatlán).

---

## 📌 Contexto Físico y Justificación Teórica

En la oceanografía física y el análisis actuarial de riesgos climáticos, la relación entre la **Velocidad del Viento ($m/s$)** y la **Altura Significativa del Oleaje ($m$)** representa un sistema acoplado complejo:
- **Mecanismo generador:** El viento transfiere energía cinética directamente a la superficie marina para originar el oleaje.
- **Asimetría temporal (*Mar de fondo*):** Vientos extremos generan incrementos inmediatos en el oleaje, pero alturas considerables persisten aun cuando el viento ha calmado.
- **Eventos extremos simultáneos:** Asumir dependencia lineal o normalidad subestima severamente las colas conjuntas. Mediante el **Teorema de Sklar**, este proyecto separa el ajuste de las marginales de la estructura de dependencia asimétrica calibrada a través de la **Cópula de Gumbel**.

---

## 🔬 Flujo Metodológico

1. **Limpieza de Datos:** Tratamiento y filtrado de registros horarios de la boya oceanográfica de la NOAA (National Data Buoy Center), descartando valores centinela nulos (`99.0`).
2. **Cópula Empírica y Pseudo-Observaciones:**
   - Transformación de observaciones continuas al hipercubo unitario $[0, 1]^2$ mediante rangos empíricos ($U_i, V_i$).
   - Cuantificación de asociación no paramétrica: $\rho_s$ de Spearman y medida de Schweizer-Wolff.
   - Validación de fronteras respecto a las Cotas de Fréchet-Hoeffding en las diagonales $C_n(t, t)$ y $C_n(t, 1-t)$.
3. **Ajuste Paramétrico de Marginales (MLE):**
   - **Viento:** Comparación de bondad de ajuste entre $\text{Gamma}(\alpha, \theta)$ y $\text{Weibull}(k, \lambda)$ por Máxima Verosimilitud.
   - **Oleaje:** Calibración de distribución $\text{Weibull}(k, \lambda)$ por Máxima Verosimilitud.
4. **Calibración de la Cópula de Gumbel por ABC:**
   - Debido a la concentración asimétrica en el cuadrante $(1, 1)$ (dependencia de cola superior $\lambda_U > 0$), se calibra el parámetro $\theta$ de Gumbel mediante un algoritmo de **Cómputo Bayesiano Aproximado (Approximate Bayesian Computation - ABC)** minimizando el error cuadrático respecto a la cópula empírica.

---

## 📂 Estructura del Repositorio

```text
├── 320088175.csv                   # Base de datos histórica procesada de la boya NOAA
├── copulas.jl                      # Implementación de C_n(u,v), Spearman y Schweizer-Wolff
├── copulasaux.jl                   # Cálculo de rangos y pseudo-observaciones empíricas
├── proyecto_final_tspe.jl               # Pipeline completo (ETL, MLE, calibración ABC y visualización)
├── dispersion_viento_oleaje.png    # Gráfica: Dispersión en unidades físicas
├── pseudo_observaciones.png        # Gráfica: Distribución en el cuadrado unitario [0, 1]²
├── ajuste_marginales.png           # Gráfica: Comparativa de marginales ajustadas
├── simulacion_gumbel.png           # Gráfica: Simulación bivariada de la cópula ajustada
├── curvas_nivel_comparacion.png    # Gráfica: Contornos Cópula Empírica vs Gumbel
└── README.md
