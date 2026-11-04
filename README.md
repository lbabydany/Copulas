# Modelación de Dependencia Oceanográfica Viento-Oleaje mediante Cópula de Gumbel

[![Julia](https://img.shields.io/badge/Language-Julia_1.9+-9558B2?style=flat&logo=julia)](https://julialang.org)
[![Data: NOAA NDBC](https://img.shields.io/badge/Data-NOAA_NDBC-003366.svg)](https://www.ndbc.noaa.gov/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Proyecto final desarrollado para la asignatura **Temas Selectos de Probabilidad y Estadística** (Licenciatura en Actuaría, UNAM FES Acatlán).

---

## 📌 Contexto Físico y Relevancia Actuarial

En oceanografía física y gestión cuantitativa de riesgos hidrometeorológicos, modelar conjuntamente la **Velocidad del Viento ($m/s$)** y la **Altura Significativa del Oleaje ($m$)** es fundamental para:
- Cuantificar riesgos de navegación e infraestructuras marítimas/costeras.
- Evaluar probabilidades de falla estructural ante eventos extremos concurrentes (huracanes, vendavales y marejadas).
- Modelar fenómenos dinámicos asimétricos: el viento transfiere energía cinética al mar con rapidez, pero el oleaje elevado persiste aun cuando el viento amaina (*mar de fondo*).

Asumir dependencia lineal o normalidad multivariada subestima severamente los escenarios en colas extremas. Por ello, se aplica el **Teorema de Sklar** para desacoplar el comportamiento de las distribuciones marginales de la estructura de dependencia asimétrica mediante la **Cópula de Gumbel**, calibrada a través de métodos de simulación estocástica.

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
