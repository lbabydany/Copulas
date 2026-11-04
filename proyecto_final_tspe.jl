# ==============================================================================
# PROYECTO FINAL: Temas Selectos de Probabilidad y Estadística
# Modelación de Dependencia Oceanográfica: Velocidad del Viento vs Altura del Oleaje
# Autor: Sánchez Benítez Luis Daniel - 320088175
# ==============================================================================

using CSV
using DataFrames
using Plots
using Statistics
using Distributions
using Copulas
using Downloads: download

# Inclusión de librerías y funciones auxiliares locales
include("copulasaux.jl")
include("copulas.jl")

println(">>> Iniciando ejecución del proyecto final...")

# ------------------------------------------------------------------------------
# 1. CARGA Y LIMPIEZA DE DATOS
# ------------------------------------------------------------------------------
url_bd = "https://raw.githubusercontent.com/lbabydany/Copulas/refs/heads/main/320088175.csv"
ruta_local_bd = "320088175.csv"

# Descarga remota si no existe localmente
if !isfile(ruta_local_bd)
    println("Descargando base de datos desde GitHub...")
    download(url_bd, ruta_local_bd)
end

datos = CSV.read(ruta_local_bd, DataFrame)

# Limpieza: La NOAA utiliza el valor centinela 99.0 para indicar datos faltantes
datos_limpios = datos[(datos[:, 1] .!= 99.0) .& (datos[:, 2] .!= 99.0), :]

viento = datos_limpios[:, 1] # Velocidad del viento (m/s)
oleaje = datos_limpios[:, 2] # Altura del oleaje (m)

println("Total de observaciones válidas tras filtrado: ", nrow(datos_limpios))

# ------------------------------------------------------------------------------
# 2. ANÁLISIS EXPLORATORIO UNIVARIADO
# ------------------------------------------------------------------------------
function calcular_estadisticas(variable, nombre)
    println("\n--- Estadísticas descriptivas: $nombre ---")
    println("Mínimo:      ", round(minimum(variable), digits=4))
    println("Cuartil 1:   ", round(quantile(variable, 0.25), digits=4))
    println("Mediana:     ", round(median(variable), digits=4))
    println("Media:       ", round(mean(variable), digits=4))
    println("Cuartil 3:   ", round(quantile(variable, 0.75), digits=4))
    println("Máximo:      ", round(maximum(variable), digits=4))
    println("------------------------------------------")
end

calcular_estadisticas(viento, "Velocidad del Viento (m/s)")
calcular_estadisticas(oleaje, "Altura del Oleaje (m)")

# Histogramas marginales exploratorios
p_hist_viento = histogram(
    viento,
    title = "Velocidad del Viento",
    xlabel = "Velocidad (m/s)",
    ylabel = "Frecuencia",
    legend = false,
    color = :skyblue,
    bins = :scott
)

p_hist_oleaje = histogram(
    oleaje,
    title = "Altura del Oleaje",
    xlabel = "Altura (m)",
    ylabel = "Frecuencia",
    legend = false,
    color = :lightgreen,
    bins = :scott
)

p_dispersion = scatter(
    viento, oleaje,
    title = "Dispersión: Viento vs Oleaje",
    xlabel = "Velocidad del Viento (m/s)",
    ylabel = "Altura del Oleaje (m)",
    legend = false,
    color = :purple,
    markersize = 2.5,
    markeralpha = 0.45,
    markerstrokewidth = 0
)

# ------------------------------------------------------------------------------
# 3. CONSTRUCCIÓN DE LA CÓPULA EMPÍRICA Y MEDIDAS DE ASOCIACIÓN
# ------------------------------------------------------------------------------
println("\n>>> Evaluando medidas de concordancia y dependencia...")
cobs = copulaPseudObs(viento, oleaje)
u = cobs[:, 1]
v = cobs[:, 2]

m_res = 100
sw_val = Schweizer(cobs, m_res)
sp_val = Spearman(cobs, m_res)

println("Dependencia de Schweizer-Wolff (m = $m_res): ", round(sw_val, digits=4))
println("Rho de Spearman empírico       (m = $m_res): ", round(sp_val, digits=4))

p_pseudo = scatter(
    u, v,
    title = "Pseudo-observaciones (Cópula Empírica)",
    xlabel = "U",
    ylabel = "V",
    legend = false,
    color = :darkblue,
    markersize = 2.5,
    markeralpha = 0.45,
    markerstrokewidth = 0,
    size = (600, 600)
)

# Diagonales de Fréchet-Hoeffding para la cópula empírica
t = range(0, 1, length=100)
diag_princ = [Cn(ti, ti, cobs) for ti in t]
diag_sec   = [Cn(ti, 1 - ti, cobs) for ti in t]

p_diag_princ = plot(
    t, diag_princ,
    title = "Diagonal Principal C_n(t, t)",
    xlabel = "t", ylabel = "C_n(t, t)",
    label = "Cópula Empírica",
    color = :blue, lw = 2,
    legend = :topleft
)
plot!(p_diag_princ, t, t, label = "Cota Sup (M)", ls = :dash, color = :black)
plot!(p_diag_princ, t, max.(2 .* t .- 1, 0), label = "Cota Inf (W)", ls = :dot, color = :red)

p_diag_sec = plot(
    t, diag_sec,
    title = "Diagonal Secundaria C_n(t, 1-t)",
    xlabel = "t", ylabel = "C_n(t, 1-t)",
    label = "Cópula Empírica",
    color = :purple, lw = 2,
    legend = :topright
)
plot!(p_diag_sec, t, min.(t, 1 .- t), label = "Cota Sup (M)", ls = :dash, color = :black)
plot!(p_diag_sec, t, zeros(length(t)), label = "Cota Inf (W)", ls = :dot, color = :red)

# ------------------------------------------------------------------------------
# 4. AJUSTE PARAMÉTRICO DE DISTRIBUCIONES MARGINALES (MLE)
# ------------------------------------------------------------------------------
println("\n>>> Ajustando marginales univariadas por MLE...")
viento_pos = viento[viento .> 0]
oleaje_pos = oleaje[oleaje .> 0]

# 4.1 Viento ~ Gamma(α, θ)
gamma_viento_fit = fit_mle(Gamma, viento_pos)
α_gamma_viento   = shape(gamma_viento_fit)
θ_gamma_viento   = scale(gamma_viento_fit)

# 4.2 Viento ~ Weibull(k, λ)
weibull_viento_fit = fit_mle(Weibull, viento_pos)
k_weibull_viento   = shape(weibull_viento_fit)
λ_weibull_viento   = scale(weibull_viento_fit)

# 4.3 Oleaje ~ Weibull(k, λ)
weibull_oleaje_fit = fit_mle(Weibull, oleaje_pos)
k_weibull_oleaje   = shape(weibull_oleaje_fit)
λ_weibull_oleaje   = scale(weibull_oleaje_fit)

println("Viento Gamma:   α = ", round(α_gamma_viento, digits=4), " | θ = ", round(θ_gamma_viento, digits=4))
println("Viento Weibull: k = ", round(k_weibull_viento, digits=4), " | λ = ", round(λ_weibull_viento, digits=4))
println("Oleaje Weibull: k = ", round(k_weibull_oleaje, digits=4), " | λ = ", round(λ_weibull_oleaje, digits=4))

# Gráficas de densidades ajustadas
x_grid_viento = range(minimum(viento_pos), maximum(viento_pos), length=500)
x_grid_oleaje = range(minimum(oleaje_pos), maximum(oleaje_pos), length=500)

p_fit_viento_gamma = histogram(viento_pos, normalize=true, bins=:scott, alpha=0.4, color=:skyblue, label="Datos")
plot!(p_fit_viento_gamma, x_grid_viento, pdf.(gamma_viento_fit, x_grid_viento), lw=2.5, color=:blue, label="Gamma", title="Viento - Gamma")

p_fit_viento_weibull = histogram(viento_pos, normalize=true, bins=:scott, alpha=0.4, color=:skyblue, label="Datos")
plot!(p_fit_viento_weibull, x_grid_viento, pdf.(weibull_viento_fit, x_grid_viento), lw=2.5, color=:darkblue, label="Weibull", title="Viento - Weibull")

p_fit_oleaje_weibull = histogram(oleaje_pos, normalize=true, bins=:scott, alpha=0.4, color=:lightgreen, label="Datos")
plot!(p_fit_oleaje_weibull, x_grid_oleaje, pdf.(weibull_oleaje_fit, x_grid_oleaje), lw=2.5, color=:darkgreen, label="Weibull", title="Oleaje - Weibull")

p_marginales_comp = plot(p_fit_viento_gamma, p_fit_viento_weibull, p_fit_oleaje_weibull, layout=(1, 3), size=(1300, 380))

# ------------------------------------------------------------------------------
# 5. CALIBRACIÓN DE CÓPULA DE GUMBEL POR ABC (APPROXIMATE BAYESIAN COMPUTATION)
# ------------------------------------------------------------------------------
println("\n>>> Ejecutando calibración ABC para la Cópula de Gumbel...")

function métodoABC_Gumbel(m, k, cobs::Matrix{<:Real})
    n = size(cobs, 1)
    CnObs = zeros(n)
    for i in 1:n
        CnObs[i] = Cn(cobs[i, 1], cobs[i, 2], cobs)
    end

    θsim = rand(Uniform(2, 5), m)
    ε = zeros(m)

    for j in 1:m
        csim = Matrix(rand(GumbelCopula(2, θsim[j]), n)')
        for i in 1:n
            ε[j] += (CnObs[i] - Cn(cobs[i, 1], cobs[i, 2], csim))^2
        end
    end

    jopt = sortperm(ε)[1:k]
    return θsim[jopt], ε[jopt]
end

# Ejecución del algoritmo ABC
@time θabc, εabc = métodoABC_Gumbel(200, 20, cobs)
θ_hat = mean(θabc)
println("Parámetro θ estimado (Gumbel): ", round(θ_hat, digits=4))

# Cópula paramétrica calibrada
cop_fit = GumbelCopula(2, θ_hat)

# Simulación de verificación
n_sim = 15000
cobs_sim = Matrix(rand(cop_fit, n_sim)')

println("Rho de Spearman en datos simulados: ", round(Spearman(cobs_sim, 100), digits=4))
println("Schweizer-Wolff en datos simulados:  ", round(Schweizer(cobs_sim, 100), digits=4))

p_sim = scatter(
    cobs_sim[:, 1], cobs_sim[:, 2],
    title = "Simulación Gumbel (θ = $(round(θ_hat, digits=3)))",
    xlabel = "U simulada", ylabel = "V simulada",
    legend = false,
    color = :darkred,
    markersize = 2,
    markeralpha = 0.35,
    markerstrokewidth = 0,
    size = (600, 600)
)

# ------------------------------------------------------------------------------
# 6. COMPARACIÓN DE CURVAS DE NIVEL (EMPÍRICA VS GUMBEL AJUSTADA)
# ------------------------------------------------------------------------------
println("\n>>> Evaluando curvas de nivel...")
Z_emp = [Cn(i, j, cobs) for i in t, j in t]
Z_fit = [cdf(cop_fit, [i, j]) for i in t, j in t]

p_emp_nivel = contour(
    t, t, Z_emp,
    title = "Curvas de Nivel\nCópula Empírica",
    xlabel = "u", ylabel = "v",
    fill = true, levels = 20,
    color = :viridis,
    colorbar_title = "C_n(u,v)"
)

p_fit_nivel = contour(
    t, t, Z_fit,
    title = "Curvas de Nivel\nGumbel Ajustada",
    xlabel = "u", ylabel = "v",
    fill = true, levels = 20,
    color = :plasma,
    colorbar_title = "C_θ(u,v)"
)

p_contour_comp = plot(p_emp_nivel, p_fit_nivel, layout=(1, 2), size=(1100, 420))

# ------------------------------------------------------------------------------
# 7. GUARDADO AUTOMÁTICO DE FIGURAS PARA EL REPOSITORIO
# ------------------------------------------------------------------------------
println("\n>>> Guardando figuras en disco para el README...")
savefig(p_dispersion, "dispersion_viento_oleaje.png")
savefig(p_pseudo, "pseudo_observaciones.png")
savefig(p_marginales_comp, "ajuste_marginales.png")
savefig(p_sim, "simulacion_gumbel.png")
savefig(p_contour_comp, "curvas_nivel_comparacion.png")

println(">>> Proceso finalizado con éxito.")