# PROYECTO FINAL - Temas Selectos de Probabilidad y Estadística

# ■ INTRODUCCIÓN.
# El presente proyecto analiza la relación de dependencia entre dos variables 
# aleatorias continuas fundamentales en la oceanografía física: la Velocidad del 
# Viento (medida en m/s) y la Altura Significativa del Oleaje (medida en metros).
#
# Justificación de la dependencia:
# Se espera una dependencia positiva y robusta entre ambas variables, dado que el 
# viento es el mecanismo físico principal que transfiere energía cinética a la 
# superficie del mar para generar el oleaje. No obstante, se anticipa una 
# estructura de dependencia no lineal y asimétrica; mientras que vientos extremos 
# suelen generar oleajes altos de forma casi inmediata, un oleaje elevado puede 
# persistir incluso cuando el viento ha calmado (fenómeno conocido como mar de fondo).
# Esta complejidad técnica justifica el uso de herramientas avanzadas como las 
# cópulas para caracterizar su dependencia subyacente de manera precisa.
#
# Relevancia y utilidad:
# Analizar esta relación es crítico para la gestión de riesgos en la navegación, 
# el diseño de infraestructura costera y la predicción de inundaciones por 
# tormentas. Modelar adecuadamente las colas de la distribución permite evaluar 
# la probabilidad de ocurrencia de eventos extremos simultáneos.
#
# Origen de los datos:
# Las observaciones se obtuvieron de los registros históricos del National Data 
# Buoy Center (NDBC) de la NOAA (National Oceanic and Atmospheric Administration).
# Los datos corresponden a mediciones horarias procesadas para eliminar valores 
# nulos o erróneos.
# Enlace de origen: https://www.ndbc.noaa.gov/
# ==============================================================================

#Sánchez Benítez Luis Daniel - 320088175

using CSV, DataFrames, Plots, Statistics, Copulas, Distributions 

include("copulasaux.jl")

include("copulas.jl")

url_bd = "https://raw.githubusercontent.com/lbabydany/Copulas/refs/heads/main/320088175.csv"
bd = download(url_bd)
datos = CSV.read(bd,  DataFrame)
# LIMPIEZA DE DATOS: 
# La NOAA utiliza el valor 99.0 para indicar datos faltantes o nulos.
# Filtramos el DataFrame conservando únicamente las filas donde TANTO la 
# columna 1 (viento) COMO la columna 2 (oleaje) sean diferentes de 99.0.
# Usamos el operador de "broadcasting" (.) para evaluar elemento por elemento.
    datos_limpios = datos[(datos[:, 1] .!= 99.0) .& (datos[:, 2] .!= 99.0), :]

# Separamos las dos variables aleatorias continuas usando el DataFrame limpio
viento = datos_limpios[:, 1] # Velocidad del viento (m/s)
oleaje = datos_limpios[:, 2] # Altura del oleaje (m)

println("Total de observaciones válidas tras la limpieza: ", nrow(datos_limpios))

# ■ ANÁLISIS EXPLORATORIO
# • Mínimo, máximo, cuartiles y media muestrales, para cada variable
# Creamos una función auxiliar para no repetir código y mantenerlo limpio
function calcular_estadisticas(variable, nombre)
    println("--- Estadísticas para: $nombre ---")
    println("Mínimo:      ", minimum(variable))
    println("Cuartil 1:   ", quantile(variable, 0.25))
    println("Mediana:     ", median(variable)) # Equivale al Cuartil 2
    println("Media:       ", mean(variable))
    println("Cuartil 3:   ", quantile(variable, 0.75))
    println("Máximo:      ", maximum(variable))
    println("----------------------------------\n")
end

# Ejecutamos la función para ambas variables
calcular_estadisticas(viento, "Velocidad del Viento (m/s)")
calcular_estadisticas(oleaje, "Altura del Oleaje (m)")


# ==============================================================================
# ■ MEDIDAS DE DEPENDENCIA (CÓPULA EMPÍRICA)
# ==============================================================================
# Para evaluar la relación asimétrica y no lineal entre el viento y el oleaje 
# sin depender de sus distribuciones marginales, construimos en este punto las 
# pseudo-observaciones empíricas.
cobs = copulaPseudObs(viento, oleaje)

# • Medida de dependencia de Schweizer-Wolff
# Captura el nivel general de dependencia (incluso si fuese no monótona) en la 
# transferencia de energía del viento hacia la superficie del mar.
m = 100
sw_val = Schweizer(cobs, m)
println("Dependencia de Schweizer-Wolff (m=$m): ", round(sw_val, digits=4))

# • Medida de concordancia Rho de Spearman
# Evalúa la monotonía de la relación. Esperamos un valor positivo alto, lo que 
# confirmará la hipótesis: a mayor velocidad de viento, mayor altura de oleaje.
sp_val = Spearman(cobs, m)
println("Rho de Spearman empírico (m=$m):       ", round(sp_val, digits=4))
println("----------------------------------\n")


# • Histogramas marginales

p_hist_viento = histogram(viento, 
    title = "Histograma Marginal: Velocidad del Viento",
    xlabel = "Velocidad (m/s)", 
    ylabel = "Frecuencia",
    legend = false,
    color = :skyblue,
    bins = :scott
)

# Generamos el histograma para la variable 2 (Altura del oleaje).
p_hist_oleaje = histogram(oleaje, 
    title = "Histograma Marginal: Altura del Oleaje",
    xlabel = "Altura (m)", 
    ylabel = "Frecuencia",
    legend = false,
    color = :lightgreen,
    bins = :scott
)

# ------------------------------------------------------------------------------
# • Diagrama de dispersión de las variables
# ------------------------------------------------------------------------------
# Graficamos la relación directa entre el viento (eje X) y el oleaje (eje Y).
# Utilizamos 'markeralpha = 0.5' para darle transparencia a los puntos; 
# esto es súper útil porque tendrás miles de observaciones de la boya de la NOAA,
# y la transparencia te permitirá visualizar en dónde hay mayor densidad de datos.

p_dispersion = scatter(viento, oleaje,
    title = "Diagrama de Dispersión: Viento vs Oleaje",
    xlabel = "Velocidad del Viento (m/s)",
    ylabel = "Altura del Oleaje (m)",
    legend = false,
    color = :purple,
    markersize = 3,
    markeralpha = 0.5,
    markerstrokewidth = 0 # Quitar el borde de los puntos para que se vea más limpio
)

# Construcción de pseudo-observaciones
cobs = copulaPseudObs(viento, oleaje)

# Extraemos las coordenadas U y V
u = cobs[:, 1]
v = cobs[:, 2]

# Visualización de las pseudo-observaciones
p_pseudo = scatter(u, v,
    title = "Pseudo-observaciones de la Cópula Empírica",
    xlabel = "U",
    ylabel = "V",
    legend = false,
    color = :darkblue,
    markersize = 2.5,
    markeralpha = 0.45,
    markerstrokewidth = 0, 
    size=(600,600)
)

# ==============================================================================
# ■ ANÁLISIS DE LA CÓPULA EMPÍRICA (Cn)
# ==============================================================================

# Definimos una cuadrícula (grid) en el intervalo [0, 1] para evaluar la cópula
t = range(0, 1, length=100)

# ------------------------------------------------------------------------------
# 1. Conjuntos de Nivel (Curvas de nivel / Contour)
# ------------------------------------------------------------------------------
# Evaluamos la cópula empírica en cada combinación de la cuadrícula
Z = [Cn(i, j, cobs) for i in t, j in t];

p_niveles = contour(t, t, Z, 
    title = "Conjuntos de Nivel: Cópula Empírica",
    xlabel = "u (Viento)", 
    ylabel = "v (Oleaje)",
    color = :viridis,
    levels = 20,  # Cantidad de curvas de nivel a dibujar
    fill = true,
    colorbar_title = "C_n(u,v)"
)

# ------------------------------------------------------------------------------
# 2. Diagonal Principal: C_n(t, t)
# ------------------------------------------------------------------------------
# Representa la probabilidad conjunta P(U ≤ t, V ≤ t)
diag_princ = [Cn(i, i, cobs) for i in t]

p_diag_princ = plot(t, diag_princ, 
    title = "Diagonal Principal de la Cópula Empírica",
    xlabel = "t", 
    ylabel = "C_n(t, t)",
    label = "Cópula Empírica",
    color = :blue,
    linewidth = 2,
    legend = :topleft
)

# Agregamos las Cotas de Fréchet-Hoeffding para comparar el grado de dependencia
# Cota Superior  M(t,t) = t
plot!(p_diag_princ, t, t, 
    label = "Cota Sup", 
    linestyle = :dash, 
    color = :black)

# Cota Inferior W(t,t) = max(2t-1, 0)
plot!(p_diag_princ, t, max.(2 .* t .- 1, 0), 
    label = "Cota Inf ", 
    linestyle = :dot, 
    color = :red)


# ------------------------------------------------------------------------------
# 3. Diagonal Secundaria: C_n(t, 1-t)
# ------------------------------------------------------------------------------
# Representa la probabilidad conjunta P(U ≤ t, V ≤ 1-t)
diag_sec = [Cn(i, 1 - i, cobs) for i in t]

p_diag_sec = plot(t, diag_sec, 
    title = "Diagonal Secundaria de la Cópula Empírica",
    xlabel = "t", 
    ylabel = "C_n(t, 1-t)",
    label = "Cópula Empírica",
    color = :purple,
    linewidth = 2,
    legend = :topright
)

# Agregamos la Cota Superior de Fréchet-Hoeffding evaluada en (t, 1-t)
# M(t, 1-t) = min(t, 1-t)
plot!(p_diag_sec, t, min.(t, 1 .- t), 
    label = "Cota Sup ", 
    linestyle = :dash, 
    color = :black)

# Agregamos la Cota Inferior de Fréchet-Hoeffding evaluada en (t, 1-t)
# W(t, 1-t) = max(t + (1-t) - 1, 0) = 0
plot!(p_diag_sec, t, zeros(length(t)), 
    label = "Cota Inf", 
    linestyle = :dot, 
    color = :red)


# ==============================================================================
# ■ MODELIZACIÓN PARAMÉTRICA
# ==============================================================================
#
# 1. Propuesta y Justificación de la Cópula Subyacente (Cópula de Gumbel)
# El análisis exploratorio de las pseudo-observaciones (cópula empírica) revela 
# una clara estructura de dependencia asimétrica. Se observa una concentración 
# de probabilidad significativamente mayor en el cuadrante superior derecho (1,1) 
# en comparación con el cuadrante inferior izquierdo (0,0). Esta evidencia gráfica 
# indica la presencia de una fuerte dependencia de cola superior.

# ==============================================================================
# ■ AJUSTE DE MARGINALES
# ==============================================================================
viento_pos = viento[viento .> 0]
oleaje_pos = oleaje[oleaje .> 0]

println("Observaciones positivas para viento: ", length(viento_pos))

println("Observaciones positivas para oleaje: ", length(oleaje_pos))

# 1. VELOCIDAD DEL VIENTO ~ Gamma(α, θ)
# Estimación clásica por Máxima Verosimilitud (MLE)
gamma_viento_fit = fit_mle(Gamma, viento_pos)
α_gamma_viento = shape(gamma_viento_fit)
θ_gamma_viento = scale(gamma_viento_fit)
println("======================================")
println("AJUSTE GAMMA PARA VELOCIDAD DEL VIENTO POR MLE")
println("======================================")
println("α estimado shape:  ", round(α_gamma_viento, digits=4))
println("θ estimado scale:  ", round(θ_gamma_viento, digits=4))
println("Media teórica:     ", round(mean(gamma_viento_fit), digits=4))
println("Varianza teórica:  ", round(var(gamma_viento_fit), digits=4))
println("Log-verosimilitud: ", round(sum(logpdf.(gamma_viento_fit, viento_pos)), digits=4))
println("======================================\n")

# Gráfica del ajuste Gamma para velocidad del viento
p_gamma_viento = histogram(
    viento_pos,
    normalize = true,
    bins = :scott,
    alpha = 0.5,
    color = :skyblue,
    label = "Datos",
    xlabel = "Velocidad del viento m/s",
    ylabel = "Densidad",
    title = "Ajuste Gamma - Velocidad del Viento"
)
x_gamma_viento = range(minimum(viento_pos), maximum(viento_pos), length = 500)
plot!(
    p_gamma_viento,
    x_gamma_viento,
    pdf.(gamma_viento_fit, x_gamma_viento),
    lw = 3,
    color = :blue,
    label = "Gamma ajustada"
)

# 2. VELOCIDAD DEL VIENTO ~ Weibull(k, λ)
# Estimación clásica por Máxima Verosimilitud (MLE) 

weibull_viento_fit = fit_mle(Weibull, viento_pos)
k_weibull_viento = shape(weibull_viento_fit)
λ_weibull_viento = scale(weibull_viento_fit)
println("========================================")
println("AJUSTE WEIBULL PARA VELOCIDAD DEL VIENTO POR MLE")
println("========================================")
println("k estimado shape:  ", round(k_weibull_viento, digits=4))
println("λ estimado scale:  ", round(λ_weibull_viento, digits=4))
println("Media teórica:     ", round(mean(weibull_viento_fit), digits=4))
println("Varianza teórica:  ", round(var(weibull_viento_fit), digits=4))
println("Log-verosimilitud: ", round(sum(logpdf.(weibull_viento_fit, viento_pos)), digits=4))
println("========================================\n")


# Gráfica del ajuste Weibull para velocidad del viento
p_weibull_viento = histogram(
    viento_pos,
    normalize = true,
    bins = :scott,
    alpha = 0.5,
    color = :skyblue,
    label = "Datos",
    xlabel = "Velocidad del viento m/s",
    ylabel = "Densidad",
    title = "Ajuste Weibull - Velocidad del Viento"
)

x_weibull_viento = range(minimum(viento_pos), maximum(viento_pos), length = 500)

plot!(
    p_weibull_viento,
    x_weibull_viento,
    pdf.(weibull_viento_fit, x_weibull_viento),
    lw = 3,
    color = :darkblue,
    label = "Weibull ajustada"
)

# 3. OLEAJE ~ Weibull(k, λ)
# Estimación clásica por Máxima Verosimilitud (MLE)
weibull_oleaje_fit = fit_mle(Weibull, oleaje_pos)
k_weibull_oleaje = shape(weibull_oleaje_fit)
λ_weibull_oleaje = scale(weibull_oleaje_fit)
println("======================================")
println("AJUSTE WEIBULL PARA EL OLEAJE POR MLE")
println("======================================")
println("k estimado shape:  ", round(k_weibull_oleaje, digits=4))
println("λ estimado scale:  ", round(λ_weibull_oleaje, digits=4))
println("Media teórica:     ", round(mean(weibull_oleaje_fit), digits=4))
println("Varianza teórica:  ", round(var(weibull_oleaje_fit), digits=4))
println("Log-verosimilitud: ", round(sum(logpdf.(weibull_oleaje_fit, oleaje_pos)), digits=4))
println("======================================\n")

# Gráfica del ajuste Weibull para oleaje
p_weibull_oleaje = histogram(
    oleaje_pos,
    normalize = true,
    bins = :scott,
    alpha = 0.5,
    color = :lightgreen,
    label = "Datos",
    xlabel = "Altura del oleaje m",
    ylabel = "Densidad",
    title = "Ajuste Weibull - Oleaje"
)

x_weibull_oleaje = range(minimum(oleaje_pos), maximum(oleaje_pos), length = 500)

plot!(
    p_weibull_oleaje,
    x_weibull_oleaje,
    pdf.(weibull_oleaje_fit, x_weibull_oleaje),
    lw = 3,
    color = :darkgreen,
    label = "Weibull ajustada"
)

# Comparación gráfica de marginales ajustadas
plot(
    p_gamma_viento,
    p_weibull_viento,
    p_weibull_oleaje,
    layout = (1, 3),
    size = (1500, 400)
)

function métodoABC_Gumbel(m, k, cobs::Matrix{<:Real})

    n = size(cobs, 1)
    CnObs = zeros(n)

    for i ∈ 1:n
        CnObs[i] = Cn(cobs[i,1],cobs[i,2],cobs)
    end

    θsim = rand(Uniform(2, 5), m)

    ε = zeros(m)

    for j ∈ 1:m
        csim = Matrix(rand(GumbelCopula(2, θsim[j]), n)')
        u_obs = cobs[:, 1]
        v_obs = cobs[:, 2]
        for i in 1:n
            ε[j] += (CnObs[i] - Cn(cobs[i, 1], cobs[i, 2], csim))^2
        end
    end

    jopt = sortperm(ε)[1:k]

    return θsim[jopt], ε[jopt]
end

@time θabc, εabc = métodoABC_Gumbel(200 ,20, cobs)

θ_hat = mean(θabc)

println("θ estimado = ", θ_hat)

cop_fit = GumbelCopula(2, θ_hat)

cobs_sim = Matrix(rand(cop_fit, 15000)')

u_sim = cobs_sim[:,1]
v_sim = cobs_sim[:,2]

Spearman(cobs_sim,200)

Schweizer(cobs_sim,200)

# Scatter plot de las simulaciones
p_sim = scatter(
    u_sim,
    v_sim,
    title = "Simulación de la Cópula de Gumbel Ajustada",
    xlabel = "U simulada",
    ylabel = "V simulada",
    legend = false,
    color = :darkred,
    markersize = 2.5,
    markeralpha = 0.45,
    markerstrokewidth = 0,
    size = (700,600)
)

# ==============================================================================
# ■ DIAGONAL PRINCIPAL DE LA CÓPULA DE GUMBEL AJUSTADA
# ==============================================================================

# Grid en [0,1]
t = range(0, 1, length = 100)

# Evaluación de la diagonal principal: C(t,t)
diag_princ_fit = [cdf(cop_fit, [tt, tt]) for tt in t]

p_diag_princ_fit = plot(
    t,
    diag_princ_fit,
    lw = 3,
    color = :darkred,
    label = "Gumbel ajustada",
    xlabel = "t",
    ylabel = "C(t,t)",
    title = "Diagonal Principal de la Cópula Ajustada",
    size = (650,500)
)

# Cota superior de Fréchet-Hoeffding
plot!(
    p_diag_princ_fit,
    t,
    t,
    lw = 2,
    ls = :dash,
    color = :black,
    label = "Cota Sup"
)

# Cota inferior
plot!(
    p_diag_princ_fit,
    t,
    max.(2 .* t .- 1, 0),
    lw = 2,
    ls = :dot,
    color = :blue,
    label = "Cota Inf"
)

plot!(t, diag_princ, 
    title = "Diagonal Principal de la Cópula Empírica",
    xlabel = "t", 
    ylabel = "C_n(t, t)",
    label = "Cópula Empírica",
    color = :blue,
    linewidth = 2,
    legend = :topleft
)

# ==============================================================================
# ■ DIAGONAL SECUNDARIA DE LA CÓPULA DE GUMBEL AJUSTADA
# ==============================================================================

# Evaluación de la diagonal secundaria: C(t,1-t)
diag_sec_fit = [cdf(cop_fit, [tt, 1-tt]) for tt in t]

p_diag_sec_fit = plot(
    t,
    diag_sec_fit,
    lw = 3,
    color = :purple,
    label = "Gumbel ajustada",
    xlabel = "t",
    ylabel = "C(t,1-t)",
    title = "Diagonal Secundaria de la Cópula Ajustada",
    size = (650,500)
)

# Cota superior
plot!(
    p_diag_sec_fit,
    t,
    min.(t, 1 .- t),
    lw = 2,
    ls = :dash,
    color = :black,
    label = "Cota Sup"
)

# Cota inferior
plot!(
    p_diag_sec_fit,
    t,
    zeros(length(t)),
    lw = 2,
    ls = :dot,
    color = :blue,
    label = "Cota Inf"
)

plot!(t, diag_sec, 
    title = "Diagonal Secundaria de la Cópula Empírica",
    xlabel = "t", 
    ylabel = "C_n(t, 1-t)",
    label = "Cópula Empírica",
    color = :purple,
    linewidth = 2,
    legend = :topright
)


# ==============================================================================
# ■ COMPARACIÓN DE CURVAS DE NIVEL:
#   CÓPULA EMPÍRICA VS GUMBEL AJUSTADA
# ==============================================================================

# Grid en [0,1]
t = range(0, 1, length = 100)

# ----------------------------------------------------------------------
# Cópula empírica
# ----------------------------------------------------------------------

Z_emp = [Cn(i, j, cobs) for i in t, j in t]

p_emp_nivel = contour(
    t,
    t,
    Z_emp,
    title = "Curvas de Nivel\nCópula Empírica",
    xlabel = "u",
    ylabel = "v",
    fill = true,
    levels = 20,
    color = :viridis,
    colorbar_title = "Cₙ(u,v)",
    size = (600,500)
)

# ----------------------------------------------------------------------
# Cópula ajustada (Gumbel)
# ----------------------------------------------------------------------

Z_fit = [cdf(cop_fit, [i, j]) for i in t, j in t]

p_fit_nivel = contour(
    t,
    t,
    Z_fit,
    title = "Curvas de Nivel\nGumbel Ajustada",
    xlabel = "u",
    ylabel = "v",
    fill = true,
    levels = 20,
    color = :plasma,
    colorbar_title = "C(u,v)",
    size = (600,500)
)

# ----------------------------------------------------------------------
# Comparación conjunta
# ----------------------------------------------------------------------

plot(
    p_emp_nivel,
    p_fit_nivel,
    layout = (1,2),
    size = (1250,400)
)