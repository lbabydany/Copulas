# ==============================================================================
# copulas.jl
# Implementación de la Cópula Empírica y Medidas de Dependencia No Paramétricas
# ==============================================================================

using Statistics

"""
    Cn(u::Real, v::Real, cobs::AbstractMatrix{<:Real})

Evalúa la función de distribución de la cópula empírica C_n(u, v):
    C_n(u, v) = (1 / n) * ∑ 𝟙_{U_i ≤ u, V_i ≤ v}
"""
function Cn(u::Real, v::Real, cobs::AbstractMatrix{<:Real})
    n = size(cobs, 1)
    conteo = 0
    @inbounds for i in 1:n
        if cobs[i, 1] <= u && cobs[i, 2] <= v
            conteo += 1
        end
    end
    return conteo / n
end

"""
    Schweizer(cobs::AbstractMatrix{<:Real}, m::Int=100)

Calcula la medida de dependencia de Schweizer-Wolff (σ) aproximada sobre una malla
uniforme de m x m en [0, 1]²:
    σ = 12 ∬ |C_n(u, v) - u * v| du dv
Captura cualquier forma de dependencia (incluso no monótona). Retorna un valor en [0, 1].
"""
function Schweizer(cobs::AbstractMatrix{<:Real}, m::Int=100)
    # Puntos medios de los subintervalos en una partición regular de tamaño m
    paso = 1.0 / m
    puntos = range(paso / 2, 1.0 - paso / 2, length=m)
    
    integral = 0.0
    for u in puntos
        for v in puntos
            cn_val = Cn(u, v, cobs)
            integral += abs(cn_val - u * v)
        end
    end
    
    # Cuadratura simple en 2D escalada por el factor 12
    return 12.0 * integral * (paso^2)
end

"""
    Spearman(cobs::AbstractMatrix{<:Real}, m::Int=100)

Calcula la medida de concordancia Rho de Spearman (ρ_s) a partir de la cópula empírica
mediante integración numérica en una malla de m x m:
    ρ_s = 12 ∬ [C_n(u, v) - u * v] du dv
Evalúa la monotonía de la relación bivariada. Retorna un valor en [-1, 1].
"""
function Spearman(cobs::AbstractMatrix{<:Real}, m::Int=100)
    paso = 1.0 / m
    puntos = range(paso / 2, 1.0 - paso / 2, length=m)
    
    integral = 0.0
    for u in puntos
        for v in puntos
            cn_val = Cn(u, v, cobs)
            integral += (cn_val - u * v)
        end
    end
    
    return 12.0 * integral * (paso^2)
end