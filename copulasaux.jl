# ==============================================================================
# copulasaux.jl
# Funciones auxiliares para el cálculo de pseudo-observaciones y rangos
# ==============================================================================

using Statistics

"""
    rangos(x::AbstractVector{<:Real})

Calcula el vector de rangos de los elementos de x.
En caso de empates, asigna el rango promedio (método estándar de la función de distribución empírica).
"""
function rangos(x::AbstractVector{<:Real})
    n = length(x)
    indices = sortperm(x)
    ranks = zeros(Float64, n)
    
    i = 1
    while i <= n
        j = i
        while j < n && x[indices[j+1]] == x[indices[i]]
            j += 1
        end
        # Rango promedio para empates
        rango_prom = (i + j) / 2.0
        for k in i:j
            ranks[indices[k]] = rango_prom
        end
        i = j + 1
    end
    return ranks
end

"""
    copulaPseudObs(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})

Calcula la matriz de pseudo-observaciones empíricas en el cuadrado unitario [0, 1]²:
    U_i = R_i / (n + 1)
    V_i = S_i / (n + 1)
donde R_i y S_i representan los rangos muestrales de x e y respectivamente.
Retorna una matriz de dimensión (n, 2).
"""
function copulaPseudObs(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) == length(y) "Los vectores x e y deben tener la misma dimensión."
    n = length(x)
    
    # Se usa denominador n + 1 para mantener las observaciones dentro de (0, 1)
    # y evitar divergencias en evaluaciones de cópulas paramétricas.
    u = rangos(x) ./ (n + 1)
    v = rangos(y) ./ (n + 1)
    
    return hcat(u, v)
end