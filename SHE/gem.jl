# https://github.com/xkDong14/double-Jpsi-Interaction/blob/master/Potential-r/quadgauss.jl
"""
    quadgauss(f, x::T, w::T) where {T<:Vector{Float64}}
Integration of function `f` using the Gaussian quadratures `x` with weights `w`.
`x` and `w` can be generated using, e.g., `gauss(N, a, b)` in the package `QuadGK`.
Using `quadgk` directly from that package causes memory allocation.
However, if the integration region `[a, b]` is fixed, this function does not lead to any allocation and thus is much faster.
"""
function quadgauss(f, x::T, w::T) where {T<:Vector{Float64}}
    res = zero(f(x[1]))  # zero of the same type as f(x[1]), to avoid type instability
    for i in eachindex(x)
        res += f(x[i]) * w[i]
    end
    return res
end


# common functions
λf(x, y, z) = x^2 + y^2 + z^2 - 2 * x * y - 2 * x * z - 2 * y * z
q0f(E, m1, m2, m3, m4) = abs(m3^2 - m4^2 - m1^2 + m2^2) / (2 * E)
μf(m1, m2) = m1 * m2 / (m1 + m2)

"""
Calculate the integral from zero to positive-infinity.
"""
function integauss1d(f, xc = 1.0)
    f1(x) = f(xc * x)
    f2(x) = 1 / x^2 * f(xc / x)
    tmp1 = quadgauss(f1, xxx, www)
    tmp2 = quadgauss(f2, xxx, www)
    return xc * (tmp1 + tmp2)
end

function double_factorial(n)
    if n == 1 || n == 2
        return n
    end
    return n * double_factorial(n - 2)
end

function N_gem(νn, l)
    return sqrt(2^(l + 2) * (2νn)^(l + 3 / 2) / (sqrt(π) * double_factorial(2l + 1)))
end

function φ_gem(r, νn, l)
    return N_gem(νn, l) * r^l * exp(-νn * r^2)
end

function NMatrix_gem(νn1, νn2, l)
    return (2sqrt(νn1 * νn2) / (νn1 + νn2))^(l + 3 / 2)
end

function TMatrix_gem(νn1, νn2, l, μ, hbar)
    return hbar^2 / μ * (2l + 3) * νn1 * νn2 / (νn1 + νn2) * (2 * sqrt(νn1 * νn2) / (νn1 + νn2))^(l + 3 / 2)
end

function VMatrix_gem(V, νn1, νn2, l)
    return N_gem(νn1, l) * N_gem(νn2, l) * integauss1d(r -> (r^(2l) * exp(-(νn1 + νn2) * r^2) * V(r) * r^2))
end


function νnf(n, r1, rnmax, nmax)
    return 1 / r1^2 * (r1 / rnmax)^((2n - 2) / (nmax - 1))
end