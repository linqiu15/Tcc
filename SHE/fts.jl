# functions for fourier transformation with monopole form factor
Y(x) = exp(-x) / x
Z(x) = (1 + 3 / x + 3 / x^2) * Y(x)
Z1(x) = (1 / x + 1 / x^2) * Y(x)
Z2(x) = (1 + x) * Y(x)

# https://link.aps.org/doi/10.1103/PhysRevD.86.014020

function H0(r, Λ, m, q0)
    u = sqrt(m^2 - q0^2)
    β = sqrt(abs(Λ^2 - m^2)) # the abs is reasonable since the form factor is squared
    λ = sqrt(Λ^2 - q0^2)
    res = Y(u * r) - λ / u * Y(λ * r) - r * β^2 / (2 * u) * Y(λ * r)

    return res
end

function H1(r, Λ, m, q0)
    u = sqrt(m^2 - q0^2)
    β = sqrt(abs(Λ^2 - m^2))
    λ = sqrt(Λ^2 - q0^2)
    res = Y(u * r) - λ / u * Y(λ * r) - r * λ^2 * β^2 / (2 * u^3) * Y(λ * r)

    return res
end

function H2(r, Λ, m, q0)
    u = sqrt(m^2 - q0^2)
    β = sqrt(abs(Λ^2 - m^2))
    λ = sqrt(Λ^2 - q0^2)
    res = Z1(u * r) - λ^3 / u^3 * Z1(λ * r) - λ * β^2 / (2 * u^3) * Y(λ * r)

    return res
end

function H3(r, Λ, m, q0)
    u = sqrt(m^2 - q0^2)
    β = sqrt(abs(Λ^2 - m^2))
    λ = sqrt(Λ^2 - q0^2)
    res = Z(u * r) - λ^3 / u^3 * Z(λ * r) - λ * β^2 / (2 * u^3) * Z2(λ * r)

    return res
end

function M0(r, Λ, m, q0)
    θ = sqrt(q0^2 - m^2)
    β = sqrt(abs(Λ^2 - m^2))
    λ = sqrt(Λ^2 - q0^2)
    # there is an extra minus sign in the literature
    res = 1 / (θ * r) * (cos(θ * r) - exp(-λ * r)) - β^2 / (2 * θ * λ) * exp(-λ * r)

    return res
end

function M1(r, Λ, m, q0)
    θ = sqrt(q0^2 - m^2)
    β = sqrt(abs(Λ^2 - m^2))
    λ = sqrt(Λ^2 - q0^2)
    res = -1 / (θ * r) * (cos(θ * r) - exp(-λ * r)) - λ * β^2 / (2 * θ^3) * exp(-λ * r)

    return res
end

function M3(r, Λ, m, q0)
    θ = sqrt(q0^2 - m^2)
    β = sqrt(abs(Λ^2 - m^2))
    λ = sqrt(Λ^2 - q0^2)
    res = (-(cos(θ * r) - 3 * sin(θ * r) / (θ * r) - 3 * cos(θ * r) / (θ^2 * r^2)) / (θ * r) - λ^3 / θ^3 * Z(λ * r) - λ * β^2 / (2 * θ^3) * Z2(λ * r))

    return res
end


# fourier transformation of δ-term
function δft(r, Λ, m, q0)
    β = sqrt(abs(Λ^2 - m^2))
    λ = sqrt(Λ^2 - q0^2)
    return β^4 * exp(-λ * r) / (8 * π * λ)
end



# fourier transformation for different tensor terms
function ft_scalar(r, Λ, m, q0) #only for heavy meson exchange
    u = sqrt(m^2 - q0^2)
    return u * H0(r, Λ, m, q0) / (4 * π)
end

function ft_square(r, Λ, m, q0)
    res = 0.0
    uθ = sqrt(abs(m^2 - q0^2))
    if m^2 >= q0^2
        res = -uθ^3 * H1(r, Λ, m, q0) / (4 * π)
    else
        res = -uθ^3 * M1(r, Λ, m, q0) / (4 * π)
    end
    return res
end

function ft_tensor(r, Λ, m, q0)
    res = 0.0
    uθ = sqrt(abs(m^2 - q0^2))
    if m^2 >= q0^2
        res = -uθ^3 / (12 * π) * H3(r, Λ, m, q0)
    else
        res = -uθ^3 / (12 * π) * M3(r, Λ, m, q0)
    end
    return res
end


# kinematic functions
λf(x, y, z) = x^2 + y^2 + z^2 - 2 * x * y - 2 * x * z - 2 * y * z
q0f(E, m1, m2, m3, m4) = abs(m3^2 - m4^2 - m1^2 + m2^2) / (2 * E)
μf(m1, m2) = m1 * m2 / (m1 + m2)