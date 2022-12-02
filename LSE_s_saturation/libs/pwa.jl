function pwa_contact(p, q, usq)::ComplexF64
    return 1.0
end

# somewhat this is better than the next one
function pwa_central(p, q, usq)::ComplexF64
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if iszero(p) || iszero(q)
        tmp = 1 / (p^2 + q^2 + usq)
    else
        tmp = -1 / (4 * p * q) * (log(ξ - 1) - log(ξ + 1))
    end
    return tmp
end

# function pwa_central(p, q, usq)::ComplexF64
#     ξ = (p^2 + q^2 + usq) / (2 * p * q)
#     tmp = zero(ComplexF64)
#     if iszero(p) || iszero(q)
#         tmp = 1 / (p^2 + q^2 + usq)
#     else
#         tmp = -1 / (4 * p * q) * log(((p-q)^2+usq)/((p+q)^2+usq))
#     end
#     return tmp
# end

function pwa_central_deform(p, q, usq; a=20.0)
    tmp = zero(ComplexF64)
    if iszero(p) || iszero(q)
        tmp = 1 / (p^2 + q^2 + usq)
    else
        tmp += 1 / 2 * quadgk(z -> 1 / (p^2 + q^2 - 2 * p * q * z + usq), -1, -1 - a * im)[1]
        tmp += 1 / 2 * quadgk(z -> 1 / (p^2 + q^2 - 2 * p * q * z + usq), -1 - a * im, 1 - a * im)[1]
        tmp += 1 / 2 * quadgk(z -> 1 / (p^2 + q^2 - 2 * p * q * z + usq), 1 - a * im, 1)[1]
    end
    return tmp
end

function pwa_tensor_ss(p, q, usq; a=0.0)::ComplexF64
    return 1 / 3 * (1 - a) * pwa_contact(p, q, usq) - usq / 3 * pwa_central(p, q, usq)
end

function pwa_tensor_ss_deform(p, q, usq; a=0.0)::ComplexF64
    return 1 / 3 * (1 - a) * pwa_contact(p, q, usq) - usq / 3 * pwa_central_deform(p, q, usq)
end