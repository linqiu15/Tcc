############################################ for heavier exchanged meson 
function pwa_contact(p, q, usq)::ComplexF64
    return 1.0
end

function pwa_central(p, q, usq)::ComplexF64
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if iszero(p) || iszero(q)
        tmp = 1 / (p^2 + q^2 + usq)
    else
        tmp = -1 / (4 * p * q) * log(((p - q)^2 + usq) / ((p + q)^2 + usq))
    end
    return tmp
end

function pwa_square(p, q, usq)::ComplexF64
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if iszero(p) || iszero(q)
        tmp = (p^2 + q^2) / (p^2 + q^2 + usq)
    else
        tmp = -1 / (4 * p * q) * (-4 * p * q + (p^2 + q^2 - 2 * p * q * ξ) * log(((p - q)^2 + usq) / ((p + q)^2 + usq)))
    end
    return tmp
end

function pwa_tensor_ss(p, q, usq; a=0.0)::ComplexF64
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if iszero(p) || iszero(q)
        tmp = 1 / 3 * (p^2 + q^2) / (p^2 + q^2 + usq)
    else
        tmp = -1 / (12 * p * q) * (-4 * p * q + (p^2 + q^2 - 2 * p * q * ξ) * log(((p - q)^2 + usq) / ((p + q)^2 + usq)))
    end
    return tmp - a / 3 * pwa_contact(p, q, usq)
end

function pwa_vvv_ss(p, q, usq; a=0.0)::ComplexF64
    return pwa_tensor_ss(p, q, usq; a=a) - pwa_square(p, q, usq)
end


################################################# for ope###########################
# function pwa_central_deform(p, q, usq; a=20.0)
#     tmp = zero(ComplexF64)
#     if iszero(p) || iszero(q)
#         tmp = 1 / (p^2 + q^2 + usq)
#     else
#         tmp += 1 / 2 * quadgk(z -> 1 / (p^2 + q^2 - 2 * p * q * z + usq), -1, -1 - a * im)[1]
#         tmp += 1 / 2 * quadgk(z -> 1 / (p^2 + q^2 - 2 * p * q * z + usq), -1 - a * im, 1 - a * im)[1]
#         tmp += 1 / 2 * quadgk(z -> 1 / (p^2 + q^2 - 2 * p * q * z + usq), 1 - a * im, 1)[1]
#     end
#     return tmp
# end

function pwa_central_deform(p, q, usq)
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if iszero(p) || iszero(q)
        tmp = 1 / (p^2 + q^2 + usq)
    else
        if imag(ξ) < 0 && abs(real(ξ)) < 1
            tmp = -1 / (4 * p * q) * (log(((p - q)^2 + usq) / ((p + q)^2 + usq)) + im * 2 * π)
        else
            tmp = -1 / (4 * p * q) * log(((p - q)^2 + usq) / ((p + q)^2 + usq))
        end
    end
    return tmp
end


function pwa_tensor_ss_deform(p, q, usq; a=0.0)::ComplexF64
    return 1 / 3 * (1 - a) * pwa_contact(p, q, usq) - usq / 3 * pwa_central_deform(p, q, usq)
end