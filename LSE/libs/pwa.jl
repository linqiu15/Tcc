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

function pwa_tensor_ss(p, q, usq; a=0.0)::ComplexF64
    return 1 / 3 * (-usq * pwa_central(p, q, usq) + a)
end

function pwa_vvv_ss(p, q, usq; a=0.0)::ComplexF64
    # -ϵ×qϵ×q/(q^2+m^2)
    return -2 / 3 * (-usq * pwa_central(p, q, usq) + a)
end


################################################# VV ###########################
function pwa_VVV_ss(p, q, usq; a=0.0, J=0)::ComplexF64
    # (ϵ1×ϵ3×q)⋅(ϵ2×ϵ4×q)/(q^2+m^2)
    tmp = 2 / 3 * (-usq * pwa_central(p, q, usq) + a)
    if J == 0
        return 2 * tmp
    elseif J == 1
        return tmp
    end
    return -1 * tmp
end

function pwa_Tensor_ss(p, q, usq; a=0.0, J=0)::ComplexF64
    # (ϵ1×ϵ3⋅q)⋅(ϵ2×ϵ4⋅q)/(q^2+m^2)
    tmp = 1 / 3 * (-usq * pwa_central(p, q, usq) + a)
    if J == 0
        return 2 * tmp
    elseif J == 1
        return tmp
    end
    return -1 * tmp
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


function pwa_tensor_ss_deform(p, q, usq; a=1.0)::ComplexF64
    return 1 / 3 * (1 - a) * pwa_contact(p, q, usq) - usq / 3 * pwa_central_deform(p, q, usq)
end