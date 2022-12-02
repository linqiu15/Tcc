# ----------------------------------------------------------DD*----------------------------------------------------#
function potential_DnDsc(E, p, q; I=0, a=0.0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDn, mDsc), q0f(E, mDn, mDsc, mDsc, mDn)
    res = zero(ComplexF64)
    # res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    # res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπc^2 - q02^2; a=a1) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * pwa_tensor_ss(p, q, mρc^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * (-2)
    return res
end

function potential_DcDsn(E, p, q; I=0, a=0.0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDc, mDsn, mDc, mDsn), q0f(E, mDc, mDsn, mDsn, mDc)
    res = zero(ComplexF64)

    # res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    # res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπc^2 - q02^2; a=a1) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * pwa_tensor_ss(p, q, mρc^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * (-2)
    return res
end

function potential_DDscross(E, p, q; I=0, a=0.0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDc, mDsn), q0f(E, mDn, mDsc, mDsn, mDc)
    res = zero(ComplexF64)

    res += 2 * β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q01^2)
    if isπ == true
        res += 2 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπn^2 - q02^2; a=a1) * (I - 1 / 2) * 2
    end
    res += -2 * g^2 / (3fπ^2) / sqrt(16) * pwa_tensor_ss(p, q, mη^2 - q02^2; a=a) * (I - 1 / 2) * 2
    # res += (-(2 * gV * gλ)^2 / sqrt(16) * pwa_tensor_ss(p, q, mρn^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    # res += ((2 * gV * gλ)^2 / sqrt(16) * pwa_tensor_ss(p, q, mω^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    #ηc,J/ψ
    res += 2 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss(p, q, mηc^2 - q02^2; a=a) * (-2) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * pwa_tensor_ss(p, q, mJψ^2 - q02^2; a=a) * (I - 1 / 2) * 2 * (-2))
    return res
end