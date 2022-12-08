# ----------------------------------------------------------DD*----------------------------------------------------#
function potential_DnDsc(E, p, q, pon, qon; I=0, a=0.0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDn, mDsc), q0f(E, mDn, mDsc, mDsc, mDn)
    res = zero(ComplexF64)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπc^2 - q02^2; a=a1) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(pon, qon, mρc^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * (-2)
    return res
end

function potential_DcDsn(E, p, q, pon, qon; I=0, a=0.0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDc, mDsn, mDc, mDsn), q0f(E, mDc, mDsn, mDsn, mDc)
    res = zero(ComplexF64)

    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπc^2 - q02^2; a=a1) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(pon, qon, mρc^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * (-2)
    return res
end

function potential_DDscross(E, p, q, pon, qon; I=0, a=0.0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDc, mDsn), q0f(E, mDn, mDsc, mDsn, mDc)
    res = zero(ComplexF64)

    res += 2 * β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q01^2)
    if isπ == true
        res += 2 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπn^2 - q02^2; a=a1) * (I - 1 / 2) * 2
    end
    res += -2 * g^2 / (3fπ^2) / sqrt(16) * pwa_tensor_ss(pon, qon, mη^2 - q02^2; a=a) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(pon, qon, mρn^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    res += ((2 * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(pon, qon, mω^2 - q02^2; a=a) * (I - 1 / 2) * 2)
    #ηc,J/ψ
    res += 2 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss(pon, qon, mηc^2 - q02^2; a=a) * (-2) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(pon, qon, mJψ^2 - q02^2; a=a) * (I - 1 / 2) * 2 * (-2))
    return res
end



#----------------------------------------------------------DD*bar----------------------------------------------------#
function potential_Xn(E, p, q, pon, qon; C=1, a=0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsn, mDn, mDsn), q0f(E, mDn, mDsn, mDsn, mDn)
    res = zero(ComplexF64)

    if isπ == true
        res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss_deform(p, q, mπn^2 - q02^2; a=a1)
    end
    res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(pon, qon, mη^2 - q02^2; a=a)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * pwa_central(p, q, mρn^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(pon, qon, mρn^2 - q02^2; a=a)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * pwa_central(p, q, mω^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(pon, qon, mω^2 - q02^2; a=a)
    res += -4 * gs^2 / 4 * pwa_central(p, q, mσ^2 - q01^2)
    #J/ψ
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / 4 * pwa_central(p, q, mJψ^2 - q01^2) * 2
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(pon, qon, mJψ^2 - q02^2; a=a) * 2
    #ηc
    res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(pon, qon, mηc^2 - q02^2; a=a) * 2
    return res
end

function potential_Xc(E, p, q, pon, qon; C=1, a=0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDc, mDsc, mDc, mDsc), q0f(E, mDc, mDsc, mDsc, mDc)
    res = zero(ComplexF64)

    if isπ == true
        res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss_deform(p, q, mπn^2 - q02^2; a=a1)
    end
    res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(pon, qon, mη^2 - q02^2; a=a)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * pwa_central(p, q, mρn^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(pon, qon, mρn^2 - q02^2; a=a)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * pwa_central(p, q, mω^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(pon, qon, mω^2 - q02^2; a=a)
    res += -4 * gs^2 / 4 * pwa_central(p, q, mσ^2 - q01^2)
    #J/ψ
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / 4 * pwa_central(p, q, mJψ^2 - q01^2) * 2
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(pon, qon, mJψ^2 - q02^2; a=a) * 2
    #ηc
    res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(pon, qon, mηc^2 - q02^2; a=a) * 2

    return res
end

function potential_Xcross(E, p, q, pon, qon; C=1, a=0, a1=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsn, mDc, mDsc), q0f(E, mDn, mDsn, mDsc, mDc)
    res = zero(ComplexF64)

    res += -(√2 * gV * β)^2 * (1 - q01^2 / mρc^2) / 4 * pwa_central(p, q, mρc^2 - q01^2)
    res += (-C) * (2 * √2 * gλ * gV)^2 / 4 * pwa_vvv_ss(pon, qon, mρc^2 - q02^2; a=a)
    if isπ == true
        res += (-C) * (2 * g / fπ)^2 / 4 * pwa_tensor_ss(p, q, mπc^2 - q02^2; a=a1)
    end
    return res
end