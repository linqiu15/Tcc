# ----------------------------------------------------------DD*----------------------------------------------------#
function potential_DnDsc(E, p, q, pon,qon; I=0, a=0, a1=0.0, cl=:ss, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDn, mDsc), q0f(E, mDn, mDsc, mDsc, mDn)
    res = zero(ComplexF64)
    # res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2; cl=cl)
    # res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2; cl=cl)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2; cl=cl)
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (pwa_tensor(pon, qon, mρc^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(pon, qon, mρc^2 - q02^2; cl=cl)
                                                    -
                                                    pwa_square(pon, qon, mρc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2; cl=cl) * (-2)
    if (cl==:dd&&abs(res)>1)
        res=0.0
    end
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor_deform(p, q, mπc^2 - q02^2; cl=cl) - a1 / 3 * pwa_contact(p, q, mπc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    end
    return res
end

function potential_DcDsn(E, p, q, pon,qon; I=0, a=0, a1=0.0, cl=:ss, isπ=true)
    q01, q02 = q0f(E, mDc, mDsn, mDc, mDsn), q0f(E, mDc, mDsn, mDsn, mDc)
    res = zero(ComplexF64)

    # res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2; cl=cl)
    # res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2; cl=cl)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2; cl=cl)
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (pwa_tensor(pon, qon, mρc^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(pon, qon, mρc^2 - q02^2; cl=cl)
                                                    -
                                                    pwa_square(pon, qon, mρc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2; cl=cl) * (-2)
    if (cl==:dd&&abs(res)>1)
        res=0.0
    end
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor_deform(p, q, mπc^2 - q02^2; cl=cl) - a1 / 3 * pwa_contact(p, q, mπc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    end
    return res
end

function potential_DDscross(E, p, q, pon,qon; I=0, a=0, a1=0.0, cl=:ss, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDc, mDsn), q0f(E, mDn, mDsc, mDsn, mDc)
    res = zero(ComplexF64)

    res += 2 * β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q01^2; cl=cl)   
    res += -2 * g^2 / (3fπ^2) / sqrt(16) * (pwa_tensor(pon, qon, mη^2 - q02^2; cl=cl) - a / 3 * pwa_contact(pon, qon, mη^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    # res += (-(2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(pon, qon, mρn^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(pon, qon, mρn^2 - q02^2; cl=cl)
    #                                        -
    #                                        pwa_square(pon, qon, mρn^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    # res += ((2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(pon, qon, mω^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(pon, qon, mω^2 - q02^2; cl=cl)
    #                                       -
    #                                       pwa_square(pon, qon, mω^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    #ηc,J/ψ
    res += 2 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(pon, qon, mηc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(pon, qon, mηc^2 - q02^2; cl=cl)) * (-2) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(pon, qon, mJψ^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(pon, qon, mJψ^2 - q02^2; cl=cl)
                                           -
                                           pwa_square(pon, qon, mJψ^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2 * (-2))
    if (cl==:dd&&abs(res)>1)
        res=0.0
    end
    if isπ == true
        res += 2 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor_deform(p, q, mπn^2 - q02^2; cl=cl) - a1 / 3 * pwa_contact(p, q, mπn^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    end
    return res
end