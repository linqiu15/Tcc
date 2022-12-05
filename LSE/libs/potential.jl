# ----------------------------------------------------------DD*----------------------------------------------------#
function potential_DnDsc(p, q; I=0, a=0, cl=:ss, isπ=true)
    q01, q02 = q0f(mDn + mDsc, mDn, mDsc, mDn, mDsc), q0f(mDn + mDsc, mDn, mDsc, mDsc, mDn)
    res = zero(ComplexF64)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2; cl=cl)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2; cl=cl)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2; cl=cl)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(p, q, mπc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mπc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (pwa_tensor(p, q, mρc^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mρc^2 - q02^2; cl=cl)
                                                    -
                                                    pwa_square(p, q, mρc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2; cl=cl) * (-2)
    return res
end

function potential_DcDsn(p, q; I=0, a=0, cl=:ss, isπ=true)
    q01, q02 = q0f(mDc + mDsn, mDc, mDsn, mDc, mDsn), q0f(mDc + mDsn, mDc, mDsn, mDsn, mDc)
    res = zero(ComplexF64)

    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2; cl=cl)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2; cl=cl)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2; cl=cl)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(p, q, mπc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mπc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (pwa_tensor(p, q, mρc^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mρc^2 - q02^2; cl=cl)
                                                    -
                                                    pwa_square(p, q, mρc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2; cl=cl) * (-2)
    return res
end

function potential_DDscross(p, q; I=0, a=0, cl=:ss, isπ=true)
    q01, q02 = q0f(mDn + mDsc, mDn, mDsc, mDc, mDsn), q0f(mDn + mDsc, mDn, mDsc, mDsn, mDc)
    res = zero(ComplexF64)

    res += 2 * β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q01^2; cl=cl)
    if isπ == true
        res += 2 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(p, q, mπn^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mπn^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    end
    res += -2 * g^2 / (3fπ^2) / sqrt(16) * (pwa_tensor(p, q, mη^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mη^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(p, q, mρn^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mρn^2 - q02^2; cl=cl)
                                           -
                                           pwa_square(p, q, mρn^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += ((2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(p, q, mω^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mω^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(p, q, mω^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    #ηc,J/ψ
    res += 2 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(p, q, mηc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mηc^2 - q02^2; cl=cl)) * (-2) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(p, q, mJψ^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mJψ^2 - q02^2; cl=cl)
                                           -
                                           pwa_square(p, q, mJψ^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2*(-2))
    return res
end


#----------------------------------------------------------DD*bar----------------------------------------------------#
function potential_Xn(p, q; C=1, a=0, cl=:ss, isπ=true)
    q01, q02 = q0f(mDn + mDsn, mDn, mDsn, mDn, mDsn), q0f(mDn + mDsn, mDn, mDsn, mDsn, mDn)
    res = zero(ComplexF64)

    if isπ == true
        res += (-C) * (√2 * g / fπ)^2 / 4 * (pwa_tensor(p, q, mπn^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mπn^2 - q02^2; cl=cl))
    end
    res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * (pwa_tensor(p, q, mη^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mη^2 - q02^2; cl=cl))
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * pwa_central(p, q, mρn^2 - q01^2; cl=cl)
    res += ((-C) * (2 * gλ * gV)^2 / 4 * (pwa_tensor(p, q, mρn^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mρn^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(p, q, mρn^2 - q02^2; cl=cl)))
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * pwa_central(p, q, mω^2 - q01^2; cl=cl)
    res += ((-C) * (2 * gλ * gV)^2 / 4 * (pwa_tensor(p, q, mω^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mω^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(p, q, mω^2 - q02^2; cl=cl)))
    res += -4 * gs^2 / 4 * pwa_central(p, q, mσ^2 - q01^2; cl=cl)
    #J/ψ
    res += -β^2 * gV^2 * (1 - q01^2 / mJΨ^2) / 4 * pwa_central(p, q, mJΨ^2 - q01^2; cl=cl) * 2
    res += ((-C) * (2 * gλ * gV)^2 / 4 * (pwa_tensor(p, q, mJΨ^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mJΨ^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(p, q, mJΨ^2 - q02^2; cl=cl)) * 2)
    #ηc
    res += (-C) * (√2 * g / fπ)^2 / 4 * (pwa_tensor(p, q, mηc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mηc^2 - q02^2; cl=cl)) * 2
    return res
end

function potential_Xc(p, q; C=1, a=0, cl=:ss, isπ=true)
    q01, q02 = q0f(mDc + mDsc, mDc, mDsc, mDc, mDsc), q0f(mDc + mDsc, mDc, mDsc, mDsc, mDc)
    res = zero(ComplexF64)

    if isπ == true
        res += (-C) * (√2 * g / fπ)^2 / 4 * (pwa_tensor(p, q, mπn^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mπn^2 - q02^2; cl=cl))
    end
    res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * (pwa_tensor(p, q, mη^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mη^2 - q02^2; cl=cl))
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * pwa_central(p, q, mρn^2 - q01^2; cl=cl)
    res += ((-C) * (2 * gλ * gV)^2 / 4 * (pwa_tensor(p, q, mρn^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mρn^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(p, q, mρn^2 - q02^2; cl=cl)))
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * pwa_central(p, q, mω^2 - q01^2; cl=cl)
    res += ((-C) * (2 * gλ * gV)^2 / 4 * (pwa_tensor(p, q, mω^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mω^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(p, q, mω^2 - q02^2; cl=cl)))
    res += -4 * gs^2 / 4 * pwa_central(p, q, mσ^2 - q01^2; cl=cl)
    #J/ψ
    res += -β^2 * gV^2 * (1 - q01^2 / mJΨ^2) / 4 * pwa_central(p, q, mJΨ^2 - q01^2; cl=cl) * 2
    res += ((-C) * (2 * gλ * gV)^2 / 4 * (pwa_tensor(p, q, mJΨ^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mJΨ^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(p, q, mJΨ^2 - q02^2; cl=cl)) * 2)
    #ηc
    res += (-C) * (√2 * g / fπ)^2 / 4 * (pwa_tensor(p, q, mηc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mηc^2 - q02^2; cl=cl)) * 2

    return res
end

function potential_Xcross(p, q; C=1, a=0, cl=:ss, isπ=true)
    q01, q02 = q0f(mDn + mDsn, mDn, mDsn, mDc, mDsc), q0f(mDn + mDsn, mDn, mDsn, mDsc, mDc)
    res = zero(ComplexF64)

    res += -(√2 * gV * β)^2 * (1 - q01^2 / mρc^2) / 4 * pwa_central(p, q, mρc^2 - q01^2; cl=cl)
    res += ((-C) * (2 * √2 * gλ * gV)^2 / 4 * (pwa_tensor(p, q, mρc^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(p, q, mρc^2 - q02^2; cl=cl)
                                               -
                                               pwa_square(p, q, mρc^2 - q02^2; cl=cl)))
    if isπ == true
        res += (-C) * (2 * g / fπ)^2 / 4 * (pwa_tensor(p, q, mπc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(p, q, mπc^2 - q02^2; cl=cl))
    end
    return res
end