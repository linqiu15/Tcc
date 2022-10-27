# ----------------------------------------------------------DD*----------------------------------------------------#
function potential_DnDsc_integrand(z, E, p, q; I=0, a=0, a1=0.0, cl=:ss)
    q01, q02 = q0f(E, mDn, mDsc, mDn, mDsc), q0f(E, mDn, mDsc, mDsc, mDn)
    res = zero(ComplexF64)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(z, p, q, mρn^2 - q01^2; cl=cl)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(z, p, q, mω^2 - q01^2; cl=cl)
    res += -4 * gs^2 / sqrt(16) * pwa_central(z, p, q, mσ^2 - q01^2; cl=cl)
    res += -4 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(z, p, q, mπc^2 - q02^2; cl=cl) - a1 / 3 * pwa_contact(z, p, q, mπc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (pwa_tensor(z, p, q, mρc^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(z, p, q, mρc^2 - q02^2; cl=cl)
                                                    -
                                                    pwa_square(z, p, q, mρc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(z, p, q, mJψ^2 - q01^2; cl=cl) * (-2)
    return res
end

function potential_DcDsn_integrand(z, E, p, q; I=0, a=0, a1=0.0, cl=:ss)
    q01, q02 = q0f(E, mDc, mDsn, mDc, mDsn), q0f(E, mDc, mDsn, mDsn, mDc)
    res = zero(ComplexF64)

    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(z, p, q, mρn^2 - q01^2; cl=cl)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(z, p, q, mω^2 - q01^2; cl=cl)
    res += -4 * gs^2 / sqrt(16) * pwa_central(z, p, q, mσ^2 - q01^2; cl=cl)
    res += -4 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(z, p, q, mπc^2 - q02^2; cl=cl) - a1 / 3 * pwa_contact(z, p, q, mπc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (pwa_tensor(z, p, q, mρc^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(z, p, q, mρc^2 - q02^2; cl=cl)
                                                    -
                                                    pwa_square(z, p, q, mρc^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(z, p, q, mJψ^2 - q01^2; cl=cl) * (-2)
    return res
end

function potential_DDscross_integrand(z, E, p, q; I=0, a=0, a1=0.0, cl=:ss)
    q01, q02 = q0f(E, mDn, mDsc, mDc, mDsn), q0f(E, mDn, mDsc, mDsn, mDc)
    res = zero(ComplexF64)

    res += 2 * β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(z, p, q, mρc^2 - q01^2; cl=cl)
    res += 2 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(z, p, q, mπn^2 - q02^2; cl=cl) - a1 / 3 * pwa_contact(z, p, q, mπn^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    res += -2 * g^2 / (3fπ^2) / sqrt(16) * (pwa_tensor(z, p, q, mη^2 - q02^2; cl=cl) - a / 3 * pwa_contact(z, p, q, mη^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(z, p, q, mρn^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(z, p, q, mρn^2 - q02^2; cl=cl)
                                           -
                                           pwa_square(z, p, q, mρn^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    res += ((2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(z, p, q, mω^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(z, p, q, mω^2 - q02^2; cl=cl)
                                          -
                                          pwa_square(z, p, q, mω^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2)
    #ηc,J/ψ
    res += 2 * g^2 / fπ^2 / sqrt(16) * (pwa_tensor(z, p, q, mηc^2 - q02^2; cl=cl) - a / 3 * pwa_contact(z, p, q, mηc^2 - q02^2; cl=cl)) * (-2) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * (pwa_tensor(z, p, q, mJψ^2 - q02^2; cl=cl) + 2 * a / 3 * pwa_contact(z, p, q, mJψ^2 - q02^2; cl=cl)
                                           -
                                           pwa_square(z, p, q, mJψ^2 - q02^2; cl=cl)) * (I - 1 / 2) * 2 * (-2))
    return res
end


function potential_DnDsc(E, p, q; I=0, a=0, a1=0.0, cl=:ss, r=100.0)
    res = zero(ComplexF64)
    res += quadgk(z -> potential_DnDsc_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl),  -1.0, -1.0 - im * r)[1]
    res += quadgk(z -> potential_DnDsc_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl),  -1.0 - im * r, 1.0 - im * r)[1]
    res += quadgk(z -> potential_DnDsc_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl),  1.0 - im * r, 1.0)[1]
    return res
end

function potential_DcDsn(E, p, q; I=0, a=0, a1=0.0, cl=:ss, r=100.0)
    res = zero(ComplexF64)
    res += quadgk(z -> potential_DcDsn_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl), -1.0, -1.0 - im * r)[1]
    res += quadgk(z -> potential_DcDsn_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl), -1.0 - im * r, 1.0 - im * r)[1]
    res += quadgk(z -> potential_DcDsn_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl),  1.0 - im * r, 1.0)[1]
    return res
end

function potential_DDscross(E, p, q; I=0, a=0, a1=0.0, cl=:ss, r=100.0)
    res = zero(ComplexF64)
    res += quadgk(z -> potential_DDscross_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl),  -1.0, -1.0 - im * r)[1]
    res += quadgk(z -> potential_DDscross_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl), -1.0 - im * r, 1.0 - im * r)[1]
    res += quadgk(z -> potential_DDscross_integrand(z, E, p, q; I=I, a=a, a1=a1, cl=cl),  1.0 - im * r, 1.0)[1]
    return res
end
