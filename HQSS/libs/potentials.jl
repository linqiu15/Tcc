# ----------------------------------------------------------DD----------------------------------------------------#
function potential_DDc(E, p, q; I=0)
    q01, q02 = q0f(E, mDn, mDc, mDn, mDc), q0f(E, mDn, mDc, mDc, mDn)
    res = zero(ComplexF64)
    res += -gV^2 * β^2 * (1 - q01^2 / mρn^2) / 4 * pwa_central(p, q, mρn^2 - q01^2)
    res += gV^2 * β^2 * (1 - q01^2 / mω^2) / 4 * pwa_central(p, q, mω^2 - q01^2)
    res += -4 * gs^2 / 4 * pwa_central(p, q, mσ^2 - q01^2)
    res += (√2 * gV * β)^2 * (1 - q02^2 / mρc^2) / 4 * pwa_central(p, q, mρc^2 - q02^2) * (I - 1 / 2) * 2
    res += -gV^2 * β^2 * (1 - q01^2 / mJψ^2) / 4 * pwa_central(p, q, mJψ^2 - q01^2)*(-2)
    return res
end