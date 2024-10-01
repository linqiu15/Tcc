# DD
function potential_DD(r, Λ; I=0, full=true, particle=["ρn"])
    r = r / hbarc
    q01, q02 = q0f(mDn + mDc, mDn, mDc, mDn, mDc), q0f(mDn + mDc, mDn, mDc, mDc, mDn)
    res = 0.0
    if full || ("ρn" in particle)
        res += (gV * β)^2 * (-1 + q01^2 / mρn^2) / sqrt(16) * ft_scalar(r, Λ, mρn, q01)
    end
    if full || ("ω" in particle)
        res += -(gV * β)^2 * (-1 + q01^2 / mω^2) / sqrt(16) * ft_scalar(r, Λ, mω, q01)
    end
    if full || ("σ" in particle)
        res += -4 * gs^2 / 4 * ft_scalar(r, Λ, mσ, q01)
    end
    if full || ("ρc" in particle)
        res += -(sqrt(2) * gV * β)^2 * (-1 + q02^2 / mρc^2) / 4 * ft_scalar(r, Λ, mρc, q02) * (I - 1 / 2) * 2
    end
    return res
end


# DD*->V->D*D potential has been modified
function potential_DnDsc(r, Λ; I=0, a=0, st=[1, 0], full=true, particle=["πc"])
    os, ot = st
    r = r / hbarc
    q01, q02 = q0f(mDn + mDsc, mDn, mDsc, mDn, mDsc), q0f(mDn + mDsc, mDn, mDsc, mDsc, mDn)
    res = 0.0

    if full || ("ρn" in particle)
        res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * ft_scalar(r, Λ, mρn, q01) * os
    end
    if full || ("ω" in particle)
        res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * ft_scalar(r, Λ, mω, q01) * os
    end
    if full || ("σ" in particle)
        res += -4 * gs^2 / sqrt(16) * ft_scalar(r, Λ, mσ, q01) * os
    end
    if full || ("πc" in particle)
        res += -4 * g^2 / fπ^2 / sqrt(16) * (ft_tensor(r, Λ, mπc, q02) * ot + 1 / 3 * ft_square(r, Λ, mπc, q02) * os - a / 3 * δft(r, Λ, mπc, q02) * os) * (I - 1 / 2) * 2
    end
    if full || ("ρc" in particle)
        res += (2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (ft_tensor(r, Λ, mρc, q02) * ot - 2 / 3 * ft_square(r, Λ, mρc, q02) * os - a / 3 * δft(r, Λ, mρc, q02) * os) * (I - 1 / 2) * 2
    end
    return res
end

function potential_DcDsn(r, Λ; I=0, a=0, st=[1, 0], full=true, particle=["πc"])
    os, ot = st
    r = r / hbarc
    q01, q02 = q0f(mDc + mDsn, mDc, mDsn, mDc, mDsn), q0f(mDc + mDsn, mDc, mDsn, mDsn, mDc)
    res = 0.0

    if full || ("ρn" in particle)
        res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * ft_scalar(r, Λ, mρn, q01) * os
    end
    if full || ("ω" in particle)
        res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * ft_scalar(r, Λ, mω, q01) * os
    end
    if full || ("σ" in particle)
        res += -4 * gs^2 / sqrt(16) * ft_scalar(r, Λ, mσ, q01) * os
    end
    if full || ("πc" in particle)
        res += -4 * g^2 / fπ^2 / sqrt(16) * (ft_tensor(r, Λ, mπc, q02) * ot + 1 / 3 * ft_square(r, Λ, mπc, q02) * os - a / 3 * δft(r, Λ, mπc, q02) * os) * (I - 1 / 2) * 2
    end
    if full || ("ρc" in particle)
        res += (2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * (ft_tensor(r, Λ, mρc, q02) * ot - 2 / 3 * ft_square(r, Λ, mρc, q02) * os - a / 3 * δft(r, Λ, mρc, q02) * os) * (I - 1 / 2) * 2
    end
    return res
end

function potential_DDscross(r, Λ; I=0, a=0, st=[1, 0], full=true, particle=["πn"])
    os, ot = st
    r = r / hbarc
    q01, q02 = q0f(mDn + mDsc, mDn, mDsc, mDc, mDsn), q0f(mDn + mDsc, mDn, mDsc, mDsn, mDc)
    res = 0.0

    if full || ("ρc" in particle)
        res += 2 * β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * ft_scalar(r, Λ, mρc, q01) * os * (I - 1 / 2) * 2
    end
    if full || ("πn" in particle)
        res += 2 * g^2 / fπ^2 / sqrt(16) * (ft_tensor(r, Λ, mπn, q02) * ot + 1 / 3 * ft_square(r, Λ, mπn, q02) * os - a / 3 * δft(r, Λ, mπn, q02) * os)
    end
    if full || ("η" in particle)
        res += -2 * g^2 / (3fπ^2) / sqrt(16) * (ft_tensor(r, Λ, mη, q02) * ot + 1 / 3 * ft_square(r, Λ, mη, q02) * os - a / 3 * δft(r, Λ, mη, q02) * os)
    end
    if full || ("ρn" in particle)
        res += -(2 * gV * gλ)^2 / sqrt(16) * (ft_tensor(r, Λ, mρn, q02) * ot - 2 / 3 * ft_square(r, Λ, mρn, q02) * os - a / 3 * δft(r, Λ, mρn, q02) * os)
    end
    if full || ("ω" in particle)
        res += (2 * gV * gλ)^2 / sqrt(16) * (ft_tensor(r, Λ, mω, q02) * ot - 2 / 3 * ft_square(r, Λ, mω, q02) * os - a / 3 * δft(r, Λ, mω, q02) * os)
    end
    return res
end