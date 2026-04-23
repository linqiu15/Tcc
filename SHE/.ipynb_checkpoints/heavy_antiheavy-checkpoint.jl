# C:1 positive C:-1 negative
# neutral channel is lower
function potential_Xn(r, Λ; C = 1, a = 0, st = [1, 0], full = true, particle = ["π"])
    os, ot = st
    r = r / hbarc
    q01, q02 = q0f(mDn + mDsn, mDn, mDsn, mDn, mDsn), q0f(mDn + mDsn, mDn, mDsn, mDsn, mDn)
    res = 0.0

    if full || ("π" in particle)
        res += (-C) * (√2 * g / fπ)^2 / 4 * (ft_tensor(r, Λ, mπn, q02) * ot + 1 / 3 * ft_square(r, Λ, mπn, q02) * os - a / 3 * δft(r, Λ, mπn, q02) * os)
    end
    if full || ("η" in particle)
        res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * (ft_tensor(r, Λ, mη, q02) * ot + 1 / 3 * ft_square(r, Λ, mη, q02) * os - a / 3 * δft(r, Λ, mη, q02) * os)
    end
    if full || ("ρ" in particle)
        res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * ft_scalar(r, Λ, mρn, q01) * os
        # res += (-C) * (-1) * (2 * gλ * gV)^2 / 4 * (ft_tensor(r, Λ, mρn, q02) * ot - 2 / 3 * ft_square(r, Λ, mρn, q02) * os - a / 3 * δft(r, Λ, mρn, q02) * os) * (-1)
    end
    if full || ("ω" in particle)
        res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * ft_scalar(r, Λ, mω, q01) * os
        # res += (-C) * (-1) * (2 * gλ * gV)^2 / 4 * (ft_tensor(r, Λ, mω, q02) * ot - 2 / 3 * ft_square(r, Λ, mω, q02) * os - a / 3 * δft(r, Λ, mω, q02) * os) * (-1)
    end
    if full || ("σ" in particle)
        res += -4 * gs^2 / 4 * ft_scalar(r, Λ, mσ, q01) * os
    end

    return res
end

function potential_Xcross(r, Λ; C = 1, a = 0, st = [1, 0], full = true, particle = ["ρ"])
    os, ot = st
    r = r / hbarc
    q01, q02 = q0f(mDn + mDsn, mDn, mDsn, mDc, mDsc), q0f(mDn + mDsn, mDn, mDsn, mDsc, mDc)
    res = 0.0

    if full || ("ρ" in particle)
        res += -(√2 * gV * β)^2 * (1 - q01^2 / mρc^2) / 4 * ft_scalar(r, Λ, mρc, q01) * os
        # res += (-C) * (-1) * (2 * √2 * gλ * gV)^2 / 4 * (ft_tensor(r, Λ, mρc, q02) * ot - 2 / 3 * ft_square(r, Λ, mρc, q02) * os - a / 3 * δft(r, Λ, mρc, q02) * os) * (-1)
    end
    if full || ("π" in particle)
        res += (-C) * (2 * g / fπ)^2 / 4 * (ft_tensor(r, Λ, mπc, q02) * ot + 1 / 3 * ft_square(r, Λ, mπc, q02) * os - a / 3 * δft(r, Λ, mπc, q02) * os)
    end

    return res
end


function potential_Xc(r, Λ; C = 1, a = 0, st = [1, 0], full = true, particle = ["π"])
    os, ot = st
    r = r / hbarc
    q01, q02 = q0f(mDc + mDsc, mDc, mDsc, mDc, mDsc), q0f(mDc + mDsc, mDc, mDsc, mDsc, mDc)
    res = 0.0

    if full || ("π" in particle)
        res += (-C) * (√2 * g / fπ)^2 / 4 * (ft_tensor(r, Λ, mπn, q02) * ot + 1 / 3 * ft_square(r, Λ, mπn, q02) * os - a / 3 * δft(r, Λ, mπn, q02) * os)
    end
    if full || ("η" in particle)
        res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * (ft_tensor(r, Λ, mη, q02) * ot + 1 / 3 * ft_square(r, Λ, mη, q02) * os - a / 3 * δft(r, Λ, mη, q02) * os)
    end
    if full || ("ρ" in particle)
        res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * ft_scalar(r, Λ, mρn, q01) * os
        # res += (-C) * (-1) * (2 * gλ * gV)^2 / 4 * (ft_tensor(r, Λ, mρn, q02) * ot - 2 / 3 * ft_square(r, Λ, mρn, q02) * os - a / 3 * δft(r, Λ, mρn, q02) * os) * (-1)
    end
    if full || ("ω" in particle)
        res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * ft_scalar(r, Λ, mω, q01) * os
        # res += (-C) * (-1) * (2 * gλ * gV)^2 / 4 * (ft_tensor(r, Λ, mω, q02) * ot - 2 / 3 * ft_square(r, Λ, mω, q02) * os - a / 3 * δft(r, Λ, mω, q02) * os) * (-1)
    end
    if full || ("σ" in particle)
        res += -4 * gs^2 / 4 * ft_scalar(r, Λ, mσ, q01) * os
    end

    return res
end