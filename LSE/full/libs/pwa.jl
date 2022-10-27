function pwa_tensor(z, p, q, usq; cl=:ss)
    res = zero(ComplexF64)
    if cl == :ss
        res = 1 / 6 * (p^2 + q^2 - 2 * p * q * z) / (p^2 + q^2 - 2 * p * q * z + usq)
    elseif cl == :sd
        res = -√2 / 6 * (p^2 + q^2 * (3z^2 - 1) / 2 - 2 * p * q * z) / (p^2 + q^2 - 2 * p * q * z + usq)
    elseif cl == :ds
        res = -√2 / 6 * (q^2 + p^2 * (3z^2 - 1) / 2 - 2 * p * q * z) / (p^2 + q^2 - 2 * p * q * z + usq)
    elseif cl == :dd
        res = 1 / 12 * (2 * (p^2 + q^2) * (3z^2 - 1) - p * q * z * (9z^2 - 1)) / (p^2 + q^2 - 2 * p * q * z + usq)
    end
    return res
end

function pwa_square(z, p, q, usq; cl=:ss)
    res = zero(ComplexF64)
    if cl == :ss
        res = 1 / 2 * (p^2 + q^2 - 2 * p * q * z) / (p^2 + q^2 - 2 * p * q * z + usq)
    elseif cl == :dd
        res = 1 / 4 * (3z^2 - 1) * (p^2 + q^2 - 2 * p * q * z) / (p^2 + q^2 - 2 * p * q * z + usq)
    end
    return res
end

function pwa_central(z, p, q, usq; cl=:ss)
    res = zero(ComplexF64)
    if cl == :ss
        res = 1 / 2 / (p^2 + q^2 - 2 * p * q * z + usq)
    elseif cl == :dd
        res = 1 / 4 * (3z^2 - 1) / (p^2 + q^2 - 2 * p * q * z + usq)
    end
    return res
end

function pwa_contact(z, p, q, usq; cl=:ss)
    res = zero(ComplexF64)
    if cl == :ss
        res = 1
    end
    return res
end