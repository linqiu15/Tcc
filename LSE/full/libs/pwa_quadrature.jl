module pwaquadrature
export pwa_tensor_deform,pwa_tensor,pwa_central,pwa_square,pwa_contact
using QuadGK
xxx,www=gauss(200,0,1)

function quadgk1(f,a::S,b::S,x::T, w::T) where {T<:Vector{Float64},S<:ComplexF64}
    res=zero(ComplexF64)
    for i in eachindex(x)
        res += f(a+(b-a)*x[i]) * w[i]
    end
    return res*(b-a)
end


# not working great
function pwa_tensor_deform(p, q, usq; cl = :ss,a=1.0)
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if cl == :ss
        if iszero(p) || iszero(q)
            tmp = 1 / 3 * (p^2 + q^2) / (p^2 + q^2 + usq)
        else
            tmp +=1/6*quadgk1(z->(p^2+q^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),-1+0.0im,-1-a*im,xxx,www)
            tmp +=1/6*quadgk1(z->(p^2+q^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),-1-a*im,1-a*im,xxx,www)
            tmp +=1/6*quadgk1(z->(p^2+q^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),1-a*im,1+0.0im,xxx,www)
        end
    elseif cl == :sd
        if iszero(p) || iszero(q)
            tmp = -√2 / 3 * p^2 / (p^2 + q^2 + usq)
        else
            tmp +=-√2/6*quadgk1(z->(q^2/2*(3*z^2-1)+p^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),-1+0.0im,-1-a*im,xxx,www)
            tmp +=-√2/6*quadgk1(z->(q^2/2*(3*z^2-1)+p^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),-1-a*im,1-a*im,xxx,www)
            tmp +=-√2/6*quadgk1(z->(q^2/2*(3*z^2-1)+p^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),1-a*im,1+0.0im,xxx,www)
        end
    elseif cl == :ds
        if iszero(p) || iszero(q)
            tmp = -√2 / 3 * q^2 / (p^2 + q^2 + usq)
        else
            tmp +=-√2/6*quadgk1(z->(p^2/2*(3*z^2-1)+q^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),-1+0.0im,-1-a*im,xxx,www)
            tmp +=-√2/6*quadgk1(z->(p^2/2*(3*z^2-1)+q^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),-1-a*im,1-a*im,xxx,www)
            tmp +=-√2/6*quadgk1(z->(p^2/2*(3*z^2-1)+q^2-2*p*q*z)/(p^2+q^2-2*p*q*z+usq),1-a*im,1+0.0im,xxx,www)
        end
    else
        if iszero(p) || iszero(q)
            tmp = 0
        else
            tmp +=1/12*quadgk1(z->(2*(p^2+q^2)*(3*z^2-1)-p*q*z*(9*z^2-1))/(p^2+q^2-2*p*q*z+usq),-1+0.0im,-1-a*im,xxx,www)
            tmp +=1/12*quadgk1(z->(2*(p^2+q^2)*(3*z^2-1)-p*q*z*(9*z^2-1))/(p^2+q^2-2*p*q*z+usq),-1-a*im,1-a*im,xxx,www)
            tmp +=1/12*quadgk1(z->(2*(p^2+q^2)*(3*z^2-1)-p*q*z*(9*z^2-1))/(p^2+q^2-2*p*q*z+usq),1-a*im,1+0.0im,xxx,www)
        end
    end
    return tmp
end


# only for exchanged mesons heavier than pion
function pwa_tensor(p, q, usq; cl = :ss)
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if cl == :ss
        if iszero(p) || iszero(q)
            tmp = 1 / 3 * (p^2 + q^2) / (p^2 + q^2 + usq)
        else
            tmp = -1 / (12 * p * q) * (-4 * p * q + (p^2 + q^2 - 2 * p * q * ξ) *(log(1-ξ)-log(-1-ξ)))
        end
    elseif cl == :sd
        if iszero(p) || iszero(q)
            tmp = -√2 / 3 * p^2 / (p^2 + q^2 + usq)
        else
            tmp = √2 / (12 * p * q) * (q * (-4 * p + 3 * q * ξ) + (p^2 - 2 * p * q * ξ + 1 / 2 * q^2 * (-1 + 3 * ξ^2)) *(log(1-ξ)-log(-1-ξ)))
        end
    elseif cl == :ds
        if iszero(p) || iszero(q)
            tmp = -√2 / 3 * q^2 / (p^2 + q^2 + usq)
        else
            tmp = √2 / (12 * p * q) * (p * (-4 * q + 3 * p * ξ) + (q^2 - 2 * p * q * ξ + 1 / 2 * p^2 * (-1 + 3 * ξ^2)) *(log(1-ξ)-log(-1-ξ)))
        end
    else
        if iszero(p) || iszero(q)
            tmp = 0
        else
            tmp = -1 / (24 * p * q) * (2 * (6 * p^2 * ξ + 6 * q^2 * ξ - p * q * (2 + 9 * ξ^2)) + (2 * q^2 * (-1 + 3 * ξ^2) + p^2 * (-2 + 6 * ξ^2) + p * q * (ξ - 9 * ξ^3)) *(log(1-ξ)-log(-1-ξ)))
        end
    end
    return tmp
end

function pwa_central(p, q, usq; cl = :ss)
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if cl == :ss
        if iszero(p) || iszero(q)
            tmp = 1 / (p^2 + q^2 + usq)
        else
            tmp = -1 / (4 * p * q) *(log(1-ξ)-log(-1-ξ))
        end
    elseif cl == :dd
        if iszero(p) || iszero(q)
            tmp = 0
        else
            tmp = -1 / (8 * p * q) * (6 * ξ + (-1 + 3 * ξ^2) *(log(1-ξ)-log(-1-ξ)))
        end
    end
    return tmp
end

function pwa_square(p, q, usq; cl = :ss)
    ξ = (p^2 + q^2 + usq) / (2 * p * q)
    tmp = zero(ComplexF64)
    if cl == :ss
        if iszero(p) || iszero(q)
            tmp = (p^2 + q^2) / (p^2 + q^2 + usq)
        else
            tmp = -1 / (4 * p * q) * (-4 * p * q + (p^2 + q^2 - 2 * p * q * ξ) *(log(1-ξ)-log(-1-ξ)))
        end
    elseif cl == :dd
        if iszero(p) || iszero(q)
            tmp = 0
        else
            tmp = -1 / (8 * p * q) * (6 * ξ * (p^2 + q^2 - 2 * p * q * ξ) + (p^2 + q^2 - 2 * p * q * ξ) * (-1 + 3 * ξ^2) *(log(1-ξ)-log(-1-ξ)))
        end
    end
    return tmp
end

function pwa_contact(p, q, usq; cl = :ss)
    tmp = 0.0
    if cl == :ss
        tmp = 1.0
    end
    return tmp
end

end