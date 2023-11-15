# ----------------------------------------------------------DD*----------------------------------------------------#
function potential_DnDsc(E, p, q; I=0, Cv=0.0, Cp=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDn, mDsc), q0f(E, mDn, mDsc, mDsc, mDn)
    res = zero(ComplexF64)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπc^2 - q02^2) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(p, q, mρc^2 - q02^2; a=Cv) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * (-2)
    return res
end

function potential_DcDsn(E, p, q; I=0, Cv=0.0, Cp=0.0, isπ=true)
    q01, q02 = q0f(E, mDc, mDsn, mDc, mDsn), q0f(E, mDc, mDsn, mDsn, mDc)
    res = zero(ComplexF64)

    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    if isπ == true
        res += -4 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπc^2 - q02^2) * (I - 1 / 2) * 2
    end
    res += ((2 * sqrt(2) * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(p, q, mρc^2 - q02^2; a=Cv) * (I - 1 / 2) * 2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * (-2)
    return res
end

function potential_DDscross(E, p, q; I=0, Cv=0.0, Cp=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsc, mDc, mDsn), q0f(E, mDn, mDsc, mDsn, mDc)
    res = zero(ComplexF64)

    res += 2 * β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q01^2)
    if isπ == true
        res += 2 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss_deform(p, q, mπn^2 - q02^2) * (I - 1 / 2) * 2
    end
    res += -2 * g^2 / (3fπ^2) / sqrt(16) * pwa_tensor_ss(p, q, mη^2 - q02^2; a=Cp) * (I - 1 / 2) * 2
    res += (-(2 * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(p, q, mρn^2 - q02^2; a=Cv) * (I - 1 / 2) * 2)
    res += ((2 * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(p, q, mω^2 - q02^2; a=Cv) * (I - 1 / 2) * 2)
    #ηc,J/ψ
    # res += 2 * g^2 / fπ^2 / sqrt(16) * pwa_tensor_ss(p, q, mηc^2 - q02^2; a=a) * (-2) * (I - 1 / 2) * 2
    # res += (-(2 * gV * gλ)^2 / sqrt(16) * pwa_vvv_ss(p, q, mJψ^2 - q02^2; a=a) * (I - 1 / 2) * 2 * (-2))
    return res
end



#----------------------------------------------------------DD*bar----------------------------------------------------#
function potential_Xn(E, p, q; C=1, Cv=0.0, Cp=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsn, mDn, mDsn), q0f(E, mDn, mDsn, mDsn, mDn)
    res = zero(ComplexF64)

    if isπ == true
        res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss_deform(p, q, mπn^2 - q02^2)
    end
    res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(p, q, mη^2 - q02^2; a=Cp)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * pwa_central(p, q, mρn^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(p, q, mρn^2 - q02^2; a=Cv)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * pwa_central(p, q, mω^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(p, q, mω^2 - q02^2; a=Cv)
    res += -4 * gs^2 / 4 * pwa_central(p, q, mσ^2 - q01^2)
    #J/ψ
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / 4 * pwa_central(p, q, mJψ^2 - q01^2) * 2
    # res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(p, q, mJψ^2 - q02^2; a=a) * 2
    # #ηc
    # res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(p, q, mηc^2 - q02^2; a=a) * 2
    return res
end

function potential_Xc(E, p, q; C=1, Cv=0.0, Cp=0.0, isπ=true)
    q01, q02 = q0f(E, mDc, mDsc, mDc, mDsc), q0f(E, mDc, mDsc, mDsc, mDc)
    res = zero(ComplexF64)

    if isπ == true
        res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss_deform(p, q, mπn^2 - q02^2)
    end
    res += (-C) * 1 / 3 * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(p, q, mη^2 - q02^2; a=Cp)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / 4 * pwa_central(p, q, mρn^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(p, q, mρn^2 - q02^2; a=Cv)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / 4 * pwa_central(p, q, mω^2 - q01^2)
    res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(p, q, mω^2 - q02^2; a=Cv)
    res += -4 * gs^2 / 4 * pwa_central(p, q, mσ^2 - q01^2)
    #J/ψ
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / 4 * pwa_central(p, q, mJψ^2 - q01^2) * 2
    # res += (-C) * (2 * gλ * gV)^2 / 4 * pwa_vvv_ss(p, q, mJψ^2 - q02^2; a=a) * 2
    # #ηc
    # res += (-C) * (√2 * g / fπ)^2 / 4 * pwa_tensor_ss(p, q, mηc^2 - q02^2; a=a) * 2

    return res
end

function potential_Xcross(E, p, q; C=1, Cv=0.0, Cp=0.0, isπ=true)
    q01, q02 = q0f(E, mDn, mDsn, mDc, mDsc), q0f(E, mDn, mDsn, mDsc, mDc)
    res = zero(ComplexF64)

    res += -(√2 * gV * β)^2 * (1 - q01^2 / mρc^2) / 4 * pwa_central(p, q, mρc^2 - q01^2)
    res += (-C) * (2 * √2 * gλ * gV)^2 / 4 * pwa_vvv_ss(p, q, mρc^2 - q02^2; a=Cv)
    if isπ == true
        res += (-C) * (2 * g / fπ)^2 / 4 * pwa_tensor_ss(p, q, mπc^2 - q02^2)
    end
    return res
end




#----------------------------------------------------------D*D*----------------------------------------------------#
function potential_DsDs(E, p, q; I=0, J=0, Cv=0.0, Cp=0.0, isπ=true)
    q01, q02 = q0f(E, mDsc, mDsn, mDsc, mDsn), q0f(E, mDsc, mDsn, mDsn, mDsc)
    res = zero(ComplexF64)

    res += -gs^2 * pwa_central(p, q, mσ^2 - q01^2)
    if isπ == true
        res += -g^2 / fπ^2 / 2 * pwa_tensor_ss(p, q, mπn^2 - q01^2)
        res += -g^2 / fπ^2 / 2 * pwa_tensor_ss(p, q, mπc^2 - q02^2) * (-2) * (I - 1 / 2) * 2
    end
    res += -g^2 / fπ^2 / 2 * pwa_tensor_ss(p, q, mη^2 - q01^2; a=Cp) * (-1 / 3)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2) + 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mρn^2 - q01^2; a=Cv, J=J)
    res += (-β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2) + 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mω^2 - q01^2; a=Cv, J=J)) * (-1)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * (-2)
    # res += -g^2 / fπ^2 / 2 * pwa_tensor_ss(p, q, mηc^2 - q01^2) * (-2)

    res += (-β^2 * gV^2 * (1 - q02^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q02^2) + 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mρc^2 - q01^2; a=Cv, J=J)) * (-2) * (I - 1 / 2) * 2
    return res
end



#----------------------------------------------------------meson-antimeson----------------------------------------------------#
function potential_DDbarn(E, p, q)
    q01 = q0f(E, mDn, mDn, mDn, mDn)
    res = zero(ComplexF64)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * 2
    return res
end
function potential_DDbarc(E, p, q)
    q01 = q0f(E, mDc, mDc, mDc, mDc)
    res = zero(ComplexF64)
    res += -4 * gs^2 / sqrt(16) * pwa_central(p, q, mσ^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * 2
    return res
end
function potential_DDbarcross(E, p, q)
    q01 = q0f(E, mDn, mDn, mDc, mDc)
    res = zero(ComplexF64)
    res += -β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q01^2) * 2
    return res
end




function potential_DsDsbarn(E, p, q; J=0, Cv=0.0, Cp=0.0)
    q01 = q0f(E, mDsn, mDsn, mDsn, mDsn)
    res = zero(ComplexF64)

    res += -gs^2 * pwa_central(p, q, mσ^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mρn^2 - q01^2; a=Cv, J=J)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mω^2 - q01^2; a=Cv, J=J)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * 2
    res += g^2 / fπ^2 / 2 * pwa_Tensor_ss(p, q, mπn^2 - q01^2; J=J)
    res += g^2 / fπ^2 / 2 * pwa_Tensor_ss(p, q, mη^2 - q01^2; a=Cp, J=J) * 1 / 3
    return res
end

function potential_DsDsbarc(E, p, q; J=0, Cv=0.0, Cp=0.0)
    q01 = q0f(E, mDsc, mDsc, mDsc, mDsc)
    res = zero(ComplexF64)

    res += -gs^2 * pwa_central(p, q, mσ^2 - q01^2)
    res += -β^2 * gV^2 * (1 - q01^2 / mρn^2) / sqrt(16) * pwa_central(p, q, mρn^2 - q01^2)
    res += 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mρn^2 - q01^2; a=Cv, J=J)
    res += -β^2 * gV^2 * (1 - q01^2 / mω^2) / sqrt(16) * pwa_central(p, q, mω^2 - q01^2)
    res += 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mω^2 - q01^2; a=Cv, J=J)
    res += -β^2 * gV^2 * (1 - q01^2 / mJψ^2) / sqrt(16) * pwa_central(p, q, mJψ^2 - q01^2) * 2
    res += g^2 / fπ^2 / 2 * pwa_Tensor_ss(p, q, mπn^2 - q01^2; J=J)
    res += g^2 / fπ^2 / 2 * pwa_Tensor_ss(p, q, mη^2 - q01^2; a=Cp, J=J) * 1 / 3
    return res
end

function potential_DsDsbarcross(E, p, q; J=0, Cv=0.0, Cp=0.0)
    q01 = q0f(E, mDsn, mDsn, mDsc, mDsc)
    res = zero(ComplexF64)

    res += -β^2 * gV^2 * (1 - q01^2 / mρc^2) / sqrt(16) * pwa_central(p, q, mρc^2 - q01^2) * 2
    res += 2 * gλ^2 * gV^2 * pwa_VVV_ss(p, q, mρc^2 - q01^2; a=Cv, J=J) * 2
    res += g^2 / fπ^2 * pwa_Tensor_ss(p, q, mπc^2 - q01^2; J=J)
    return res
end