@kwdef @concrete struct ReactionWheel
    J
    w
    wmax = 628.318530
    wmaxmid = wmax - 10
    wdeadzone = 52.359877
    wdeadzonemid = wdeadzone + 10
    saturationα
    deadzoneα
    maxtorque = 0.007
    Fc = 0.001
    b = 0.0005
    Fs = 0.0005
    w0 = 1
end

function stribeck(Fc, w, b, w0)
    τcoulomb = Fc * sign(w)
    τviscous = b * w
    τstribeck = Fc * exp(-(w / w0)^2) * sign(w)
    return τcoulomb + τviscous + τstribeck
end

function stribeck(RW::ReactionWheel)
    return stribeck.(RW.Fc, RW.w, RW.b, RW.w0)
end

function deadzone_compensation(w, wdeadzone, wdeadzonemid, deadzoneα, maxtorque)
    abs(w) <= wdeadzone || return zero(w)
    return -sign(w) * maxtorque / (1 + exp(-deadzoneα * (w - wdeadzonemid)))
end

function deadzone_compensation(RW::ReactionWheel)
    return deadzone_compensation.(RW.w, RW.wdeadzone, RW.wdeadzonemid, RW.deadzoneα, RW.maxtorque)
end

function saturation_compensation(w, wmax, saturationα, wmaxmid, maxtorque)
    wmax <= abs(w) || return zero(w)
    return -sign(w) * maxtorque / (1 + exp(-saturationα * (w - wmaxmid)))
end

function saturation_compensation(RW::ReactionWheel)
    return saturation_compensation.(RW.w, RW.wmax, RW.saturationα, RW.wmaxmid, RW.maxtorque)
end
