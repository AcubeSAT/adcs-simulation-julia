@kwdef @concrete struct ReactionWheel
    J
    w
    wmax = 628.318530
    wmaxmid = wmax - 10
    wdeadzone = 52.359877
    wdeadzonemid = wdeadzone + 10
    saturationα
    deadzoneα
    maxtorque = 0.001
    Fc = 0.001
    b = 0.0005
    Fs = 0.0005
    w0 = 1
end

function stribeck(RW::ReactionWheel)
    τcoulomb = RW.Fc * sign(RW.w)
    τviscous = RW.b * RW.w
    τstribeck = RW.Fs * exp(-(RW.w / RW.w0)^2) * sign(RW.w)
    return τcoulomb + τviscous + τstribeck
end

function deadzone_compensation(RW::ReactionWheel)
    abs(RW.w) <= RW.wdeadzone || return zero(RW.w)
    return -sign(RW.w) * RW.maxtorque / (1 + exp(-RW.deadzoneα * (RW.w - RW.wdeadzonemid)))
end

function saturation_compensation(RW::ReactionWheel)
    RW.wmax <= abs(RW.w) || return zero(RW.w)
    return -sign(RW.w) * RW.maxtorque / (1 + exp(-RW.saturationα * (RW.w - RW.wmaxmid)))
end
