@kwdef @concrete struct ReactionWheel
    J
    w
    wdeadzone = 52.359877
    wmid = wdeadzone + 10
    α
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
    return RW.maxtorque / (1 + exp(-RW.α * (RW.w - RW.wmid)))
end
