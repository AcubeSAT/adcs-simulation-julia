@concrete struct PDController
    Kp
    Kd
end

function calculate_torque(PD::PDController, qtarget, qestimated, w, wtarget)
    qrel = qestimated * conj(qtarget)
    werr = w - wtarget
    return -sign(real(qrel)) * PD.Kp * vec(qrel) - PD.Kd * werr
end

# FIXME: This is for a 3RW system! It is not applicable for AcubeSAT
function decompose_torque(τ, b, msaturation)
    τw = dot(τ, b) / dot(b, b) * b
    τm = τ - τw
    m = cross(b, τm) / dot(b, b)
    mprime = m / norm(τm)
    ksm = min(abs(msaturation / mprime[1]),
        abs(msaturation / mprime[2]),
        abs(msaturation / mprime[3]))
    mtrue = ksm * mprime
    τsm = ksm * cross(b, mprime)
    τw = τ - τsm
    return τw, τsm, mtrue
end
