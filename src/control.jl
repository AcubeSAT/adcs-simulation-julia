@concrete struct PDController
    Kp
    Kd
end

# TODO: quaternion frame?
function calculate_torque(PD::PDController, qtarget, qestimated, w, wtarget)
    qrel = qtarget * conj(qestimated) # TODO: should it be conj(qtarget)?
    werr = w - wtarget
    return -sign(real(qrel)) * PD.Kp * vec(qrel) - PD.Kd * werr
end

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
