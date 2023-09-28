@concrete struct PDController
    Kp
    Kd
end

function calculate_torque(PD::PDController, qtarget, qestimated, w, wtarget)
    qrel = qestimated * conj(qtarget)
    werr = w - wtarget
    return -sign(real(qrel)) * PD.Kp * vec(qrel) - PD.Kd * werr
end
