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
