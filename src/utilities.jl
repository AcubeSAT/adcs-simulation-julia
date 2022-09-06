using Quaternions

function rotateVector(q::Quaternions.Quaternion, v::Array{Float64})
    rotatedQuaternion = conj(q) * Quaternions.Quaternion(0, v[1], v[2], v[3]) * q
    rotateVector = [rotatedQuaternion.v1, rotatedQuaternion.v2, rotatedQuaternion.v3]
end

function DCM2Quat(M::Matrix{Float64})
    qw = √(max(0, 1 + M[1, 1] + M[2, 2] + M[3, 3])) / 2
    qx = √(max(0, 1 + M[1, 1] - M[2, 2] - M[3, 3])) / 2
    qy = √(max(0, 1 - M[1, 1] + M[2, 2] - M[3, 3])) / 2
    qz = √(max(0, 1 - M[1, 1] - M[2, 2] + M[3, 3])) / 2
    normalize(Quaternions.Quaternion(qw, qx, qy, qz))
end

function julianDay2GST(julianDay::Float64)

    deg2rad = pi / 180.0
    tut1 = (julianDay - 2451545.0) / 36525.0
    temp =
        -6.2e-6 * tut1 * tut1 * tut1 +
        0.093104 * tut1 * tut1 +
        (876600.0 * 3600.0 + 8640184.812866) * tut1 +
        67310.54841
    # 360/86400 = 1/240, to deg, to rad
    temp = rem(temp * deg2rad / 240.0, 2π)

    if (temp < 0.0)
        temp = temp + twopi
    end

    timeGST = temp
end
