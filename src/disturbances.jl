residual_dipole(m, b) = cross(m, b)

function epoch_to_datetime(epcs)
    function x(v)
        v[end] = floor(v[end] / 1.0e6)
        return v
    end
    return [DateTime(x([caldate(epc)...])...) for epc in epcs]
end

function gravity_gradient_tensor(model, p, t, max_degree, P, dP; δ=1e-3)
    T = @MMatrix zeros(3, 3)
    for i in 1:3
        g_forward = CartesianFromSpherical()(Spherical(GravityModels.gravitational_field_derivative(model,
            [p[1] + (i == 1 ? δ : 0), p[2] + (i == 2 ? δ : 0), p[3] + (i == 3 ? δ : 0)],
            t, max_degree=max_degree, P=P, dP=dP)...))
        g_backward = CartesianFromSpherical()(Spherical(GravityModels.gravitational_field_derivative(model,
            [p[1] - (i == 1 ? δ : 0), p[2] - (i == 2 ? δ : 0), p[3] - (i == 3 ? δ : 0)],
            t, max_degree=max_degree, P=P, dP=dP)...))

        @. T[:, i] = (g_forward - g_backward) / (2δ)
    end
    return T
end

function gravity_torque(G_ecef, R_ecef_to_body, I)
    G_body = transpose(R_ecef_to_body) * G_ecef * R_ecef_to_body
    T1 = G_body[2,3] * (I[3,3] - I[2,2]) + G_body[1,3] * I[1,2] - G_body[1,2] * I[1,3] + I[2,3] * (G_body[3,3] - G_body[2,2])
    T2 = G_body[1,3] * (I[1,1] - I[3,3]) - G_body[2,3] * I[1,2] + G_body[1,2] * I[2,3] + I[1,3] * (G_body[1,1] - G_body[3,3])
    T3 = G_body[1,2] * (I[2,2] - I[1,1]) + G_body[2,3] * I[1,3] - G_body[1,3] * I[2,3] + I[1,2] * (G_body[2,2] - G_body[1,1])
    return SVector(3*T1, 3*T2, 3*T3)
end
