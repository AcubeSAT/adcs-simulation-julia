using Quaternions
using LinearAlgebra
using Infiltrator

function simulateStep(parameterStruct, ω::Vector, q::Quaternion, T::Vector, β::Vector)
    Δt = parameterStruct.Δt
    σu = parameterStruct.σu
    In = parameterStruct.In
    q = q * exp(Δt * Quaternion(ω) / 2)
    Δω = In \ (T - ω × (In * ω))
    ω = ω + Δt * Δω
    β = β + σu * sqrt(Δt) * randn(3)
    (q, ω, β)
end
