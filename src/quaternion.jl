# Inspired by @moble and Base.Complex

# Signify standard arithmetic operations should be implemented on it
abstract type AbstractQuaternion{T<:Number} <: Number end

struct Quaternion{T<:Number} <: AbstractQuaternion{T}
    coeffs::SVector{4,T}
    Quaternion{T}(v::SVector{4,T}) where {T<:Number} = new{T}(v)
    # Always use SVector internally
    Quaternion{T}(v::A) where {T<:Number, A<:AbstractVector} = new{T}(SVector{4,T}(v))
end

# Given a numeric type, return a QT type (not instance) specialized on T
(::Type{QT})(::Type{T}) where {T<:Number, QT<:AbstractQuaternion} = QT{T}
# Given a QT specialized on T, return another QT also specialized on T
(::Type{QT})(::Type{<:AbstractQuaternion{T}}) where {T<:Number, QT<:AbstractQuaternion} = QT{T}

# Construct without having to specify T
Quaternion(v::SVector{4,T}) where {T<:Number} = Quaternion{T}(v)
Quaternion(v::AbstractVector{T}) where {T<:Number} = Quaternion{T}(v)

# Explicitly state type, use SVector promote_rules for mixed types
(::Type{QT})(w, x, y, z) where {T<:Number, QT<:AbstractQuaternion{T}} = QT(SVector{4,T}(w, x, y, z))
# Rely on type inference, use SVector promote_rules for mixed types
function (::Type{QT})(w, x, y, z) where {QT<:AbstractQuaternion}
    v = SVector{4}(w, x, y, z)
    return QT{eltype(v)}(v)
end

coeffs(q::AbstractQuaternion) = getfield(q, :coeffs)

# Type-preserving copy constructor
(::Type{QT})(q::QT) where {T<:Number, QT<:AbstractQuaternion{T}} = QT(coeffs(q)...)
# Type conversion copy constructor
(::Type{QT})(q::AbstractQuaternion{S}) where {T<:Number, S<:Number, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(coeffs(q)...))

const QuaternionF64 = Quaternion{Float64}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF16 = Quaternion{Float16}
