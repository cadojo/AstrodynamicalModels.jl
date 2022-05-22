#
# Functions and types to interpolate ephemeris data with
# a Cubic Hermite interpolator.
# 
# This functionality is ~entirely~ enabled by Interpolations.jl!
#

export ContinuousEphemeris

abstract type AbstractEphemerisInterpolation end

struct ContinuousEphemeris{T1<:Number, T2<:Number, TX<:Interpolations.CubicHermite, TY<:Interpolations.CubicHermite, TZ<:Interpolations.CubicHermite} <:AbstractEphemerisInterpolation
    T::Tuple{T1,T2}
    X::TX
    Y::TY
    Z::TZ
end

function ContinuousEphemeris(data::DataFrame)

    timespan = (
        data.dⱼ[1],
        data.dⱼ[end],
    )

    X = Interpolations.CubicHermite(data.dⱼ, data.x, data.ẋ)
    Y = Interpolations.CubicHermite(data.dⱼ, data.y, data.ẏ)
    Z = Interpolations.CubicHermite(data.dⱼ, data.z, data.ż)

    return ContinuousEphemeris{typeof(timespan[1]), typeof(timespan[2]), typeof(X), typeof(Y), typeof(Z)}(
        timespan, X, Y, Z
    )
end

function (ephemeris::ContinuousEphemeris{F})(timepoint::Number; type = SVector{6,F}) where F <: AbstractFloat
    t = timepoint isa F ? timepoint : F(timepoint)
    return (
        ephemeris.X(t),
        ephemeris.Y(t),
        ephemeris.Z(t),
        Interpolations.gradient(ephemeris.X, t),
        Interpolations.gradient(ephemeris.Y, t),
        Interpolations.gradient(ephemeris.Z, t),
    ) |> type
end

function Base.show(io::IO, ephemeris::ContinuousEphemeris{F}) where F <: AbstractFloat
    print(io, "Cubic Hermite interpolation of ephemeris data with eltype $F")
end
