#
# Functions and types to interpolate ephemeris data with
# a Cubic Hermite interpolator.
#
# This functionality is ~entirely~ enabled by Interpolations.jl!
#

export ContinuousEphemeris

"""
An abstract supertype for all possible interpolations for ephemeris data.
"""
abstract type AbstractEphemerisInterpolation end

"""
A callable structure which holds ephemeris data for a body's positions and velocities between the timepoints specified by `T`.
Call the object with a Julian day input to interpolate the body's positions and velocities using `Interpolations.CubicHermite`.
"""
struct ContinuousEphemeris{
            D<:Number,
            TX<:Interpolations.CubicHermite,
            TY<:Interpolations.CubicHermite,
            TZ<:Interpolations.CubicHermite
        } <:AbstractEphemerisInterpolation
    T::Tuple{D,D}
    X::TX
    Y::TY
    Z::TZ
end

"""
Construct a `ContinuousEphemeris` instance from an appropriately formatted
`DataFrame`. This is mostly an internal package method. See `ephemeris(args...; continuous = Val{true}, kwargs...)`
for a more user friendly constructor!
"""
function ContinuousEphemeris(data::DataFrame)

    T = eltype(data.dⱼ)
    timespan = (
        data.dⱼ[1],
        data.dⱼ[end],
    )

    X = Interpolations.CubicHermite(data.dⱼ, data.x, data.ẋ)
    Y = Interpolations.CubicHermite(data.dⱼ, data.y, data.ẏ)
    Z = Interpolations.CubicHermite(data.dⱼ, data.z, data.ż)

    return ContinuousEphemeris{T, typeof(X), typeof(Y), typeof(Z)}(
        timespan, X, Y, Z
    )
end

"""
Interpolate the body's position and velocity by a Julian day input.
"""
function (ephemeris::ContinuousEphemeris)(t::Number; bias = false)
    time = bias ? t + ephemeris.T[1] : t

    return SVector{6}(
        ephemeris.X(time),
        ephemeris.Y(time),
        ephemeris.Z(time),
        Interpolations.gradient(ephemeris.X, time),
        Interpolations.gradient(ephemeris.Y, time),
        Interpolations.gradient(ephemeris.Z, time),
    )

end

function Base.show(io::IO, ephemeris::ContinuousEphemeris)
    print(io, "Cubic Hermite interpolation of ephemeris data")
end
