#
# Functions, and types for downloading and interpolating 
# solar system ephemeric data.
#
# This functionality is enabled by JPL HORIZONS and 
# HORIZONS.jl!
#

export ephemeris

"""
Given a celestial body's name, or the NAIF ID of the body,
return a DataFrame of ephemeris data for the specified 
`timespan` and `intervol` between time steps.

!!! warning
    The [license file](../../LICENSE) for this package applies to all code 
    in this package, but note 
    that the ephemeris data handling functionality of this package
    (e.g. the `ephemeris` function) uses copyrighted tools with their 
    own usage and code sharing restrictions. These copyrighted tools are 
    owned by the Jet Propulsion Laboratory at the California Institute 
    of Technology. These copyrighted tools are, in part, located at 
    ftp://ssd.jpl.nasa.gov/pub/ssd/SCRIPTS/. 

    Open source wrappers for these copyrighted tools exist, e.g. 
    `HORIZONS.jl` and AstroPy's `astroquery`. The ephemeris handling 
    functionality of this package uses `HORIZONS.jl`.

    It is the responsibility of the user to verify they meet the necessary
    requirements, as specifed by JPL in their copyrighted tools (e.g. `vec_tbl`), 
    before they share or use those tools.
"""
function ephemeris(
        body::Union{<:AbstractString, <:Integer}, 
        timespan::Tuple{<:Dates.AbstractDateTime, <:Dates.AbstractDateTime} = (now() - Year(5), now() + Year(5)),
        intervol::Dates.Period = Day(1);
        continuous  = Val{false},
        type        = Float64,
        email       = "", 
        wrt         = "@ssb", 
        epoch       = "J2000", 
    )

    if continuous == Val{true}
        if "Interpolations" ∉ map(string, Base.loaded_modules_array())
            throw(ArgumentError("The `continuous` keyword argument requires Interpolations.jl to be loaded!"))
        end
    elseif continuous != Val{false}
        throw(ArgumentError("The `continuous` keyword argument must be set to `Val{true}` or `Val{false}`!"))
    end

    if isnothing(email)
        EMAIL_ADDR = ""
    elseif isempty(email) 
        @warn """No email was provided throuh the `email` keyword argument! The JPL HORIZONS servers accept an email address with each request. Provide an email address using `email = "my-email@provider.com"`, or set `email = nothing` to suppress this warning.""" 
        EMAIL_ADDR = ""
    else
        EMAIL_ADDR = email
    end

    ID = body isa AbstractString ? naifcode(body) : body

    options = (
        EMAIL_ADDR = EMAIL_ADDR, 
        CENTER = wrt,
        REF_SYSTEM = epoch,
        REF_PLANE = "FRAME", 
        CSV_FORMAT = true, 
        VEC_TABLE = 2, 
        VEC_CORR = 1, 
        OUT_UNITS = 1, 
        VEC_LABELS = false, 
        VEC_DELTA_T = false, 
    )

    data, = HORIZONS.vec_tbl_csv(ID, timespan[1], timespan[2], intervol; options...)

    if data[1,:] == ["JDTDB", "Calendar_Date_TDB", "X", "Y", "Z", "VX", "VY", "VZ"]
        # Units in Julian Days since J2000, km, km/s
        ephemeris = DataFrame(
            :dⱼ => type.(data[2:end,1]),
            :x  => type.(data[2:end,3]),
            :y  => type.(data[2:end,4]),
            :z  => type.(data[2:end,5]),
            :ẋ  => type.(data[2:end,6]),
            :ẏ  => type.(data[2:end,7]),
            :ż  => type.(data[2:end,8]),
        )
    else
        @warn "Invalid columns! $(data[1,:])"
        throw(ErrorException("Ephemeris data downloaded successfully, but the column format is different than expected!"))
    end
    
    if continuous == Val{true}
        return ContinuousEphemeris(ephemeris)
    else
        return ephemeris
    end
end

"""
Returns the Cartesian state vector of the requested body at the provided date and time!
All `kwargs` are passed directly to the `timespan` method for `ephemeris`.

!!! warning
    The [license file](../../LICENSE) for this package applies to all code 
    in this package, but note 
    that the ephemeris data handling functionality of this package
    (e.g. the `ephemeris` function) uses copyrighted tools with their 
    own usage and code sharing restrictions. These copyrighted tools are 
    owned by the Jet Propulsion Laboratory at the California Institute 
    of Technology. These copyrighted tools are, in part, located at 
    ftp://ssd.jpl.nasa.gov/pub/ssd/SCRIPTS/. 

    Open source wrappers for these copyrighted tools exist, e.g. 
    `HORIZONS.jl` and AstroPy's `astroquery`. The ephemeris handling 
    functionality of this package uses `HORIZONS.jl`.

    It is the responsibility of the user to verify they meet the necessary
    requirements, as specifed by JPL in their copyrighted tools (e.g. `vec_tbl`), 
    before they share or use those tools.
"""
function ephemeris(body::Union{<:Integer, <:AbstractString}, time::Dates.AbstractDateTime; kwargs...)
    override = (; continuous = Val{false})
    options  = merge(kwargs, override)
    data     = ephemeris(body, (time, time + Year(1)), Year(1); options...)
    return DataFrame(data[1,2:end])
end

