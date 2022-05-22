#
# Functions to convert between body names and NAIF ID codes,
# and vice versa. This functionality is enabled by the 
# CSPICE library, CSPICE_jll, and SPICE.jl!
#
# Mass parameter and radii values are also provided through
# functions. This functionality is enabled by JPL HORIZONS
# and the IAU!
#

export naifcode, naifname

"""
Return the NAIF ID for the provided celestial body name. If no ID for 
the provided body is found, a `KeyError` exception is thrown.

!!! note
    This is a simple wrapper around `SPICE.bodn2c`, which is itself 
    a wrapper around the CSPICE library's `bod2nc` function! This 
    integer code lookup is robust to different capitalizations. 
    For example, `"Earth"`, `"earth"`, and `"eArTH"` should all return 
    integer code `399`.
"""
function naifcode(name::AbstractString)
    code = SPICE.bodn2c(name)

    if code isa Integer
        return code
    else
        throw(KeyError("No NAIF ID for $name could be found!"))
    end
end

"""
Return the celestial body name for the provided NANIF ID. If no name for 
the provided ID is found, a `KeyError` exception is thrown.

!!! note
    This is a simple wrapper around `SPICE.bodc2n`, which is itself 
    a wrapper around the CSPICE library's `bodc2n` function!
"""
function naifname(code::Int)
    name = SPICE.bodc2n(code)

    if name isa AbstractString
        return name
    else
        throw(KeyError("No body name for $code could be found!"))
    end
end