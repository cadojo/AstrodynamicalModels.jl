"""
Provides common astrodynamical models through `ModelingToolkit`,
physical attributes for solar system bodies through `SPICE`,
and ephemeris model downloading and interpolation through
`HORIZONS` and `Interpolations`!

!!! tip
    Check out the `ModelingToolkit` docs to learn how to use these systems for orbit
    propagation with `DifferentialEquations`, or see `GeneralAstrodynamics` for some
    convenient orbit propagation wrappers.

# Usage

## Physical Attributes

```julia
naifcode("earth") == 399
naifname("mars") == 499

μₑ = massparameter("earth")
Rₘ = meanradius("mars")
```

## Dynamics

```julia
using ModelingToolkit

length(states(R2BP())) == 6
length(states(NBP(2))) == 12
```

## Ephemeris

!!! warning
    Before usage, see the warning about copyrighted tools in this project's README, or in
    the Extended Help section of this docstring!

```julia
using HORIZONS
ephemeris("earth") isa DataFrame

using Dates
ephemeris("earth", (now()-Year(5), now()+Year(5)), Day(1)) isa DataFrame

using Interpolations: CubicHermite
ContinuousEphemeris("earth") isa ContinuousEphemeris{<:Number, <:Number, <:CubicHermite, <:CubicHermite, <:CubicHermite}
```

# Extended help

## License

!!! warning
    The license in this docstring applies to all code in this package, but note
    that the ephemeris data downloading functionality of this package
    (e.g. the `ephemeris` function) uses copyrighted tools with their
    own usage and code sharing restrictions. These copyrighted tools are
    owned by the Jet Propulsion Laboratory at the California Institute
    of Technology. These copyrighted tools are, in part, located at
    ftp://ssd.jpl.nasa.gov/pub/ssd/SCRIPTS/.

    Open source wrappers for these copyrighted tools exist, e.g.
    `HORIZONS.jl` and AstroPy's `astroquery`. The ephemeris handling
    functionality of this package uses `HORIZONS.jl`.

    It is the responsibility of the user to verify they meet the necessary
    requirements, as specifed by JPL in their scripts (e.g. `vec_tbl`),
    before they share or use the copyrighted tools. This includes the
    ephemeris handling functionality in this package (e.g. `ephemeris`).

$(LICENSE)

## Exports

### Always

* `naifcode`
* `naifname`
* `axialradii`
* `meanradius`
* `massparameter`

### Requires `ModelingToolkit`

* `NBP`
* `R2BP`
* `CR3BP`
* `R2BPFunction`
* `CR3BPFunction`
* `NBPFunction`

### Requires `HORIZONS`

* `ephemeris`

#### Requires `Interpolations` (and `HORIZONS`)

* `ContinuousEphemeris`

## Imports

### Always

* `Requires`, `SPICE`, `LinearAlgebra`, `DocStringExtensions`

### Requires `ModelingToolkit`

* `Memoize`

### Requires `HORIZONS`

* `CSV`, `DataFrames`

"""
module AstrodynamicalModels

import SPICE

using Requires
using LinearAlgebra
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)

    $(DOCSTRING)
    """

@template (TYPES, CONSTANTS) =
    """
    $(TYPEDEF)

    $(DOCSTRING)
    """

function __init__()
    @require ModelingToolkit="961ee093-0014-501f-94e3-6117800e7a78" begin

        using Memoize: @memoize

        export R2BP, CR3BP, NBP
        export R2BPFunction, CR3BPFunction, NBPFunction

        include(joinpath(@__DIR__, "models", "R2BP.jl"))
        include(joinpath(@__DIR__, "models", "CR3BP.jl"))
        include(joinpath(@__DIR__, "models", "NBP.jl"))
    end

    @require HORIZONS="5a3ac768-beb4-554a-9c98-3342fe3377f5" begin

        using Dates
        using DataFrames
        using StaticArraysCore

        include(joinpath("ephemeris", "ephemeris.jl"))
        @require Interpolations="a98d9a8b-a2ab-59e6-89dd-64a1c18fca59" include(joinpath("ephemeris", "interpolations.jl"))
    end
end

include(joinpath("properties", "NAIF.jl"))

end # module
