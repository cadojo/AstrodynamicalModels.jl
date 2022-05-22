[![Tests](https://github.com/cadojo/AstrodynamicalModels.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/AstrodynamicalModels.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/AstrodynamicalModels.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/AstrodynamicalModels.jl/dev)


# AstrodynamicalModels.jl
_Solar system body physical attributes, dynamical models, and ephemeris interpolation made simple!_

## Physical Attributes
_Fetch physical attributes for major and minor bodies in our solar system!_

### Currently Implemented

* NAIF ID to solar system body name conversion, and vice versa.

### Soon to be Implemented

* Body mass parameter querying
* Body axial radii querying

## Dynamical Models
_Requires [`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl)!_

### Currently Implemented

Note – for all models below, you can optionally append state transition matrix dynamics.

* Restricted Two-body Problem
* Circular Restricted Three-body Problem
* N-body Problem

### Soon to be Implemented

* Aspherical Restricted Two-body Problem

### Way in the Future

* Solar radiation pressure dynamics
* _Others? Let me know about, or submit a PR with, your desired astrodynamics models!_

## Ephemeris Interpolation
_Requires [`HORIZONS.jl`](https://github.com/PerezHz/HORIZONS.jl)!_

### Currently Implemented

* Ephemeris fetching and formatting through a `DataFrame` 
* Cubic Hermite Interpolation of fetched ephemeris data through `Interpolations`

### Soon to be Implemented

* Model building at any (terrestrial) epoch

### Way in the Future

* Model building at any epoch, using [`AstroTime`](https://github.com/JuliaAstro/AstroTime.jl)

## License & Usage

The [license file](LICENSE) applies to all code in this package, but note 
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
