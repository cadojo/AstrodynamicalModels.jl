var documenterSearchIndex = {"docs":
[{"location":"docstrings/#Documentation","page":"Docstrings","title":"Documentation","text":"","category":"section"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"All docstrings!","category":"page"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"Modules = [\n    AstrodynamicalModels\n]\nOrder = [:module, :type, :function]","category":"page"},{"location":"docstrings/#AstrodynamicalModels.AstrodynamicalModels","page":"Docstrings","title":"AstrodynamicalModels.AstrodynamicalModels","text":"Provides astrodynamical models as AstrodynamicalModels.ODESystems. Check out the ModelingToolkit docs to learn how to use these systems for orbit propagation with DifferentialEquations, or see GeneralAstrodynamics for some convenient orbit propagation wrappers.\n\nExtended help\n\nLicense\n\nMIT License\n\nCopyright (c) 2023 Joseph D Carpinelli\n\nPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n\nExports\n\nAttitudeFunction\nAttitudeParameters\nAttitudeState\nAttitudeSystem\nCR3BFunction\nCR3BOrbit\nCR3BParameters\nCR3BState\nCR3BSystem\nCartesianOrbit\nCartesianState\nKeplerianOrbit\nKeplerianParameters\nKeplerianState\nNBFunction\nNBSystem\nOrbit\nOrbitalElements\nPlanarEntryFunction\nPlanarEntryParameters\nPlanarEntryState\nPlanarEntrySystem\nR2BFunction\nR2BOrbit\nR2BParameters\nR2BState\nR2BSystem\ndynamics\nparameters\nstate\nsystem\n\nImports\n\nBase\nCore\nDocStringExtensions\nLinearAlgebra\nMemoize\nModelingToolkit\nSciMLBase\nStaticArrays\nSymbolics\n\n\n\n\n\n","category":"module"},{"location":"docstrings/#AstrodynamicalModels.AstrodynamicalOrbit","page":"Docstrings","title":"AstrodynamicalModels.AstrodynamicalOrbit","text":"abstract type AstrodynamicalOrbit{U, P}\n\nAn abstract supertype for all orbits. \n\nExtended Help\n\nTo support the AstrodynamicalOrbit interface, you must implement the following methods.\n\nAstrodynamicalModels.states(orbit)\nAstrodynamicalModels.parameters(orbit)\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.AstrodynamicalParameters","page":"Docstrings","title":"AstrodynamicalModels.AstrodynamicalParameters","text":"abstract type AstrodynamicalParameters{F, N} <: StaticArraysCore.FieldVector{N, F}\n\nAn abstract supertype for all astrodynamical parameter vectors.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.AstrodynamicalState","page":"Docstrings","title":"AstrodynamicalModels.AstrodynamicalState","text":"abstract type AstrodynamicalState{F, N} <: StaticArraysCore.FieldVector{N, F}\n\nAn abstract supertype for all astrodynamical state vectors.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.AttitudeParameters","page":"Docstrings","title":"AstrodynamicalModels.AttitudeParameters","text":"struct AttitudeParameters{F} <: AstrodynamicalModels.AstrodynamicalParameters{F, 15}\n\nA parameter vector for attitude dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.AttitudeState","page":"Docstrings","title":"AstrodynamicalModels.AttitudeState","text":"mutable struct AttitudeState{F} <: AstrodynamicalModels.AstrodynamicalState{F, 7}\n\nA mutable state vector for attitude dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.CR3BOrbit","page":"Docstrings","title":"AstrodynamicalModels.CR3BOrbit","text":"struct Orbit{var\"#s19\"<:CartesianState, var\"#s18\"<:CR3BParameters} <: AstrodynamicalModels.AstrodynamicalOrbit{var\"#s19\"<:CartesianState, var\"#s18\"<:CR3BParameters}\n\nAn Orbit which exists within CR3BP dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.CR3BParameters","page":"Docstrings","title":"AstrodynamicalModels.CR3BParameters","text":"struct CR3BParameters{F} <: AstrodynamicalModels.AstrodynamicalParameters{F, 1}\n\nA paremeter vector for CR3BP dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.CR3BState","page":"Docstrings","title":"AstrodynamicalModels.CR3BState","text":"mutable struct CartesianState{F} <: AstrodynamicalModels.AstrodynamicalState{F, 6}\n\nCartesianState\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.CartesianOrbit","page":"Docstrings","title":"AstrodynamicalModels.CartesianOrbit","text":"struct Orbit{var\"#s19\"<:CartesianState, P<:(AbstractVector)} <: AstrodynamicalModels.AstrodynamicalOrbit{var\"#s19\"<:CartesianState, P<:(AbstractVector)}\n\nAn Orbit which exists within R2BP dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.CartesianSTM","page":"Docstrings","title":"AstrodynamicalModels.CartesianSTM","text":"mutable struct CartesianSTM{F} <: StaticArraysCore.FieldMatrix{6, 6, F}\n\nA mutable matrix, with labels, for a 6DOF Cartesian state transition matrix.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.CartesianState","page":"Docstrings","title":"AstrodynamicalModels.CartesianState","text":"mutable struct CartesianState{F} <: AstrodynamicalModels.AstrodynamicalState{F, 6}\n\nA mutable vector, with labels, for 6DOF Cartesian states.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.KeplerianOrbit","page":"Docstrings","title":"AstrodynamicalModels.KeplerianOrbit","text":"struct Orbit{var\"#s19\"<:OrbitalElements, var\"#s18\"<:KeplerianParameters} <: AstrodynamicalModels.AstrodynamicalOrbit{var\"#s19\"<:OrbitalElements, var\"#s18\"<:KeplerianParameters}\n\nAn Orbit which exists within Keplerian dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.KeplerianParameters","page":"Docstrings","title":"AstrodynamicalModels.KeplerianParameters","text":"struct KeplerianParameters{F} <: AstrodynamicalModels.AstrodynamicalParameters{F, 1}\n\nA parameter vector for Keplerian dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.KeplerianState","page":"Docstrings","title":"AstrodynamicalModels.KeplerianState","text":"mutable struct OrbitalElements{F} <: AstrodynamicalModels.AstrodynamicalState{F, 6}\n\nOrbitalElements\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.Orbit","page":"Docstrings","title":"AstrodynamicalModels.Orbit","text":"struct Orbit{U<:(AbstractVector), P<:(AbstractVector)} <: AstrodynamicalModels.AstrodynamicalOrbit{U<:(AbstractVector), P<:(AbstractVector)}\n\nA full representation of an orbit, including a numerical state, and the parameters of the system.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.OrbitalElements","page":"Docstrings","title":"AstrodynamicalModels.OrbitalElements","text":"mutable struct OrbitalElements{F} <: AstrodynamicalModels.AstrodynamicalState{F, 6}\n\nA mutable vector, with labels, for 6DOF Keplerian states.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.PlanarEntryParameters","page":"Docstrings","title":"AstrodynamicalModels.PlanarEntryParameters","text":"struct PlanarEntryParameters{F} <: AstrodynamicalModels.AstrodynamicalParameters{F, 7}\n\nA parameter vector for planar entry dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.PlanarEntryState","page":"Docstrings","title":"AstrodynamicalModels.PlanarEntryState","text":"mutable struct PlanarEntryState{F} <: AstrodynamicalModels.AstrodynamicalState{F, 4}\n\nA state vector for planar entry dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.R2BOrbit","page":"Docstrings","title":"AstrodynamicalModels.R2BOrbit","text":"struct Orbit{var\"#s19\"<:CartesianState, var\"#s18\"<:R2BParameters} <: AstrodynamicalModels.AstrodynamicalOrbit{var\"#s19\"<:CartesianState, var\"#s18\"<:R2BParameters}\n\nAn Orbit which exists within R2BP dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.R2BParameters","page":"Docstrings","title":"AstrodynamicalModels.R2BParameters","text":"struct R2BParameters{F} <: AstrodynamicalModels.AstrodynamicalParameters{F, 1}\n\nA parameter vector for R2BP dynamics.\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.R2BState","page":"Docstrings","title":"AstrodynamicalModels.R2BState","text":"mutable struct CartesianState{F} <: AstrodynamicalModels.AstrodynamicalState{F, 6}\n\nCartesianState\n\n\n\n\n\n","category":"type"},{"location":"docstrings/#AstrodynamicalModels.AttitudeFunction-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.AttitudeFunction","text":"AttitudeFunction(; stm, name, kwargs...)\n\n\nReturns an ODEFunction for spacecraft attitude dynamics.\n\nExtended Help\n\nUsage\n\nThe stm and name keyword arguments are passed to Attitude. All other keyword arguments are passed directly to SciMLBase.ODEFunction.\n\nf = AttitudeFunction()\nlet u = randn(7), p = randn(15), t = NaN # time invariant\n    f(u, p, t)\nend\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.AttitudeSystem-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.AttitudeSystem","text":"AttitudeSystem(; stm, name, defaults, kwargs...)\n\n\nA ModelingToolkit.ODESystem for atmospheric entry. Currently, only exponential atmosphere models are provided! The output model is cached with Memoize.jl. Planet-specific parameters default to Earth values.\n\nThe order of the states follows: [q₁, q₂, q₃, q₄, ω₁, ω₂, ω₃].\n\nThe order of the parameters follows: []\n\nExtended Help\n\nThis model describes how an object moves through an exponential atmosphere, above a spherical planet.\n\nStates\n\nq: scalar-last attitude quaternion\nω: body rates (radians per second)\n\nParameters\n\nJ: inertial matrix\nL: lever arm where input torque is applied\nf: torques on the vehicle body (Newton-meters)\n\nUsage\n\nmodel = Attitude()\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.CR3BFunction-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.CR3BFunction","text":"CR3BFunction(; stm, name, kwargs...)\n\n\nReturns an ODEFunction for CR3B dynamics.\n\nThe order of the states follows: [μ].\n\nThe order of the parameters follows: [μ].\n\nExtended Help\n\nUsage\n\nThe stm, and name keyword arguments are passed to CR3B. All other keyword arguments are passed directly to SciMLBase.ODEFunction.\n\nf = CR3BFunction(; stm=false, jac=true)\nlet u = randn(6), p = randn(1), t = 0\n    f(u, p, t)\nend\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.CR3BProblem-Tuple{Any, Any, Any}","page":"Docstrings","title":"AstrodynamicalModels.CR3BProblem","text":"CR3BProblem(u0, tspan, p; kwargs...)\n\n\nReturn an ODEProblem for the provided CR3B system.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.CR3BSystem-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.CR3BSystem","text":"CR3BSystem(; stm, name, defaults, kwargs...)\n\n\nA ModelingToolkit.ODESystem for the Circular Restricted Three-body Problem.\n\nThe order of the states follows: [x, y, z, ẋ, ẏ, ż].\n\nThe order of the parameters follows: [μ].\n\nExtended Help\n\nThe Circular Restricted Three-body Problem is a simplified dynamical model describing one small body (spacecraft, etc.) and two celestial bodies moving in a circle about their common center of mass. This may seem like an arbitrary simplification, but this assumption holds reasonably well for the Earth-Moon, Sun-Earth, and many other systems in our solar system.\n\nUsage\n\nmodel = CR3BSystem(; stm=true)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.NBFunction-Tuple{Int64}","page":"Docstrings","title":"AstrodynamicalModels.NBFunction","text":"NBFunction(N; stm, name, kwargs...)\n\n\nReturns an ODEFunction for NBP dynamics. The order of states and parameters in the ODEFunction arguments are equivalent to the order of states and parameters for the system produced with NBP(N). As a general rule, the order of the states follows: [x₁, y₁, z₁, ..., xₙ, yₙ, zₙ, ẋ₁, ẏ₁, ż₁, ..., ẋₙ, ẏₙ, żₙ].\n\nnote: Note\nUnlike R2BP and CR3BP, jac is set to false by default. The number of states for NBP systems can be very large for relatively small numbers of bodies (N). Enabling jac=true by default would cause unnecessarily long waiting times for this @memoize function to return for N ≥ 3 or so. If N=2 and stm=true, setting jac=true could still result in several minutes of calculations, depending on the computer you're using.\n\nwarning: Warning\nBe careful about specifying stm=true for systems with N ≥ 3! If state transition matrix dynamics are enabled, you can calculate the total number of system states with N*6 + (N*6)^2. Note that this increases exponentially as N grows! For N == 6, unless you're using parallelization, your computer may run for several hours.\n\nExtended Help\n\nUsage\n\nThe stm, and name keyword arguments are passed to NBP. All other keyword arguments are passed directly to SciMLBase.ODEFunction.\n\nf = NBFunction(3; stm=false, name=:NBP, jac=false, sparse=false)\nlet u = randn(3*6), p = randn(1 + 3), t = 0\n    f(u, p, t)\nend\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.NBSystem-Tuple{Int64}","page":"Docstrings","title":"AstrodynamicalModels.NBSystem","text":"NBSystem(N; stm, name, defaults, kwargs...)\n\n\nA ModelingToolkit.ODESystem for the Newtonian N-body Problem.\n\nThe order of the states follows: [x₁, y₁, z₁, ..., xₙ, yₙ, zₙ, ẋ₁, ẏ₁, ż₁, ..., ẋₙ, ẏₙ, żₙ].\n\nThe order of the parameters follows: [G, m₁, m₂, ..., mₙ].\n\nwarning: Warning\nBe careful about specifying stm=true for systems with N ≥ 3! If state transition matrix dynamics are enabled, you can calculate the total number of system states with N*6 + (N*6)^2. Note that this increases exponentially as N grows! For N == 6, unless you're using parallelization, your computer may run for several hours.\n\nExtended Help\n\nThe N-body problem is a model which describes how N bodies will move with respect to a common origin. This problem typically involves many bodies which act due to one force: electromagentism, gravity, etc. This model applies most closely to many celestial bodies moving due to gravity. That's about right for a model in a package called AstrodynamicalModels!\n\nUsage\n\n# One model for ALL the planets in our solar system 😎\nmodel = NBSystem(9)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.PlanarEntryFunction-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.PlanarEntryFunction","text":"PlanarEntryFunction(; name, kwargs...)\n\n\nReturns an ODEFunction for Planar Entry dynamics. Results are cached with Memoize.jl.\n\nThe order of the states follows: [γ, v, r, θ].\n\nThe order of the parameters follows: [R, P, H, m, A, C, μ]\n\nExtended Help\n\nUsage\n\nThe name keyword argument is ]passed to PlanarEntry. All other keyword arguments are passed directly to SciMLBase.ODEFunction.\n\nf = PlanarEntryFunction()\nlet u = randn(4), p = randn(7), t = NaN # time invariant\n    f(u, p, t)\nend\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.PlanarEntrySystem-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.PlanarEntrySystem","text":"PlanarEntrySystem(; name, defaults, kwargs...)\n\n\nA ModelingToolkit.ODESystem for atmospheric entry. Currently, only exponential atmosphere models are provided! The output model is cached with Memoize.jl. Planet-specific parameters default to Earth values.\n\nThe order of the states follows: [γ, v, r, θ].\n\nThe order of the parameters follows: [R, P, H, m, A, C, μ]\n\nExtended Help\n\nThis model describes how an object moves through an exponential atmosphere, above a spherical planet.\n\nUsage\n\nmodel = PlanarEntrySystem()\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.R2BFunction-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.R2BFunction","text":"R2BFunction(; stm, name, kwargs...)\n\n\nReturns an ODEFunction for R2B dynamics.\n\nThe order of the states follows: [x, y, z, ẋ, ẏ, ż].\n\nThe order of the parameters follows: [μ].\n\nExtended Help\n\nUsage\n\nThe stm, and name keyword arguments are passed to R2B. All other keyword arguments are passed directly to SciMLBase.ODEFunction.\n\nf = R2BFunction(; stm=false, name=:R2B, jac=true)\nlet u = randn(6), p = randn(1), t = 0\n    f(u, p, t)\nend\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.R2BProblem-Tuple{Any, Any, Any}","page":"Docstrings","title":"AstrodynamicalModels.R2BProblem","text":"R2BProblem(u0, tspan, p; kwargs...)\n\n\nReturn an ODEProblem for the provided R2B system.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.R2BSystem-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.R2BSystem","text":"R2BSystem(; stm, name, defaults, kwargs...)\n\n\nA ModelingToolkit.ODESystem for the Restricted Two-body Problem.\n\nThe order of the states follows: [x, y, z, ẋ, ẏ, ż].\n\nThe order of the parameters follows: [μ].\n\nExtended Help\n\nThe Restricted Two-body Problem is a simplified dynamical model describing one small body (spacecraft, etc.) and one celestial body. The gravity of the celestial body exhibits a force on the small body. This model is commonly used as a simplification to descibe our solar systems' planets orbiting our sun, or a spacecraft orbiting Earth.\n\nUsage\n\nmodel = R2BSystem()\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.dynamics-Tuple{AstrodynamicalModels.AstrodynamicalOrbit, Vararg{Any}}","page":"Docstrings","title":"AstrodynamicalModels.dynamics","text":"dynamics(orbit, args; kwargs...)\n\n\nReturn the underlying dynamics of the system in the form of a ModelingToolkit.ODEFunction.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.parameters-Tuple{Orbit}","page":"Docstrings","title":"AstrodynamicalModels.parameters","text":"parameters(orbit)\n\n\nReturn the parameter vector for an Orbit.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.state-Tuple{Orbit}","page":"Docstrings","title":"AstrodynamicalModels.state","text":"state(orbit)\n\n\nReturn the state vector for an Orbit.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.system-Tuple{AstrodynamicalModels.AstrodynamicalOrbit, Vararg{Any}}","page":"Docstrings","title":"AstrodynamicalModels.system","text":"system(orbit, args; kwargs...)\n\n\nReturn the underlying dynamics of the system in the form of a ModelingToolkit.ODESystem.\n\n\n\n\n\n","category":"method"},{"location":"#AstrodynamicalModels.jl","page":"Getting Started","title":"AstrodynamicalModels.jl","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Common models within astrodynamics!","category":"page"},{"location":"#Overview","page":"Getting Started","title":"Overview","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"This package extends ModelingToolkit to represent common astrodynamical models. All available models are shown on the Docstrings page. Consult the Models pages for more detail about each model in this package!","category":"page"},{"location":"#Usage","page":"Getting Started","title":"Usage","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"If you're familiar with ModelingToolkit.jl, then you'll be able to use this package! Some AstrodynamicalModels-specific usage instructions are provided here. Please don't be shy about making Discourse posts, or filing issues on GitHub!","category":"page"},{"location":"#Installation-and-Setup","page":"Getting Started","title":"Installation & Setup","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"This package can be installed just like any other registered Julia package.","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"# To install wherever Julia code runs...\nimport Pkg\nPkg.add(\"AstrodynamicalModels\") # or ]add AstrodynamicalModels in Julia's REPL","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"To load the package, simply enter using AstrodynamicalModels.","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using AstrodynamicalModels","category":"page"},{"location":"#Retrieving-a-Model","page":"Getting Started","title":"Retrieving a Model","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Each model within this package is implemented with a function – each function returns some AbstractSystem from ModelingToolkit.jl. Typically, this will be an ODESystem. If you're worried about overhead from calling each function every time you need a particular model, don't! Each function is implemented with @memoize, so all results are cached the first time you call a model's function with a particular function signature.","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"R2BPModel = R2BP() # Restricted Two-body Problem dynamics\n\nCR3BPModel = CR3BP() # Circular Restricted Three-body Problem dynamics\n\nCR3BPModelWithSTM = CR3BP(; stm=true) # Optionally include state transition matrix dynamics","category":"page"},{"location":"#Using-a-Model","page":"Getting Started","title":"Using a Model","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"To actually use each model, you probably also want to load ModelingToolkit (and any other SciML packages of your choice).","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using ModelingToolkit","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Now you can use any method defined for ModelingToolkit.AbstractSystem instances. Once again, the ModelingToolkit Documentation are the best place to learn how to interact with AbstractSystem instances! Some quick examples are shown below.","category":"page"},{"location":"#Check-the-Equations-of-Motion","page":"Getting Started","title":"Check the Equations of Motion","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"eqs = equations(R2BSystem())","category":"page"},{"location":"#List-the-States-and-Parameters","page":"Getting Started","title":"List the States and Parameters","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"x = states(R2BSystem())\np = parameters(R2BSystem())","category":"page"},{"location":"#Calculate-the-Jacobian","page":"Getting Started","title":"Calculate the Jacobian","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"J = calculate_jacobian(R2BSystem())","category":"page"},{"location":"#Generate-Code-to-Replicate-the-Model","page":"Getting Started","title":"Generate Code to Replicate the Model","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"print(build_function(R2BSystem()))","category":"page"},{"location":"#Generate-Code-which-Implements-the-Dynamics","page":"Getting Started","title":"Generate Code which Implements the Dynamics","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"print(R2BFunction())","category":"page"},{"location":"#Generate-C/C-and-MATLAB-Code","page":"Getting Started","title":"Generate C/C++ and MATLAB Code","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"print(build_function([eq.rhs for eq in equations(R2BSystem())], states(R2BSystem()), parameters(R2BSystem()); target=Symbolics.CTarget()))\nprint(build_function([eq.rhs for eq in equations(R2BSystem())], states(R2BSystem()), parameters(R2BSystem()); target=Symbolics.MATLABTarget()))","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"If you're interested in learning a bit about each astrodynamical model, or you'd like more specific examples which show how to use each model, consult the Models pages!","category":"page"},{"location":"NBP/#N-body-Problem-Dynamics","page":"NBP","title":"N-body Problem Dynamics","text":"","category":"section"},{"location":"NBP/","page":"NBP","title":"NBP","text":"Also known as NBP dynamics!","category":"page"},{"location":"NBP/","page":"NBP","title":"NBP","text":"using AstrodynamicalModels\nusing ModelingToolkit\nusing Latexify\nLatexify.auto_display(true)","category":"page"},{"location":"NBP/#Overview","page":"NBP","title":"Overview","text":"","category":"section"},{"location":"NBP/","page":"NBP","title":"NBP","text":"In an astrodynamical context, the N-body problem assumes N celestial bodies which move with respect to some common origin. A body i moves due to the cumulative gravity of every other body in the system. This problem is notoriously difficult because it cannot be solved analytically for Ngeq3!","category":"page"},{"location":"NBP/#Examples","page":"NBP","title":"Examples","text":"","category":"section"},{"location":"NBP/","page":"NBP","title":"NBP","text":"All NBSystem calls require the number of bodies to be specified as the first argument, like so. As always, use the stm argument at your leisure. Beware, though! using stm=true for N-body systems with more than 5 bodies may cause NBSystem to compute for a really, really long time! True story – v1.0.1 of this package had a version of these docs which tried to compute NBSystem(30; stm=true). That resulted in GitHub and JuliaHub failing the job on a timeout after several hours – at first I was surprised, until I realized that appending state transition matrix dynamics to a 30-body system results in a total of 32580 states!","category":"page"},{"location":"NBP/","page":"NBP","title":"NBP","text":"If you're curious, the number of states of an N-body system with state transition matrix dynamics appended is equivalent to N*6 + (N*6)^2.","category":"page"},{"location":"NBP/","page":"NBP","title":"NBP","text":"model = NBSystem(2; stm=true)","category":"page"},{"location":"NBP/","page":"NBP","title":"NBP","text":"Like other models, we can compute the Jacobian for these dynamics.","category":"page"},{"location":"NBP/","page":"NBP","title":"NBP","text":"using SparseArrays\nJ = sparse(calculate_jacobian(NBSystem(4)))","category":"page"},{"location":"NBP/","page":"NBP","title":"NBP","text":"Finally, let's construct a Julia function which implements these dynamics!","category":"page"},{"location":"NBP/","page":"NBP","title":"NBP","text":"f = NBFunction(2)\nlet u = randn(12), m = randn(2), G = rand(), t = 0\n    f(u, [G, m...], t)\nend","category":"page"},{"location":"R2BP/#Restricted-Two-body-Dynamics","page":"R2BP","title":"Restricted Two-body Dynamics","text":"","category":"section"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Also known as R2BP dynamics!","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"using AstrodynamicalModels\nusing ModelingToolkit\nusing Latexify\nLatexify.auto_display(true)","category":"page"},{"location":"R2BP/#Overview","page":"R2BP","title":"Overview","text":"","category":"section"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"The Restricted Two-body Problem (R2BP) assumes a massless spacecraft which moves due to the gravity of one celestial body: one star, or one planet, or one moon, or one asteroid. The equations of motion for R2BP dynamics are shown below.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"beginaligned\nfracdx(t)dt = ẋleft( t right) \nfracdy(t)dt = ẏleft( t right) \nfracdz(t)dt = żleft( t right) \nfracdẋ(t)dt = frac - mu xleft( t right)left( sqrtx^2left(tright) + y^2left(tright) + z^2left(tright) right)^3 \nfracdẏ(t)dt = frac - mu yleft( t right)left( sqrtx^2left(tright) + y^2left(tright) + z^2left(tright) right)^3 \nfracdż(t)dt = frac - mu zleft( t right)left( sqrtx^2left(tright) + y^2left(tright) + z^2left(tright) right)^3\nendaligned","category":"page"},{"location":"R2BP/#Examples","page":"R2BP","title":"Examples","text":"","category":"section"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"model = R2BSystem()","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Every model also offers optional state transition matrix dynamics. Use stm=true to append the state transition matrix dynamics to your model's equations of motion. State transition dynamics can also be thought of the model's local linearization.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"note: Note\nThe state transition dynamics for R2BSystem are not nearly as useful as the state transition dynamics within CR3BP models. Within CR3BP dynamics, a spacecraft's local linearization offers stability characteristics for periodic orbits, and provides stable and unstable directions (in state-space) for invariant manifolds about periodic orbits and Lagrange points.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"model = R2BSystem(; stm=true)","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Let's compute the Jacobian for these dynamics.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"J = calculate_jacobian(R2BSystem())","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Finally, let's construct a Julia function which implements these dynamics!","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"f = R2BFunction()\nlet u = randn(6), p = [3e6], t = 0\n    f(u, p, t)\nend","category":"page"},{"location":"CR3BP/#Circular-Restricted-Three-body-Dynamics","page":"CR3BP","title":"Circular Restricted Three-body Dynamics","text":"","category":"section"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"Also known as CR3BP dynamics!","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"using AstrodynamicalModels\nusing ModelingToolkit\nusing Latexify\nLatexify.auto_display(true)","category":"page"},{"location":"CR3BP/#Overview","page":"CR3BP","title":"Overview","text":"","category":"section"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"The Circular Restricted Three-body Problem (CR3BP) assumes a massless spacecraft which moves due to the gravity of two celestial bodies which orbit their common center of mass. This may seem like an arbitrary model, but it's actually a pretty decent approximation for how  a spacecraft moves nearby the Earth and the Sun, the Earth and the Moon, the Sun and  Jupiter, and other systems in our solar system! The equations of motion  are provided below.","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"beginaligned\nfracdx(t)dt = ẋleft( t right) \nfracdy(t)dt = ẏleft( t right) \nfracdz(t)dt = żleft( t right) \nfracdẋ(t)dt = 2 ẏleft( t right) - left( frac1sqrtleft( mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 left( 1 - mu right) left( mu + xleft( t right) right) - left( frac1sqrtleft( -1 + mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 mu left( -1 + mu + xleft( t right) right) + xleft( t right) \nfracdẏ(t)dt =  - 2 ẋleft( t right) - left( left( frac1sqrtleft( mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 left( 1 - mu right) + left( frac1sqrtleft( -1 + mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 mu right) yleft( t right) + yleft( t right) \nfracdż(t)dt = left(  - left( frac1sqrtleft( mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 left( 1 - mu right) - left( frac1sqrtleft( -1 + mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 mu right) zleft( t right)\nendaligned","category":"page"},{"location":"CR3BP/#Examples","page":"CR3BP","title":"Examples","text":"","category":"section"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"State transition dynamics are particularly valuable for CR3BP models. Recall that the state transition matrix is simply the local linearization of a spacecraft within CR3BP dynamics. Let's look at the Jacobian (another word for \"local linearization\") below, evaluated at some random state.","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"f = CR3BPFunction(; jac=true)\nlet x = randn(6), p = rand((0.0, 0.5)), t = 0\n    f.jac(x, p, t)\nend","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"The Jacobian will always have this form (zeros in the top-left, the identity matrix in the top-right, a dense matrix in the  bottom-left, and the same sparse \"-2, 2\" matrix in the bottom-right). We can include the state transition dynamics in our model with  stm=true, initialize the state transition matrix states to the  identity matrix, and propagate our spacecraft for one periodic orbit: the result is known as the Monodromy Matrix! The Monodromy Matrix provides stability characteristics for the entire periodic orbit.","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"model = CR3BSystem(; stm=true)","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"Note that periodic orbits are not easy to find within CR3BP dynamics. Various algorithms have been developed to analytically approximate,  and numerically refine, periodic CR3BP orbits. Some of those  algorithms have already been implemented in Julia! See  OrbitalTrajectories and GeneralAstrodynamics.","category":"page"},{"location":"Attitude/#Attitude-Dynamics","page":"Attitude","title":"Attitude Dynamics","text":"","category":"section"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"Quaternion kinematics, and dynamics!","category":"page"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"using AstrodynamicalModels\nusing ModelingToolkit\nusing Latexify\nLatexify.auto_display(true)","category":"page"},{"location":"Attitude/#Overview","page":"Attitude","title":"Overview","text":"","category":"section"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"The Attitude model assumes a spacecraft with some orientation described by a scalar-last quaternion, and body rates which are small enough such that they appear constant for small numerical integration tolerance values.","category":"page"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"danger: Danger\nYou should normalize the quaternion vector at each time step using a ManifoldCallback or DiscreteCallback when simulating this model! Without normalizing, the solution will drift such that the quaternion state vector is no longer a unit quaternion. The dynamics in this  model assume a unit quaternion norm!","category":"page"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"beginaligned\n    dotq = frac12 beginbmatrix\n        0  omega_3  -omega_2  omega_1 \n        -omega_3  0  omega_1  omega_2 \n        omega_2  -omega_1  0  omega_3 \n        -omega_1  -omega_2  -omega_3  0\n    endbmatrix q \n    dotomega = -J^-1 (omegatimes) J omega + J^-1 L + u \nendaligned","category":"page"},{"location":"Attitude/#Examples","page":"Attitude","title":"Examples","text":"","category":"section"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"model = AttitudeSystem()","category":"page"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"Let's compute the Jacobian for these dynamics.","category":"page"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"J = calculate_jacobian(AttitudeSystem())","category":"page"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"Finally, let's construct a Julia function which implements these dynamics!","category":"page"},{"location":"Attitude/","page":"Attitude","title":"Attitude","text":"f = AttitudeFunction()\nlet u = randn(7), p = randn(15), t = 0\n    f(u, p, t)\nend","category":"page"},{"location":"Entry/#Planar-Entry-Dynamics","page":"Entry","title":"Planar Entry Dynamics","text":"","category":"section"},{"location":"Entry/","page":"Entry","title":"Entry","text":"Also known as canonical entry dynamics!","category":"page"},{"location":"Entry/","page":"Entry","title":"Entry","text":"using AstrodynamicalModels\nusing ModelingToolkit\nusing Latexify\nLatexify.auto_display(true)","category":"page"},{"location":"Entry/#Overview","page":"Entry","title":"Overview","text":"","category":"section"},{"location":"Entry/","page":"Entry","title":"Entry","text":"The Planar Entry model assumes a spacecraft moving in an exponential atmosphere about a spherical planet. Acceleration due to gravity is ignored. The equations of motion are shown below.","category":"page"},{"location":"Entry/","page":"Entry","title":"Entry","text":"beginaligned\n  dotgamma = frac1v left( L_m - (1 - fracv^2v_c^2) g cosgamma right) \n  dotv = -D_m - g singamma \n  dotr = v singamma \n  dottheta = fracvr cosgamma \nendaligned","category":"page"},{"location":"Entry/#Examples","page":"Entry","title":"Examples","text":"","category":"section"},{"location":"Entry/","page":"Entry","title":"Entry","text":"model = PlanarEntrySystem()","category":"page"},{"location":"Entry/","page":"Entry","title":"Entry","text":"Let's compute the Jacobian for these dynamics.","category":"page"},{"location":"Entry/","page":"Entry","title":"Entry","text":"J = calculate_jacobian(PlanarEntrySystem())","category":"page"},{"location":"Entry/","page":"Entry","title":"Entry","text":"Finally, let's construct a Julia function which implements these dynamics!","category":"page"},{"location":"Entry/","page":"Entry","title":"Entry","text":"f = PlanarEntryFunction()\nlet u = abs.(randn(4)), p = abs.(randn(7)), t = 0\n    f(u, p, t)\nend","category":"page"}]
}
