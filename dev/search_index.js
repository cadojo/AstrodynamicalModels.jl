var documenterSearchIndex = {"docs":
[{"location":"R2BP/#Restricted-Two-body-Dynamics","page":"R2BP","title":"Restricted Two-body Dynamics","text":"","category":"section"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Also known as R2BP dynamics!","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"using AstrodynamicalModels\nusing ModelingToolkit\nusing Latexify\nLatexify.auto_display(true)","category":"page"},{"location":"R2BP/#Overview","page":"R2BP","title":"Overview","text":"","category":"section"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"The Restricted Two-body Problem (R2BP) assumes a massless spacecraft which moves due to the gravity of one celestial body: one star, or one planet,  or one moon, or one asteroid. The equations of motion for R2BP dynamics are shown below.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"beginalign*\nfracdx(t)dt = ẋleft( t right) \nfracdy(t)dt = ẏleft( t right) \nfracdz(t)dt = żleft( t right) \nfracdẋ(t)dt = frac - mu xleft( t right)left( sqrtx^2left(tright) + y^2left(tright) + z^2left(tright) right)^3 \nfracdẏ(t)dt = frac - mu yleft( t right)left( sqrtx^2left(tright) + y^2left(tright) + z^2left(tright) right)^3 \nfracdż(t)dt = frac - mu zleft( t right)left( sqrtx^2left(tright) + y^2left(tright) + z^2left(tright) right)^3\nendalign*","category":"page"},{"location":"R2BP/#Examples","page":"R2BP","title":"Examples","text":"","category":"section"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"By default, ModelingToolkit.structural_simplify is called on every model. This typically makes solving the equations more efficient! If you really want to, you can disable this by specifying structural_simplify=false.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"model = R2BP(; structural_simplify=false)","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Every model also offers optional state transition matrix dynamics. Use stm=true to append the state transition matrix dynamics to your  model's equations of motion. State transition dynamics can also  be thought of the model's local linearization.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"note: Note\nThe state transition dynamics for R2BP systems are not  nearly as useful as the state transition dynamics within   CR3BP models. Within CR3BP dynamics,  a spacecraft's local linearization offers stability   characteristics for periodic orbits, and provides   stable and unstable directions (in state-space)  for invariant manifolds about periodic orbits and Lagrange   points.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"model = R2BP(; stm=true)","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Let's compute the Jacobian for these dynamics.","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"J = calculate_jacobian(R2BP())","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"Finally, let's construct a Julia function which implements these dynamics!","category":"page"},{"location":"R2BP/","page":"R2BP","title":"R2BP","text":"f = R2BPFunction()\nlet u = randn(6), p = [3e6], t = 0\n    f(u, p, t)\nend","category":"page"},{"location":"CR3BP/#Circular-Restricted-Three-body-Dynamics","page":"CR3BP","title":"Circular Restricted Three-body Dynamics","text":"","category":"section"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"Also known as CR3BP dynamics!","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"using AstrodynamicalModels\nusing ModelingToolkit\nusing Latexify\nLatexify.auto_display(true)","category":"page"},{"location":"CR3BP/#Overview","page":"CR3BP","title":"Overview","text":"","category":"section"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"The Circular Restricted Three-body Problem (CR3BP) assumes a massless spacecraft which moves due to the gravity of two celestial bodies which orbit their common center of mass. This may seem like an arbitrary model, but it's actually a pretty decent approximation for how  a spacecraft moves nearby the Earth and the Sun, the Earth and the Moon, the Sun and  Jupiter, and other systems in our solar system! The equations of motion  are provided below.","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"beginalign*\nfracdx(t)dt = ẋleft( t right) \nfracdy(t)dt = ẏleft( t right) \nfracdz(t)dt = żleft( t right) \nfracdẋ(t)dt = 2 ẏleft( t right) - left( frac1sqrtleft( mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 left( 1 - mu right) left( mu + xleft( t right) right) - left( frac1sqrtleft( -1 + mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 mu left( -1 + mu + xleft( t right) right) + xleft( t right) \nfracdẏ(t)dt =  - 2 ẋleft( t right) - left( left( frac1sqrtleft( mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 left( 1 - mu right) + left( frac1sqrtleft( -1 + mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 mu right) yleft( t right) + yleft( t right) \nfracdż(t)dt = left(  - left( frac1sqrtleft( mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 left( 1 - mu right) - left( frac1sqrtleft( -1 + mu + xleft( t right) right)^2 + left( yleft( t right) right)^2 + left( zleft( t right) right)^2 right)^3 mu right) zleft( t right)\nendalign*","category":"page"},{"location":"CR3BP/#Examples","page":"CR3BP","title":"Examples","text":"","category":"section"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"State transition dynamics are particularly valuable for CR3BP models. Recall that the state transition matrix is simply the local linearization of a spacecraft within CR3BP dynamics. Let's look at the Jacobian (another word for \"local linearization\") below, evaluated at some random state.","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"f = CR3BPFunction(; jac=true)\nlet x = randn(6), p = rand((0.0, 0.5)), t = 0\n    f(Val{:jac}, x, p, t)\nend","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"The Jacobian will always have this form (zeros in the top-left, the identity matrix in the top-right, a dense matrix in the  bottom-left, and the same sparse \"-2, 2\" matrix in the bottom-right). We can include the state transition dynamics in our model with  stm=true, initialize the state transition matrix states to the  identity matrix, and propagate our spacecraft for one periodic orbit: the result is known as the Monodromy Matrix! The Monodromy Matrix provides stability characteristics for the entire periodic orbit.","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"model = CR3BP(; stm=true)","category":"page"},{"location":"CR3BP/","page":"CR3BP","title":"CR3BP","text":"Note that periodic orbits are not easy to find within CR3BP dynamics. Various algorithms have been developed to analytically approximate,  and numerically refine, periodic CR3BP orbits. Some of those  algorithms have already been implemented in Julia! See  OrbitalTrajectories and GeneralAstrodynamics.","category":"page"},{"location":"#AstrodynamicalModels.jl","page":"Getting Started","title":"AstrodynamicalModels.jl","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Common models within astrodynamics!","category":"page"},{"location":"#Overview","page":"Getting Started","title":"Overview","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"This package extends ModelingToolkit to represent common astrodynamical models. All available models are shown on the Docstrings page. Consult the Models pages for more detail about each model in this package!","category":"page"},{"location":"#Usage","page":"Getting Started","title":"Usage","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"If you're familiar with ModelingToolkit.jl, then you'll be able to use this package! Some  AstrodynamicalModels-specific usage instructions are provided  here. Please don't be shy about making Discourse posts, or filing issues on GitHub!","category":"page"},{"location":"#Installation-and-Setup","page":"Getting Started","title":"Installation & Setup","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"This package can be installed just like any other  registered Julia package.","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"# To install wherever Julia code runs...\nimport Pkg\nPkg.add(\"AstrodynamicalModels\") # or ]add AstrodynamicalModels in Julia's REPL","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"To load the package, simply enter using AstrodynamicalModels.","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using AstrodynamicalModels","category":"page"},{"location":"#Retrieving-a-Model","page":"Getting Started","title":"Retrieving a Model","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Each model within this package is implemented with a function – each function returns some AbstractSystem from ModelingToolkit.jl. Typically, this will be an ODESystem. If you're worried about overhead from calling each function every time you need a particular model, don't! Each function is implemented with  @memoize, so all  results are cached the first time you call a model's function  with a particular function signature. ","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"R2BPModel = R2BP() # Restricted Two-body Problem dynamics\n\nCR3BPModel = CR3BP() # Circular Restricted Three-body Problem dynamics\n\nCR3BPModelWithSTM = CR3BP(; stm=true) # Optionally include state transition matrix dynamics","category":"page"},{"location":"#Using-a-Model","page":"Getting Started","title":"Using a Model","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"To actually use each model, you probably also want to load  ModelingToolkit (and any other SciML  packages of your choice).","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"using ModelingToolkit","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"Now you can use any method defined for ModelingToolkit.AbstractSystem instances. Once again, the ModelingToolkit Documentation are the best place to learn how to interact with AbstractSystem instances! Some quick examples are shown below. ","category":"page"},{"location":"#Check-the-Equations-of-Motion","page":"Getting Started","title":"Check the Equations of Motion","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"eqs = equations(R2BP())","category":"page"},{"location":"#List-the-States-and-Parameters","page":"Getting Started","title":"List the States and Parameters","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"x = states(R2BP())\np = parameters(R2BP())","category":"page"},{"location":"#Calculate-the-Jacobian","page":"Getting Started","title":"Calculate the Jacobian","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"J = calculate_jacobian(R2BP())","category":"page"},{"location":"#Generate-Code-to-Replicate-the-Model","page":"Getting Started","title":"Generate Code to Replicate the Model","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"print(build_function(R2BP()))","category":"page"},{"location":"#Generate-Code-which-Implements-the-Dynamics","page":"Getting Started","title":"Generate Code which Implements the Dynamics","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"print(ODEFunction(R2BP()))","category":"page"},{"location":"#Generate-C/C-and-MATLAB-Code","page":"Getting Started","title":"Generate C/C++ and MATLAB Code","text":"","category":"section"},{"location":"","page":"Getting Started","title":"Getting Started","text":"print(build_function([eq.rhs for eq in equations(R2BP())], states(R2BP()), parameters(R2BP()); target=Symbolics.CTarget()))\nprint(build_function([eq.rhs for eq in equations(R2BP())], states(R2BP()), parameters(R2BP()); target=Symbolics.MATLABTarget()))","category":"page"},{"location":"","page":"Getting Started","title":"Getting Started","text":"If you're interested in learning a bit about each astrodynamical model, or you'd like more specific examples which show how to use each model, consult the Models  pages! ","category":"page"},{"location":"docstrings/#Documentation","page":"Docstrings","title":"Documentation","text":"","category":"section"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"All docstrings!","category":"page"},{"location":"docstrings/","page":"Docstrings","title":"Docstrings","text":"Modules = [\n    AstrodynamicalModels\n]\nOrder = [:module, :type, :function]","category":"page"},{"location":"docstrings/#AstrodynamicalModels.AstrodynamicalModels","page":"Docstrings","title":"AstrodynamicalModels.AstrodynamicalModels","text":"Provides astrodynamical models as ModelingToolkit.ODESystems.  Check out the ModelingToolkit docs to learn how to use these  systems for orbit propagation with DifferentialEquations, or see GeneralAstrodynamics for some convenient orbit propagation  wrappers.\n\nExtended help\n\nExports\n\nCR3BP\nCR3BPFunction\nR2BP\nR2BPFunction\n\nImports\n\nBase\nCore\nDocStringExtensions\nLinearAlgebra\nModelingToolkit\nStaticArrays\nSymbolics\n\n\n\n\n\n","category":"module"},{"location":"docstrings/#AstrodynamicalModels.CR3BP-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.CR3BP","text":"CR3BP(; stm, structural_simplify, name)\n\n\nA ModelingToolkit.ODESystem for the Circular Restricted Three-body Problem. \n\nExtended Help\n\nThe Circular Restricted Three-body Problem is a simplified dynamical model  describing one small body (spacecraft, etc.) and two celestial  bodies moving in a circle about their common center of mass.  This may seem like an arbitrary simplification, but this assumption holds reasonably well for the Earth-Moon, Sun-Earth, and many other  systems in our solar system.\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.CR3BPFunction-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.CR3BPFunction","text":"CR3BPFunction(; stm, structural_simplify, name, kwargs...)\n\n\nReturns an ODEFunction for CR3BP dynamics.  Results are cached with Memoize.jl.\n\nExtended Help\n\nUsage\n\nThe stm, structural_simplify, and name keyword arguments  are passed to CR3BP. All other keyword arguments are passed directly to SciMLBase.ODEFunction.\n\nf = CR3BPFunction(; stm=false, jac=true)\n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.R2BP-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.R2BP","text":"R2BP(; stm, structural_simplify, name)\n\n\nA ModelingToolkit.ODESystem for the Restricted Two-body Problem. \n\nExtended Help\n\nThe Restricted Two-body Problem is a simplified dynamical model  describing one small body (spacecraft, etc.) and one celestial  body. The gravity of the celestial body exhibits a force on the  small body. This model is commonly used as a simplification to  descibe our solar systems' planets orbiting our sun, or a  spacecraft orbiting Earth. \n\n\n\n\n\n","category":"method"},{"location":"docstrings/#AstrodynamicalModels.R2BPFunction-Tuple{}","page":"Docstrings","title":"AstrodynamicalModels.R2BPFunction","text":"R2BPFunction(; stm, structural_simplify, name, kwargs...)\n\n\nReturns an ODEFunction for R2BP dynamics.  Results are cached with Memoize.jl.\n\nExtended Help\n\nUsage\n\nThe stm, structural_simplify, and name keyword arguments  are passed to R2BP. All other keyword arguments are passed directly to SciMLBase.ODEFunction.\n\nf = R2BPFunction(; stm=false, structural_simplify=true, name=:R2BP, jac=true)\n\n\n\n\n\n","category":"method"}]
}
