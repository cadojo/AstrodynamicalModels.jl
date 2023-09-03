#
# Circular Restricted Three-body Problem models
#

"""
A `ModelingToolkit.ODESystem` for the Circular Restricted Three-body Problem.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help
The Circular Restricted Three-body Problem is a simplified dynamical model
describing one small body (spacecraft, etc.) and two celestial
bodies moving in a circle about their common center of mass.
This may seem like an arbitrary simplification, but this assumption
holds reasonably well for the Earth-Moon, Sun-Earth, and many other
systems in our solar system.

### Usage

```julia
model = CR3BP(; stm=true)
```
"""
function CR3BP(; stm=false, name=:CR3BP)

    @parameters μ
    @variables t x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = [x, y, z]
    v = [ẋ, ẏ, ż]

    eqs = vcat(
        δ.(r) .~ v,
        δ(ẋ) ~ x + 2ẏ - (μ * (x + μ - 1) * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3)) - ((x + μ) * (sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ)),
        δ(ẏ) ~ y - (2ẋ) - (y * (μ * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) + (sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ))),
        δ(ż) ~ z * (-μ * (sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) - ((sqrt(y^2 + z^2 + (x + μ)^2)^-3) * (1 - μ)))
    )

    if stm
        @variables (Φ(t))[1:6, 1:6] [description = "state transition matrix estimate"]
        Φ = Symbolics.scalarize(Φ)
        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r, v))

        LHS = map(δ, Φ)
        RHS = map(simplify, A * Φ)

        eqs = vcat(eqs, [LHS[i] ~ RHS[i] for i in eachindex(LHS)])
    end

    if string(name) == "CR3BP" && stm
        modelname = Symbol("CR3BPWithSTM")
    else
        modelname = name
    end

    if stm
        return ODESystem(
            eqs, t, vcat(r, v, vec(Φ)), [μ];
            name=modelname,
            defaults=Dict(vec(Φ .=> I(6)))
        )
    else
        return ODESystem(
            eqs, t, vcat(r, v), [μ]; name=modelname
        )
    end
end


"""
A `ModelingToolkit.ODESystem` for **Dimensioned** Circular Restricted Three-body Problem dynamics.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help
The Circular Restricted Three-body Problem is a simplified dynamical model
describing one small body (spacecraft, etc.) and two celestial
bodies moving in a circle about their common center of mass.
This may seem like an arbitrary simplification, but this assumption
holds reasonably well for the Earth-Moon, Sun-Earth, and many other
systems in our solar system.

### Usage

```julia
model = DCR3BP()
```
"""
function DCR3BP(; name=:DCR3BP)

    @parameters μ₁ μ₂ a
    @variables begin
        t, [description = "dimensioned time"],
        τ(t), [description = "normalized time"],
        T(t), [description = "dimensioned synodic period"],
        μ(t), [description = "normalized (reduced) mass parameter"],
        x(t), [description = "dimensioned x position"],
        xₙ(t), [description = "normalized x position"],
        y(t), [description = "dimensioned y position"],
        yₙ(t), [description = "normalized y position"],
        z(t), [description = "dimensioned z position"],
        zₙ(t), [description = "normalized z position"],
        ẋ(t), [description = "dimensioned x velocity"], 
        ẋₙ(t), [description = "normalized x velocity"], 
        ẏ(t), [description = "dimensioned y velocity"], 
        ẏₙ(t), [description = "normalized y velocity"], 
        ż(t), [description = "dimensioned z velocity"], 
        żₙ(t), [description = "normalized z velocity"], 
        ẍ(t), [description = "dimensioned x acceleration"], 
        ẍₙ(t), [description = "normalized x acceleration"], 
        ÿ(t), [description = "dimensioned y acceleration"], 
        ÿₙ(t), [description = "normalized y acceleration"], 
        z̈(t), [description = "dimensioned z acceleration"], 
        z̈ₙ(t), [description = "normalized z acceleration"], 
        rₚ(t), [description = "dimensioned distance to the primary body"], 
        rₛ(t), [description = "dimensioned distance to the secondary body"]
    end

    δ = Differential(t)

    n = CR3BP()

    eqs = [
        T ~ 2π * sqrt((a^3) / (μ₁ + μ₂)),
        μ ~ min(μ₁, μ₂) / (μ₁ + μ₂),
        rₚ ~ sqrt((xₙ + μ)^2 + yₙ^2 + zₙ^2) * a,
        rₛ ~ sqrt((xₙ - 1 + μ)^2 + yₙ^2 + zₙ^2) * a,
        τ ~ t / T,
        xₙ ~ x / a,
        yₙ ~ y / a,
        zₙ ~ z / a,
        ẋₙ ~ ẋ / a * T,
        ẏₙ ~ ẏ / a * T,
        żₙ ~ ż / a * T,
        δ(xₙ) ~ ẋₙ,
        δ(yₙ) ~ ẏₙ,
        δ(zₙ) ~ żₙ,
        δ(ẋₙ) ~ xₙ + 2ẏ - (μ * (xₙ + μ - 1) * (sqrt(yₙ^2 + zₙ^2 + (xₙ + μ - 1)^2)^-3)) - ((xₙ + μ) * (sqrt(yₙ^2 + zₙ^2 + (xₙ + μ)^2)^-3) * (1 - μ)),
        δ(ẏₙ) ~ yₙ - (2ẋ) - (yₙ * (μ * (sqrt(yₙ^2 + zₙ^2 + (xₙ + μ - 1)^2)^-3) + (sqrt(yₙ^2 + zₙ^2 + (xₙ + μ)^2)^-3) * (1 - μ))),
        δ(żₙ) ~ zₙ * (-μ * (sqrt(yₙ^2 + zₙ^2 + (xₙ + μ - 1)^2)^-3) - ((sqrt(yₙ^2 + zₙ^2 + (xₙ + μ)^2)^-3) * (1 - μ)))
    ]

    return ODESystem(
        eqs, t, [x, y, z, ẋ, ẏ, ż, xₙ, yₙ, zₙ, ẋₙ, ẏₙ, żₙ, T, μ, rₚ, rₛ, τ], [μ₁, μ₂, a]; name=name
    )
end

"""
Returns an `ODEFunction` for CR3BP dynamics.
Results are cached with `Memoize.jl`.

The order of the states follows: `[x, y, z, ẋ, ẏ, ż]`.

The order of the parameters follows: `[μ]`.

# Extended Help

### Usage

The `stm`, and `name` keyword arguments
are passed to `CR3BP`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = CR3BPFunction(; stm=false, jac=true)
let u = randn(6), p = randn(1), t = 0
    f(u, p, t)
end
```
"""
function CR3BPFunction(; stm=false, name=:CR3BP, kwargs...)
    defaults = (; jac=true)
    options = merge(defaults, kwargs)
    return ODEFunction(
        CR3BP(; stm=stm, name=name);
        options...
    )
end
