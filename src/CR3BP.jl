#
# Circular Restricted Three-body Problem models
# 

"""
A `ModelingToolkit.ODESystem` for the Circular Restricted Three-body Problem. 

# Extended Help
The Circular Restricted Three-body Problem is a simplified dynamical model 
describing one small body (spacecraft, etc.) and two celestial 
bodies moving in a circle about their common center of mass. 
This may seem like an arbitrary simplification, but this assumption
holds reasonably well for the Earth-Moon, Sun-Earth, and many other 
systems in our solar system.
"""
@memoize function CR3BP(; stm=false, structural_simplify=true, name=:CR3BP)

    @parameters t μ 
    @variables x(t) y(t) z(t) ẋ(t) ẏ(t) ż(t)
    δ = Differential(t)
    r = @SVector [x,y,z]
    v = @SVector [ẋ,ẏ,ż]
              
    eqs = vcat(
        δ.(r) .~ v,
        δ(ẋ)   ~ x + 2ẏ - (μ*(x + μ - 1)*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3)) - ((x + μ)*(sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ)),
        δ(ẏ)   ~ y - (2ẋ) - (y*(μ*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) + (sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ))),
        δ(ż)   ~ z*(-μ*(sqrt(y^2 + z^2 + (x + μ - 1)^2)^-3) - ((sqrt(y^2 + z^2 + (x + μ)^2)^-3)*(1 - μ)))
    )

    if stm 
        @variables Φ[1:6,1:6](t)
        Φ = Symbolics.scalarize(Φ)
        A = Symbolics.jacobian(map(el -> el.rhs, eqs), vcat(r,v))
    
        LHS = map(δ, Φ)
        RHS = map(simplify, A * Φ)

        eqs = vcat(eqs, [LHS[i] ~ RHS[i] for i in 1:length(LHS)])
    end

    if string(name) == "CR3BP" && stm 
        modelname = Symbol("CR3BPWithSTM")
    else
        modelname = name
    end
    sys = ODESystem(
        eqs, t, stm  ? vcat(r,v,Φ...) : vcat(r,v), [μ]; 
        name = modelname
    )

    return structural_simplify ? ModelingToolkit.structural_simplify(sys) : sys
end

"""
Returns an `ODEFunction` for CR3BP dynamics. 
Results are cached with `Memoize.jl`.

# Extended Help

### Usage

The `stm`, `structural_simplify`, and `name` keyword arguments 
are passed to `CR3BP`. All other keyword arguments are passed
directly to `SciMLBase.ODEFunction`.

```julia
f = CR3BPFunction(; stm=false, jac=true)
```
"""
@memoize function CR3BPFunction(; stm=false, structural_simplify=true, name=:CR3BP, kwargs...)
    defaults = (; jac=true)
    options  = merge(defaults, kwargs)
    return ODEFunction(
        CR3BP(; stm=stm, structural_simplify=structural_simplify, name=name);
        options...
    )
end