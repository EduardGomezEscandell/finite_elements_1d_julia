include("time_integrators/time_integrator_steady.jl")
include("time_integrators/time_integrator_quasistatic.jl")
include("time_integrators/time_integrator_theta_method.jl")
include("time_integrators/time_integrator_forward_euler.jl")
include("time_integrators/time_integrator_RK3.jl")
include("time_integrators/time_integrator_RK4.jl")

struct LazyConstructor{T}
    construct::Function
end

function LazyConstructor(T::Type; kwargs...)::LazyConstructor{T}
    fun = (spaceint; future_kwargs...) -> T(spaceint; kwargs..., future_kwargs...)
    return LazyConstructor{T}(fun)
end

function time_integrator_factory(name::String, space_integrator::SpaceIntegrator; kwargs...)::TimeIntegrator

    curry(f, x...) = (spaceint, kwargs...) -> f(spaceint, x...; kwargs...)

    integrators = Dict{String, LazyConstructor}([
        ("steady", LazyConstructor(TimeIntegratorSteady))
        ("quasistatic", LazyConstructor(TimeIntegratorQuasiStatic))
        ("theta-method", LazyConstructor(TimeIntegratorThetaMethod))
        ("crank-nicholson", LazyConstructor(TimeIntegratorThetaMethod; θ=0.5))
        ("backward-euler", LazyConstructor(TimeIntegratorThetaMethod; θ=0.5))
        ("forward-euler", LazyConstructor(TimeIntegratorForwardEuler))
        ("rk3", LazyConstructor(TimeIntegratorRK3))
        ("rk4", LazyConstructor(TimeIntegratorRK4))
    ])

    searchterm = lowercase(name)

    if searchterm in keys(integrators)
        return integrators[searchterm].construct(space_integrator; kwargs...)
    end

    error("Unknown time integrator: $(searchterm). Try any of: $(keys(integrators))")
end
