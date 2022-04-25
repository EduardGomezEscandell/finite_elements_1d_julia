#=
This time integrator solves the transient equation using the
standard fourth order Runge-Kutta method.

Assumes that the mass martix remains constant throughout the time-step
=#

include("time_integrator.jl")

struct TimeIntegratorRK4 <: TimeIntegrator
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
end

function TimeIntegratorRK4(space_integrator::SpaceIntegrator; t_end::Float64, n_steps::Integer, kwargs...)::TimeIntegratorRK4
    return TimeIntegratorRK4(space_integrator, t_end, n_steps)
end

function solve_step(self::TimeIntegratorRK4, u_old::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    u0_f = extract_free_dof_solution(self, u_old)
    Minv = invert_mass_matrix(self, u_old, t_old+0.5*Δt)

    θ = 0.0
    (k_1, _) = explicit_substep(self, Minv, u0_f, t_old, θ*Δt; u=u_old, kwargs...)
    uf = u0_f + θ*Δt*k_1

    θ = 0.5
    (k_2, _) = explicit_substep(self, Minv, uf, t_old, θ*Δt; u=u_old, kwargs...)

    θ = 0.5
    uf = u0_f + θ*Δt*k_2
    (k_3, _) = explicit_substep(self, Minv, uf, t_old, θ*Δt; u=u_old, kwargs...)

    θ = 1.0
    uf = u0_f + θ*Δt*k_3
    (k_4, u_l) = explicit_substep(self, Minv, uf, t_old, θ*Δt; u=u_old, kwargs...)

    Δu_f = Δt * (k_1 + 2*k_2 + 2*k_3 + k_4) / 6
    return (u0_f + Δu_f, u_l)
end
