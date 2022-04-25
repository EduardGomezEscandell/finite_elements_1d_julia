#=
This time integrator solves the transient equation using the
strong-stability preserving third order Runge-Kutta method.

Assumes that the mass martix remains constant throughout the time-step
=#

include("time_integrator.jl")

mutable struct TimeIntegratorRK3 <: TimeIntegrator
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
end

function TimeIntegratorRK3(space_integrator::SpaceIntegrator; t_end::Float64, n_steps::Integer, kwargs...)::TimeIntegratorRK3
    return TimeIntegratorRK3(space_integrator, t_end, n_steps)
end

function get_butcher_tableau_A(_::TimeIntegratorRK3)::Array{Array{Float64}}
    return [[],
            [1.0],
            [0.5 0.5]]
end

function get_butcher_tableau_B(_::TimeIntegratorRK3)::Array{Float64}
    return [1/6 1/6 2/3]
end

function get_butcher_tableau_C(_::TimeIntegratorRK3)::Array{Float64}
    return [0.0 1.0 0.5]
end

function solve_step(self::TimeIntegratorRK3, u_old::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    nsubsteps = 3

    u0_f = extract_free_dof_solution(self, u_old)
    Minv = invert_mass_matrix(self, u_old, t_old+0.5*Δt)

    K = zeros(size(u_old, 1), nsubsteps)

    A = get_butcher_tableau_A(self)
    C = get_butcher_tableau_C(self)

    u_l = zeros(Float64, 0)

    for (substep, (A_s, θ)) = enumerate(zip(A, C))
        uf = u0_f
        if substep != 1
            uf += θ*Δt* reduce(+, K[:, i] * a for (i, a) = enumerate(A_s))
        end
        (K[:, substep], ul) = explicit_substep(self, Minv, u0_f, t_old, θ*Δt; u=uf, kwargs...)
        if substep == nsubsteps
            u_l = ul
            break
        end
    end

    B = get_butcher_tableau_B(self)
    Δu_f = Δt * reduce(+, K[:, i]*b for (i, b) = enumerate(B))
    return (u0_f + Δu_f, u_l)
end
