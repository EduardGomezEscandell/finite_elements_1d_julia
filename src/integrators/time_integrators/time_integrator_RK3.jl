#=
This time integrator solves the transient equation using the
strong-stability preserving third order Runge-Kutta method.

Assumes that the mass martix remains constant throughout the time-step
=#

include("time_integrator.jl")

struct TimeIntegratorRK3 <: TimeIntegrator
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
end

function TimeIntegratorRK3(space_integrator::SpaceIntegrator; t_end::Float64, n_steps::Integer, kwargs...)::TimeIntegratorRK3
    return TimeIntegratorRK3(space_integrator, t_end, n_steps)
end

function get_butcher_tableau_A(self::TimeIntegratorRK3)::Matrix{Float64}
    return [[1.0],
            [0.5  0.5]]
end

function get_butcher_tableau_B(self::TimeIntegratorRK3)::Matrix{Float64}
    return [1/6 1/6 2/3]
end

function get_butcher_tableau_C(self::TimeIntegratorRK3)::Matrix{Float64}
    return [0.0 1.0 0.5]
end

function solve_step(self::TimeIntegratorRK3, u_old::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    nsubsteps = 3

    u0_f = extract_free_dof_solution(self, u_old)
    Minv = invert_mass_matrix(self, u_old, t_old+0.5*Δt)

    K = zeros(size(u_old, 1), nsubsteps)
    U = zeros(size(u_old, 1) - 1, nsubsteps)

    A = get_butcher_tableau_A(self)
    C = get_butcher_tableau_C(self)
    for substep=1:nsubsteps
        θ = C[substep]
        (K[:, substep], ul) = explicit_substep(self, Minv, u0_f, t_old, θ*Δt; u=u_old, kwargs...)
        if substep != nsubsteps
            uf = u0_f + θ*Δt* reduce(+, K[:, i] * A[i, :] for i in 1:substep)
            U[:, substep] = reconstruct_solution(self.space_integrator, uf, ul)
        end
    end

    B = get_butcher_tableau_B(self)
    Δu_f = Δt * reduce(+, B[i]*k[:, i] for i=1:nsubsteps)
    return (u0_f + Δu_f, u_l)
end
