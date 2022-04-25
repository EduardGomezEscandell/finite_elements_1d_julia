#=
This time integrator solves the transient equation using the
standard fourth order Runge-Kutta method.

Assumes that the mass martix remains constant throughout the time-step
=#

include("time_integrator.jl")

struct TimeIntegratorForwardEuler <: TimeIntegrator
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
end

function TimeIntegratorForwardEuler(space_integrator::SpaceIntegrator; t_end::Float64, n_steps::Integer, kwargs...)::TimeIntegratorForwardEuler
    return TimeIntegratorForwardEuler(space_integrator, t_end, n_steps)
end

function solve_step(self::TimeIntegratorForwardEuler, u_old::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    u_f = extract_free_dof_solution(self, u_old)
    Minv = invert_mass_matrix(self, u_old, t_old+0.5*Δt)
    (K_ff, K_fl, F_f, _, u_l) = build_diferential_and_source(self.space_integrator; t=t_old, Δt=Δt, u=u_f, kwargs...)
    return u_f + Δt * Minv * (F_f - K_ff*u_f - K_fl*u_l), u_l
end

function integrate(self::TimeIntegratorForwardEuler, end_of_step_hook::Function = (u::Vector{Float64}; kwargs...) -> nothing ; kwargs...)::Vector{Float64}
    # Assembly
    initialize(self.space_integrator)
    Δt = self.time_duration / self.number_of_steps

    u = kwargs[:u0]([node.x for node in self.space_integrator.mesh.nodes])
    if isa(u, Number)
        u = u.* ones(size(self.space_integrator.mesh.nodes))
    end

    for (step, time) in enumerate(LinRange(0, self.time_duration, self.number_of_steps))
        # Assembly
        (u_f, u_l) = solve_step(self, u, time, Δt; u=u, kwargs...)
        u = reconstruct_solution(self.space_integrator, u_f, u_l)

        end_of_step_hook(u; time=time, step=step)
    end

    return u
end