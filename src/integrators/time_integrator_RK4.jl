#=
This time integrator solves the transient equation using the
standard fourth order Runge-Kutta method.

Assumes that the mass martix remains constant throughout the time-step
=#

include("space_integrator.jl")

struct TimeIntegratorRK4
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
end

function runge_kutta_4_step(self::TimeIntegratorRK4, u_old::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    u0_f = extract_free_dof_solution(self, u_old)
    Minv = invert_mass_matrix(self, u_old, t_old+0.5*Δt)

    θ = 0.0
    (k_1, _) = runge_kutta_substep(self, Minv, u0_f, t_old, θ*Δt; u=u_old, kwargs...)
    uf = u0_f + θ*Δt*k_1

    θ = 0.5
    (k_2, _) = runge_kutta_substep(self, Minv, uf, t_old, θ*Δt; u=u_old, kwargs...)

    θ = 0.5
    uf = u0_f + θ*Δt*k_2
    (k_3, _) = runge_kutta_substep(self, Minv, uf, t_old, θ*Δt; u=u_old, kwargs...)

    θ = 1.0
    uf = u0_f + θ*Δt*k_3
    (k_4, u_l) = runge_kutta_substep(self, Minv, uf, t_old, θ*Δt; u=u_old, kwargs...)

    Δu_f = Δt * (k_1 + 2*k_2 + 2*k_3 + k_4) / 6
    return (u0_f + Δu_f, u_l)
end

function extract_free_dof_solution(self::TimeIntegratorRK4, u::Vector{Float64})::Vector{Float64}
    free_nodes = map(node -> node.id, filter(node -> node.dof.free, self.space_integrator.mesh.nodes))
    return u[free_nodes]
end

function runge_kutta_substep(self::TimeIntegratorRK4, Minv::Matrix{Float64}, uf::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    (K_ff, K_fl, F_f, _, U_l) = build_diferential_and_source(self.space_integrator; t=t_old + Δt, Δt=Δt, u=uf, kwargs...)
    return Minv * (F_f - K_ff*uf - K_fl*U_l), U_l
end

function invert_mass_matrix(self::TimeIntegratorRK4, u::Vector{Float64}, t::Float64; kwargs...)::Matrix{Float64}
    M, _ = build_lumped_mass_matrix(self.space_integrator; t=t, u=u, kwargs...)
    dropzeros(M)
    [M.nzval[i] =1/M.nzval[i] for (i,_) in enumerate(M.nzval)]
    return M
end


function integrate(self::TimeIntegratorRK4, end_of_step_hook::Function = (u::Vector{Float64}; kwargs...) -> nothing ; kwargs...)::Vector{Float64}
    # Assembly
    initialize(self.space_integrator)
    Δt = self.time_duration / self.number_of_steps

    u = kwargs[:u0]([node.x for node in self.space_integrator.mesh.nodes])
    if isa(u, Number)
        u = u.* ones(size(self.space_integrator.mesh.nodes))
    end

    for (step, time) in enumerate(LinRange(0, self.time_duration, self.number_of_steps))
        # Assembly
        (u_f, u_l) = runge_kutta_4_step(self, u, time, Δt; kwargs...)
        u = reconstruct_solution(self.space_integrator, u_f, u_l)

        end_of_step_hook(u; time=time, step=step)
    end

    return u
end