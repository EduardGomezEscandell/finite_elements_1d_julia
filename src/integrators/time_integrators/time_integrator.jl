include("../space_integrator.jl")

abstract type TimeIntegrator end

function extract_free_dof_solution(self::TimeIntegrator, u::Vector{Float64})::Vector{Float64}
    free_nodes = map(node -> node.id, filter(node -> node.dof.free, self.space_integrator.mesh.nodes))
    return u[free_nodes]
end

function explicit_substep(self::TimeIntegrator, Minv::Matrix{Float64}, uf::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    (K_ff, K_fl, F_f, _, U_l) = build_diferential_and_source(self.space_integrator; t=t_old + Δt, Δt=Δt, u=uf, kwargs...)
    return Minv * (F_f - K_ff*uf - K_fl*U_l), U_l
end

function invert_mass_matrix(self::TimeIntegrator, u::Vector{Float64}, t::Float64; kwargs...)::Matrix{Float64}
    M, _ = build_lumped_mass_matrix(self.space_integrator; t=t, u=u, kwargs...)
    dropzeros(M)
    [M.nzval[i] =1/M.nzval[i] for (i,_) in enumerate(M.nzval)]
    return M
end

function solve_step(self::TimeIntegrator, u_old::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    error("Calling abstract base class explicit_step")
end

function integrate(self::TimeIntegrator, end_of_step_hook::Function = (u::Vector{Float64}; kwargs...) -> nothing ; kwargs...)::Vector{Float64}
    # Assembly
    initialize(self.space_integrator)
    Δt = self.time_duration / self.number_of_steps

    u = kwargs[:u0]([node.x for node in self.space_integrator.mesh.nodes])
    if isa(u, Number)
        u = u.* ones(size(self.space_integrator.mesh.nodes))
    end

    for (step, time) in enumerate(LinRange(0, self.time_duration, self.number_of_steps))
        # Assembly
        (u_f, u_l) = solve_step(self, u, time, Δt; kwargs...)
        u = reconstruct_solution(self.space_integrator, u_f, u_l)

        end_of_step_hook(u; time=time, step=step)
    end

    return u
end