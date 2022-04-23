#=
This time integrator solves the transient equation using the theta method.
=#

include("space_integrator.jl")

struct TimeIntegratorThetaMethod
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
    θ::Float64
end

function get_time_integration_function(self::TimeIntegratorThetaMethod)::Function
    ε = 1e-8    # Floating point tolerance
    θ = self.θ

    if(θ < -ε) || ((θ-1) > ε)
        error("Theta must be: θ ∈ [0, 1], currently θ = $(θ)")
    end

    if θ <= ε  # Forward Euler: θ ≈ 0
        @info "TimeIntegratorThetaMethod: Selected Forward Euler"
        θ = 0
        return (self::TimeIntegratorThetaMethod, t_old::Float64, Δt::Float64, Uf_old::Vector{Float64}; kwargs...) -> begin
            (K_ff, K_fl, F_f, U_f, U_l) = build_diferential_and_source(self.space_integrator; t=t_old, Δt=Δt, kwargs...)
            rhs_f  = Δt*(F_f - K_ff*Uf_old - K_fl*U_l)

            (M_ff, M_fl) = build_mass_matrix(self.space_integrator; t=t_old, Δt=Δt, kwargs...)
            lhs_ff = M_ff
            lhs_fl = M_fl
            rhs_f  += M_ff*Uf_old + M_fl*U_l
            return SystemOfEquations(lhs_ff, lhs_fl, rhs_f, U_f, U_l)
        end
    end

    if θ-1 >= -ε  # Backward Euler: θ ≈ 1
        @info "TimeIntegratorThetaMethod: Selected Backward Euler"
        θ = 1
        return (self::TimeIntegratorThetaMethod, t_old::Float64, Δt::Float64, Uf_old::Vector{Float64}; kwargs...) -> begin
            (K_ff, K_fl, F_f, U_f, U_l) = build_diferential_and_source(self.space_integrator; t=t_old+Δt, Δt=Δt, kwargs...)
            lhs_ff = Δt*K_ff
            lhs_fl = Δt*K_fl
            rhs_f  = Δt*F_f

            (M_ff, M_fl) = build_mass_matrix(self.space_integrator; t=t_old+Δt, Δt=Δt, kwargs...)
            lhs_ff += M_ff
            lhs_fl += M_fl
            rhs_f  += M_ff*Uf_old + M_fl*U_l
            return SystemOfEquations(lhs_ff, lhs_fl, rhs_f, U_f, U_l)
        end
    end

    #Crank-Nicolson if θ = 0.5, arbitrary θ-method otherwise
    @info "TimeIntegratorThetaMethod: Selected θ-method with θ=$(θ)"
    return (self::TimeIntegratorThetaMethod, t_old::Float64, Δt::Float64, Uf_old::Vector{Float64}; kwargs...) -> begin
        (K_ff, K_fl, F_f, U_f, U_l) = build_diferential_and_source(self.space_integrator; t=t_old+Δt, Δt=Δt, kwargs...)
        χ = θ*Δt
        lhs_ff = χ*K_ff
        lhs_fl = χ*K_fl
        rhs_f  = χ*F_f

        (K_ff, K_fl, F_f, _, Ul_old) = build_diferential_and_source(self.space_integrator; t=t_old, Δt=Δt, kwargs...)
        χ = (1-θ)*Δt
        rhs_f  += χ*(F_f - K_ff*Uf_old - K_fl*Ul_old)

        (M_ff, M_fl) = build_mass_matrix(self.space_integrator; t=t_old+0.5*Δt, Δt=Δt, kwargs...)
        lhs_ff += M_ff
        lhs_fl += M_fl
        rhs_f  += M_ff*Uf_old + M_fl*Ul_old
        return SystemOfEquations(lhs_ff, lhs_fl, rhs_f, U_f, U_l)
    end
end

function extract_free_dof_solution(self::TimeIntegratorThetaMethod, u::Vector{Float64})::Vector{Float64}
    free_nodes = map(node -> node.id, filter(node -> node.dof.free, self.space_integrator.mesh.nodes))
    return u[free_nodes]
end

function integrate(self::TimeIntegratorThetaMethod, end_of_step_hook::Function = (u::Vector{Float64}; kwargs...) -> nothing ; kwargs...)::Vector{Float64}
    # Assembly
    initialize(self.space_integrator)
    Δt = self.time_duration / self.number_of_steps

    u = kwargs[:u0]([node.x for node in self.space_integrator.mesh.nodes])
    if isa(u, Number)
        u = u.* ones(size(self.space_integrator.mesh.nodes))
    end

    time_integration_function = get_time_integration_function(self)

    for (step, time) in enumerate(LinRange(0, self.time_duration, self.number_of_steps))
        # Assembly
        u_f = extract_free_dof_solution(self, u)
        system = time_integration_function(self, time, Δt, u_f; u=u, kwargs...)

        # Solution
        solve(system)
        u = reconstruct_solution(self.space_integrator, system)

        end_of_step_hook(u; time=time, step=step)
    end

    return u
end