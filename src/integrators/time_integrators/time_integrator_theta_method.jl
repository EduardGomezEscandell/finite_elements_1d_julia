#=
This time integrator solves the transient equation using the theta method.
=#

include("time_integrator.jl")

struct TimeIntegratorThetaMethod <: TimeIntegrator
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
    θ::Float64
end

function TimeIntegratorThetaMethod(space_integrator::SpaceIntegrator; t_end::Float64, n_steps::Integer, θ::Float64, kwargs...)::TimeIntegratorThetaMethod
    return TimeIntegratorThetaMethod(space_integrator, t_end, n_steps, θ)
end

function get_time_integration_function(self::TimeIntegratorThetaMethod)::Function
    ε = 1e-8    # Floating point tolerance
    θ = self.θ

    if(θ < -ε) || ((θ-1) > ε)
        error("Theta must be: θ ∈ [0, 1], currently θ = $(θ)")
    end

    if θ <= ε  # Forward Euler: θ ≈ 0
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

function solve_step(self::TimeIntegratorThetaMethod, u_old::Vector{Float64}, t_old::Float64, Δt::Float64; kwargs...)::Tuple{Vector{Float64}, Vector{Float64}}
    time_integration_function = get_time_integration_function(self)

    # Assembly
    u_f = extract_free_dof_solution(self, u_old)
    system = time_integration_function(self, t_old, Δt, u_f; u=u_old, kwargs...)

    # Solution
    solve(system)
    return (system.u_f, system.u_l)
end
