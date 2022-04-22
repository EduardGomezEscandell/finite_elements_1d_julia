#=
This time integrator solves the steady-state equation (i.e. no inertial terms) for various values of time.
=#

include("space_integrator.jl")

struct TimeIntegratorQuasiStatic
    space_integrator::SpaceIntegrator
    time_duration::Float64
    number_of_steps::Integer
end

function integrate(self::TimeIntegratorQuasiStatic, end_of_step_hook::Function = (u::Vector{Float64}; kwargs...) -> nothing ; kwargs...)::Vector{Float64}
    # Assembly
    initialize(self.space_integrator)
    Δt = self.time_duration / self.number_of_steps
    u = zeros(Float64, 0)

    for (step, time) in enumerate(LinRange(0, self.time_duration, self.number_of_steps))
        (K_ff, K_fl, U_l, U_f, F_b) = build_diferential_and_source(self.space_integrator; t=time, Δt=Δt, kwargs...)

        system = SystemOfEquations(K_ff, K_fl, U_l, U_f, F_b)

        # Solution
        solve(system)
        u = reconstruct_solution(self.space_integrator, system)

        end_of_step_hook(u; time=time, step=step)
    end

    return u
end