include("space_integrator.jl")

struct TimeIntegratorSteady
    space_integrator::SpaceIntegrator
end

function integrate(self::TimeIntegratorSteady, end_of_step_hook::Function = (u::Vector{Float64}) -> nothing ; kwargs...)::Vector{Float64}
    # Assembly
    initialize(self.space_integrator)
    (K_ff, K_fl, U_l, U_f, F_b) = build_diferential_and_source(self.space_integrator; t=time, kwargs...)

    system = SystemOfEquations(K_ff, K_fl, U_l, U_f, F_b)

    # Solution
    u = solve(self.space_integrator, system)

    end_of_step_hook(u)
    return u
end