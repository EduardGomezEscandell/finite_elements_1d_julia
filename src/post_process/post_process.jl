using Gaston

mutable struct Plotter
    mesh::Mesh
    shape_fun::ShapeFunctions
    gauss_data::GaussData
end

function Plotter(mesh::Mesh, n_gauss::Integer)
    shape_funs = ShapeFunctions(mesh.shape_functions)
    gauss_data = get_gauss_quadrature(n_gauss)
    cache(shape_funs, gauss_data)
    return Plotter(mesh, shape_funs, gauss_data)
end

function plot_step(self::Plotter, solution::Vector{Float64}; axes_kwargs...)
    # Plots the solution at the nodes and gauss points
    X = []
    U = []
    for e in self.mesh.elems
        datapoints::Array{Tuple{Float64, Float64}} = [] # {x, u}

        nnodes = size(e.nodes, 1)
        x0 = (e.nodes[nnodes].x + e.nodes[1].x) / 2
        jacobian = (e.nodes[nnodes].x - e.nodes[1].x) / 2
        u_at_nodes = transpose([solution[n.id] for n in e.nodes])

        # Adding gauss points values
        for (ξ, N) in zip(self.gauss_data.points, eachcol(self.shape_fun.N))
            x = x0 + jacobian*ξ
            u = u_at_nodes * N
            push!(datapoints, (x, u))
        end

        # Adding nodal values
        for (n, u) in zip(e.nodes, u_at_nodes)
            push!(datapoints, (n.x, u))
        end

        sort!(datapoints, by=first)
        X = vcat(X, [d[1] for d in datapoints])
        U = vcat(U, [d[2] for d in datapoints])
    end

    return plot(X, U, w  = :l,
        Axes(grid = :on, xlabel = "'x'" , ylabel = "'u(x)'", key = :off; axes_kwargs...),
    )
end