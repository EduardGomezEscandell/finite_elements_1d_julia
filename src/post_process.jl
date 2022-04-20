using Gaston

function plot_at_nodes(mesh::Mesh, solution::Vector{Float64})
    # Plots the solution at the nodes
    x = [node.x for node in mesh.nodes]
    y = solution
    display(
        plot(x, y, w  = :lp, marker = "ecircle", legend = :Solution_at_nodes)
    )
end

function plot_solution(mesh::Mesh, solution::Vector{Float64}, shape_fun::ShapeFunctions, gauss_data::GaussData)
    # Plots the solution at the nodes and gauss points
    X = []
    U = []
    for e in mesh.elems
        datapoints::Array{Tuple{Float64, Float64}} = [] # {x, u}

        nnodes = size(e.nodes, 1)
        x0 = (e.nodes[nnodes].x + e.nodes[1].x) / 2
        jacobian = (e.nodes[nnodes].x - e.nodes[1].x) / 2
        u_at_nodes = transpose([solution[n.id] for n in e.nodes])

        # Adding gauss points values
        for (ξ, N) in zip(gauss_data.points, eachcol(shape_fun.N))
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
        Axes(grid = :on, title = "'Solution'", xlabel = "'x'" , ylabel = "'u(x)'", key = :off),
    )
end