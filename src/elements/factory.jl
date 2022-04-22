include("element.jl")
include("laplacian_element.jl")
include("unsteady_laplacian_element.jl")

function ElementFactory(
    name::String,
    id::Int64,
    nodes::Vector{Node})::Element

    if name == "Laplacian"
        return LaplacianElement(id, nodes)
    elseif name == "UnsteadyLaplacian"
        return UnsteadyLaplacianElement(id, nodes)
    end

    error("Unknown element type: $(name)")
end