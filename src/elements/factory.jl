include("element.jl")
include("laplacian_element.jl")
include("convection_diffusion_element.jl")
include("burgers_element.jl")

function ElementFactory(
    name::String,
    id::Int64,
    nodes::Vector{Node})::Element

    if name == "Laplacian"
        return LaplacianElement(id, nodes)
    elseif name=="ConvectionDiffusion"
        return ConvectionDiffusionElement(id, nodes)
    elseif name=="Burgers"
        return BurgersElement(id, nodes)
    end

    error("Unknown element type: $(name)")
end