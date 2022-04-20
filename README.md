## Finite element implementation in Julia
This is a toy project as a first contact with the Julia programming language.

To run it, use:
```bash
julia -i demo.jl
```
or
```bash
bash demo.jl
```
After either of these, a julia interactive session will stay open, and you can subsequently call
```julia
julia> include("demo.jl")
```
to avoid recompiling the plotting library.


You can also open the demo file and edit the settings.

### Missing features
- Neumann boundary conditions
- Non-uniform diffusivity
- Proper output (à la Matplotlib)