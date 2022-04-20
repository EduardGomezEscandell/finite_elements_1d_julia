# Finite element implementation in Julia
This is a toy project as a first contact with the Julia programming language.

### Dependencies
You need to install the `Gaston` package and `GNUplot`:
```bash
julia -e 'using Pkg; Pkg.add("Gaston")'
sudo apt install gnuplot
```
You can also use `bash install_dependencies.sh` and this will be don automatically for you.

## How to run
Open the demo file and edit the settings to your liking. Then use:
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
to avoid recompiling the plotting library Gaston.
