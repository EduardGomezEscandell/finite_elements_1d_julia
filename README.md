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
The following figure should pop up:

![image](https://user-images.githubusercontent.com/47142856/164203009-7d13d0b8-6a17-40e2-b03d-a21111d446b0.png)


After either of these, a julia interactive session will stay open, and you can subsequently call
```julia
julia> include("demo.jl")
```
after changing any of the settings to avoid recompiling Gaston, the plotting library.
