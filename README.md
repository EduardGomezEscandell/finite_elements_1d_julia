# Finite element implementation in Julia
This is a toy project as a first contact with the Julia programming language.

### Dependencies
You need to install the `Gaston` package and `GNUplot`:

```bash
julia -e 'using Pkg; Pkg.add("Gaston")'
sudo apt install gnuplot
```

## How to run
Open the `demos/demo_steady.jl` file and edit the settings to your liking. Then use:
```bash
julia -i demos/demo_steady.jl
```
The following figure should pop up:

![demos/demo_steady.jl](https://user-images.githubusercontent.com/47142856/164203009-7d13d0b8-6a17-40e2-b03d-a21111d446b0.png)

If no figure pops up, you may be missing `gnuplot` (a warning will show in the console) or you might have run the
program without the `-i` flag. A julia interactive session will stay open, and you can subsequently call
```julia
julia> include("demos/demo_steady.jl")
```
after changing any of the settings to avoid recompiling Gaston, the plotting library.

The same process works for the other demos. For instance, here is the result for `demos/demo_unsteady.jl`:

![demos/demo_unsteady.jl](https://user-images.githubusercontent.com/47142856/164744229-e9387896-5b54-42c4-b86e-17ea35298d0c.gif)
