# setup.jl: run once per machine to install the exact package versions in Manifest.toml
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
println("Setup complete! You can now run the simulations.")