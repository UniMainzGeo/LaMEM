using GeophysicalModelGenerator

# Creates the initial topography file for t34_TopoDiffusion.
# The topography is a semicircular dome (half-sphere in 2D cross-section):
#   z_topo(x) = sqrt(R^2 - x^2)  for |x| < R
#   z_topo(x) = 0                 otherwise
# where R = 10 km and the sphere centre is at the surface level z = 0.
# The dome peak is 10 km high; diffusion will visibly smooth its flanks.

function t34_CreateSetup(dir="./", ParamFile="t34_TopoDiffusion.dat";
        NumberCores=1, mpiexec="mpiexec", is64bit=false)

    cur_dir = pwd()
    cd(dir)

    R      = 10.0                         # dome radius [km]
    x_vals = collect(-100.0:1.0:100.0)    # 201 points covering model x-extent
    y_vals = collect(-1.0:2.0:1.0)        # 2 points covering model y-extent

    X_t, Y_t, Z_t = xyz_grid(x_vals, y_vals, 0)

    Topo_z = zeros(length(x_vals), length(y_vals), 1)
    for i in eachindex(x_vals)
        z_surf = abs(x_vals[i]) < R ? sqrt(R^2 - x_vals[i]^2) : 0.0
        Topo_z[i, :, 1] .= z_surf
    end

    Topo = CartData(X_t, Y_t, Z_t, (Topography=Topo_z,))
    save_LaMEM_topography(Topo, "t34_topo.bin")

    cd(cur_dir)
end
