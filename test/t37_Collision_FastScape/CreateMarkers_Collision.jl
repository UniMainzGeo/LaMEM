# Marker setup for the t37_Collision_FastScape test.
#
# Continent-continent collision above a subducting oceanic/continental slab,
# with a FastScape-driven free surface (surf_mode = 2). This is a direct
# translation of the original MATLAB setup script CreateMarkers_Collision.m
# (subduction-collision example) to Julia/GeophysicalModelGenerator: the phase
# and temperature fields are built with the same half-space cuts on the
# (dipping) slab boundaries X = Z/k + x0, in the same order, so the resulting
# geometry matches the original setup.

using GeophysicalModelGenerator, SpecialFunctions

function CreateMarkers_Collision(dir="./", ParamFile="Collision_FastScape.dat"; NumberCores=1, mpiexec="mpiexec", is64bit=false)

    cur_dir = pwd()
    cd(dir)

    # Load LaMEM particles grid
    Grid = read_LaMEM_inputfile(ParamFile)
    X, Y, Z = Grid.X, Grid.Y, Grid.Z

    # ----------------------------------------------------------------------
    # DOMAIN PARAMETERS (DIMENSIONAL, in km)
    # ----------------------------------------------------------------------
    ThickOC   = 8.0
    ThickSP   = 90.0
    ThickUCC  = 20.0
    ThickLCC  = 20.0
    ThickCML  = 80.0

    z_air  = 0.0
    z_oc   = z_air - ThickOC
    z_sp   = z_oc  - ThickSP
    z_ucc  = z_air - ThickUCC
    z_lcc  = z_ucc - ThickLCC
    z_left = z_lcc - 20.0
    z_cml  = z_lcc - ThickCML

    # ----------------------------------------------------------------------
    # PHASES
    # ----------------------------------------------------------------------
    Air    = 0
    WZ     = 1
    asth   = 2
    ucc1   = 3
    lcc1   = 4
    cml1   = 5
    basalt = 6
    sp     = 7
    ucc2   = 8
    lcc2   = 9
    cml2   = 10

    Phase = fill(asth, size(X))     # initialize phases (background: asthenosphere)
    Temp  = zeros(Float64, size(X))

    # ----------------------------------------------------------------------
    # SETUP GEOMETRY
    # ----------------------------------------------------------------------
    xop         = 1000.0
    wwz         = 15.0
    xl          = xop          # left side of the weak zone
    xr          = xop + wwz    # right side of the weak zone
    x_ocean     = 650.0        # left side of oceanic plate
    theta       = 20.0         # initial subducting dip angle
    theta_trans = 20.0
    k           = -tand(theta)
    k_trans     = -tand(theta_trans)

    wridge = 50.0   # distance between right slab and right boundary (asthenosphere-only strip)

    # Air
    ind = findall(Z .> z_air)
    Phase[ind] .= Air

    # Subducting continent
    # crust
    ind = findall((Z .> z_ucc) .& (Z .<= z_air) .& (X .<= x_ocean))
    Phase[ind] .= ucc1
    ind = findall((Z .> z_lcc) .& (Z .<= z_ucc) .& (X .<= x_ocean))
    Phase[ind] .= lcc1
    # mantle
    ind = findall((Z .> z_left) .& (Z .< z_lcc) .& (X .<= x_ocean))
    Phase[ind] .= cml1
    ind = findall((Z .> z_sp) .& (Z .< z_lcc) .& (X .<= x_ocean) .& (X .> wridge))
    Phase[ind] .= cml1
    ind = findall((Z .> z_cml) .& (Z .< z_sp) .& (X .<= (x_ocean .- (Z .- z_sp) ./ k_trans)) .& (X .> wridge))
    Phase[ind] .= cml1

    # Overriding continent
    # crust
    ind = findall((Z .> z_ucc) .& (Z .<= z_air) .& (X .>= (Z ./ k .+ xr)))
    Phase[ind] .= ucc2
    ind = findall((Z .> z_lcc) .& (Z .<= z_ucc) .& (X .>= (Z ./ k .+ xr)))
    Phase[ind] .= lcc2
    # right continental mantle lithosphere
    ind = findall((Z .> z_cml) .& (Z .< z_lcc) .& (X .>= (Z ./ k .+ xr)))
    Phase[ind] .= cml2

    # Subducting ocean
    # oceanic crust
    ind = findall((Z .> z_oc) .& (Z .<= z_air) .& (X .< (Z ./ k .+ xl)) .& (X .> x_ocean))
    Phase[ind] .= basalt
    # oceanic mantle lithosphere
    ind = findall((Z .> z_sp) .& (Z .<= z_oc) .& (X .<= (Z ./ k .+ xl)) .& (X .> x_ocean))
    Phase[ind] .= sp

    # weak zone
    ind = findall((X .>= (Z ./ k .+ xl)) .& (X .<= (Z ./ k .+ xr)) .& (Z .> z_sp) .& (Z .< z_air))
    Phase[ind] .= WZ

    # ----------------------------------------------------------------------
    # TEMPERATURE - in Celsius
    # ----------------------------------------------------------------------
    Ttop  = 0.0            # [C]
    Tpot  = 1300.0         # potential temperature [C]
    dTdz  = 0.5            # adiabatic temperature gradient [K/km]
    z_bot = minimum(Z)     # bottom of the domain (dimensional, km)
    Tbot  = Tpot + dTdz * (z_air - z_bot)

    Temp .= Tpot .+ (Tbot - Tpot) ./ (z_bot - z_air) .* (Z .- z_air)

    kappa   = 1.0e-6
    Myr     = 1.0e6 * 3600.0 * 24.0 * 365.25
    TspAge  = 70.0 * Myr        # slab thermal age
    Tsp     = Tpot + ThickSP * dTdz    # slab bottom temperature
    T_young = 10.0 * Myr
    w_young = 50.0              # width of the leftmost gradual-variation zone
    Thick_young = 10.0
    z_young = z_air - Thick_young

    # subducting plate: leftmost side with gradual age variation
    ind = findall((X .< w_young) .& (X .> w_young ./ (z_sp - z_young) .* (Z .- z_young)) .& (Z .> z_air - 20.0) .& (Z .< z_air))
    Temp[ind] .= (Tsp - Ttop) .* erf.((z_air .- Z[ind]) .* 1.0e3 ./ 2.0 ./ sqrt.(kappa .* (TspAge .- (TspAge - T_young) ./ w_young .* X[ind]))) .+ Ttop

    # subducting plate: constant plate age
    ind = findall((X .> w_young) .& (X .< (Z ./ k .+ xr)) .& (Z .> z_sp) .& (Z .< z_air))
    Temp[ind] .= (Tsp - Ttop) .* erf.((z_air .- Z[ind]) .* 1.0e3 ./ 2.0 ./ sqrt(kappa * TspAge)) .+ Ttop

    # continental plate temperature field, linear
    Tmoho = 550.0
    Tcml  = Tpot + (z_air - z_cml) * dTdz
    ind = findall(in.(Phase, Ref([ucc1, lcc1, ucc2, lcc2])))
    Temp[ind] .= Ttop .+ (Tmoho - Ttop) ./ (z_lcc - z_air) .* (Z[ind] .- z_air)
    ind = findall(in.(Phase, Ref([cml1, cml2])))
    Temp[ind] .= Tmoho .+ (Tcml - Tmoho) ./ (z_cml - z_lcc) .* (Z[ind] .- z_lcc)

    ind = findall((Phase .== sp) .& (Temp .> 1300.0))
    Phase[ind] .= asth

    # air
    ind = findall(Z .>= z_air)
    Temp[ind] .= Ttop

    # Save julia setup
    Model3D = CartData(Grid, (Phases=Phase, Temp=Temp))
    write_paraview(Model3D, "LaMEM_ModelSetup_Collision", verbose=false)

    # Save LaMEM markers
    if NumberCores == 1
        save_LaMEM_markers_parallel(Model3D, directory="./markers", verbose=false)
    else
        PartFile = CreatePartitioningFile_local(ParamFile, NumberCores; LaMEM_dir="../../bin/", mpiexec=mpiexec)
        save_LaMEM_markers_parallel(Model3D, PartitioningFile=PartFile, directory="./markers", verbose=false, is64bit=is64bit)
    end

    cd(cur_dir)
end
