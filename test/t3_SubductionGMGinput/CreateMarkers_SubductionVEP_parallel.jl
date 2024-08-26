# Load package that contains LaMEM I/O routines
using GeophysicalModelGenerator, SpecialFunctions  

function CreateMarkers_SubductionVEP(dir="./", ParamFile="Subduction_VEP.dat"; NumberCores=8, mpiexec="mpiexec", is64bit=false)

    cur_dir = pwd()
    cd(dir)

    # Load LaMEM particles grid
    #ParamFile_2 =   "Subduction_MATLAB_Particles.dat"
    Grid        =   read_LaMEM_inputfile(ParamFile)

    # Geometry- related parameters
    ThickAir=   20;
    
    ThickOC     =   8;  # thickness of oceanic crust [km]
    ThickSP     =   80; # thickness of subducting plate [km]
    ThickOP     =   80; # thickness of overriding plate [km]
    ThickSML    =   ThickSP - ThickOC; # thickness of subducting mantle lithosphere
    ThickOML    =   ThickOP - ThickOC; # thickness of overriding mantle lithosphere
    
    z_air       =   0.0;
    z_oc        =   z_air - ThickOC;
    z_sp        =   z_air - ThickSP;
    z_op        =   z_air - ThickOP;
    
    xtrench     =   0; # initial trench position
    theta       =   30 * pi/180; # initial subducting angle [rad]
    
    
    # Weak zone geometry
    Hweak = ThickSP + 20;
    x1 = xtrench; z1 = z_oc;
    x2 = x1-50;   z2 = z_oc;
    x3 = x2-Hweak/tan(theta); z3 = z_oc - Hweak;

    Phase       =   ones(Int64, size(Grid.X));     # Rock numbers
    Temp        =   ones(Float64,size(Grid.X))*1350;     # Temperature in C    
    
    # ==========================================================================
    # TEMPERATURE - in Celcius
    # ==========================================================================
    Z = Grid.Z;
    X = Grid.X;
    H = Grid.H;

    # Temperature - in Celsius
    Ttop        =   0
    Tmantle     =   1280  # mantle potential temperature
    dTdz        =   0.3   # adiabatic gradient [oC/km]
    Tbottom     =   Tmantle + (H-ThickAir)*dTdz  # bottom temperature

    SecYear     =   3600*24*365.25  # 1 year in seconds
    kappa       =   1e-6  # thermal diffusivity
    Tage_SP     =   80e6*SecYear  # thermal age of subducting plate
    Tage_OP     =   40e6*SecYear  # thermal age of overriding plate

    # Thermal structure in the mantle
    Temp        =   Tmantle .+ abs.(Z .- z_air)*dTdz

    # Thermal structure in subducting plate
    TSP         =   Tmantle .+ abs.(z_sp - z_air)*dTdz
    dT          =   (TSP - Ttop)*(1 - erf.(abs(z_sp-z_air)*1000/2/sqrt(Tage_SP*kappa)) )
    ind         =   findall((Z .<= z_air) .& (Z .>= z_sp))
    Temp[ind]   =   Ttop .+ (TSP+dT-Ttop).*erf.(abs.(Z[ind] .- z_air).*1000/2/sqrt(Tage_SP*kappa))

    # Thermal structure in overriding plate
    TOP         =   Tmantle .+ abs(z_op-z_air)*dTdz
    dT          =   (TOP - Ttop)*(1 - erf(abs(z_op-z_air)*1000/2/sqrt(Tage_OP*kappa)) )
    ind         =   findall((Z .>= z_op) .& (Z .<= z_air) .& (X .<= (x2 .+ (Z .- z2)./tan(theta)) ) )
    Temp[ind]   =   Ttop .+ (TOP+dT-Ttop)*erf.(abs.(Z[ind] .- z_air).*1000/2/sqrt(Tage_OP*kappa))

    # Constrain air to Ttop
    ind         =   findall(Z .>= z_air)
    Temp[ind]  .=   Ttop


    # ==========================================================================
    # PHASES
    # ==========================================================================
    Air             =   0;
    Mantle          =   1;
    WeakZone        =   2;
    OceanicCrust    =   3;
    SubductPlate    =   4;
    OverridingPlate =   5;

    # Air
    ind         =   findall(Z .> z_air);
    Phase[ind]  .=   Air;

    # Oceanic crust
    ind         =   findall( (Z .>= z_oc) .& (Z .<= z_air));
    Phase[ind]  .=   OceanicCrust;

    # Subducting plate
    ind         =   findall( (Z .>= z_sp) .& (Z .<= z_oc) .& (X .>= (x1 .+ (Z .- z1) .* (x3 .- x1) ./ (z3 .- z1))) .& (Temp .<= 1200));
    Phase[ind]  .=   SubductPlate;

    # Overriding plate
    ind         =   findall((Z .>= z_op) .& (Z .<= z_oc) .& (X .<= x2 .+ (Z .- z2) ./ tan(theta)) .& (Temp .<= 1200) );
    Phase[ind]  .=   OverridingPlate;

    # Weak zone
    ind         =   findall( (Z .<= z_oc) .& (X .>= x2 .+ (Z .- z2) ./ tan(theta)) .& (X .<= x1 .+ (Z .- z1) .* (x3 .- x1) ./ (z3 .- z1)));
    Phase[ind]  .=   WeakZone;

    

    # Save julia setup 
    Model3D     =   CartData(Grid, (Phases=Phase,Temp=Temp))   # Create LaMEM model:
    write_paraview(Model3D,"LaMEM_ModelSetup_VEP", verbose=false)              # Save model to paraview   (load with opening LaMEM_ModelSetup.vts in paraview)  

    # Save LaMEM markers
    if NumberCores==1
        # 1 core
        save_LaMEM_markers_parallel(Model3D, directory="./markers", verbose=false)                      # Create LaMEM marker input on 1 core
    else
        #> 1 cores; create partitioning file first
        #PartFile = CreatePartitioningFile(ParamFile,NumberCores, LaMEM_dir="../../bin/opt/", verbose=false);
        PartFile = CreatePartitioningFile_local(ParamFile, NumberCores; LaMEM_dir="../../bin/", mpiexec=mpiexec)
        save_LaMEM_markers_parallel(Model3D, PartitioningFile=PartFile,  directory="./markers", verbose=false, is64bit=is64bit)     
    end

    cd(cur_dir)
end
