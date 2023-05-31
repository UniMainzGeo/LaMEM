# Load package that contains LaMEM I/O routines
using GeophysicalModelGenerator, SpecialFunctions  

function t24_CreateMarkers(dir="./", ParamFile="test.dat"; NumberCores=1, is64bit=false, mpiexec="mpiexec")

    cur_dir = pwd()
    cd(dir)

    # Load LaMEM particles grid
    Grid        =   ReadLaMEM_InputFile(ParamFile)

    # Geometry- related parameters
    ThickAir=   20;
    ThickOC     =   8;  # thickness of oceanic crust [km]
    ThickCUC    =   15;  # thickness of continental upper crust [km]
    ThickCLC    =   20;  # thickness of continental lower crust [km]
    ThickSP     =   80; # thickness of subducting plate [km]
    ThickOCP    =   120; # thickness of overriding continental plate [km]
    ThickSML    =   ThickSP - ThickOC; # thickness of subducting mantle lithosphere
    ThickOML    =   ThickOCP - ThickCUC - ThickCLC; # thickness of overriding continental mantle lithosphere
    
    z_air       =   0;
    z_oc        =   z_air - ThickOC;
    z_cuc       =   z_air - ThickCUC;
    z_clc       =   z_air - ThickCUC-ThickCLC;
    z_sp        =   z_air - ThickSP;
    z_ocp       =   z_air - ThickOCP;
    
    xtrench     =   0; # initial trench position
    theta       =   30 * pi/180; # initial subducting angle [rad]
    
    # Slab parameters
    w_op                =   310;
    w_max_op            =   1750;
    w_min_op            =   -1000;

    # ==========================================================================
    # PHASES
    # ==========================================================================
    H,W = Grid.H, Grid.W
    Z = Grid.Z;
    X = Grid.X;
    Phase       =   ones(Int64, size(Grid.X))*3;     # Rock numbers
    Temp        =   ones(Float64,size(Grid.X))*1350;     # Temperature in C    
   
    Air              =   0;
    Mantle           =   1;
    WeakZone         =   2;
    OceanicCrust     =   3;
    ContinentalUC    =   4;
    ContinentalLC    =   5;
    SubductPlate     =   6;
    OverridingPlate  =   7;

    # Weak zone geometry
    Hweak = ThickSP + 20;
    x1 = xtrench; z1 = z_oc;
    x2 = x1-50;   z2 = z_oc;
    x3 = x2-Hweak/tan(theta); z3 = z_oc - Hweak;
        
    # Temperature - in Celsius
    Ttop        =   0
    Tmantle     =   1280  # mantle potential temperature
    dTdz        =   0.3   # adiabatic gradient [oC/km]
    Tbottom     =   Tmantle + (H-ThickAir)*dTdz  # bottom temperature
    
    SecYear     =   3600*24*365.25  # 1 year in seconds
    kappa       =   1e-6  # thermal diffusivity
    Tage_SP     =   80e6*SecYear  # thermal age of subducting plate
    Tage_OCP     =   40e6*SecYear  # thermal age of overriding plate
    
    # Thermal structure in the mantle
    Temp        =   Tmantle .+ abs.(Z .- z_air)*dTdz
    
    # Thermal structure in subducting plate
    TSP         =   Tmantle .+ abs.(z_sp - z_air)*dTdz
    dT          =   (TSP - Ttop)*(1 - erf.(abs(z_sp-z_air)*1000/2/sqrt(Tage_SP*kappa)) )
    ind         =   findall((Z .<= z_air) .& (Z .>= z_sp))
    Temp[ind]   =   Ttop .+ (TSP+dT-Ttop).*erf.(abs.(Z[ind] .- z_air).*1000/2/sqrt(Tage_SP*kappa))
    
    # Thermal structure in overriding plate
    TOCP        =   Tmantle .+ abs(z_ocp-z_air)*dTdz
    dT          =   (TOCP - Ttop)*(1 - erf(abs(z_ocp-z_air)*1000/2/sqrt(Tage_OCP*kappa)) )
    ind         =   findall((Z .>= z_ocp) .& (Z .<= z_air) .& (X .<= (x2 .+ (Z .- z2)./tan(theta)) ) )
    Temp[ind]   =   Ttop .+ (TOCP+dT-Ttop)*erf.(abs.(Z[ind] .- z_air).*1000/2/sqrt(Tage_OCP*kappa))
    
    # Constrain air to Ttop
    ind         =   findall(Z .>= z_air)
    Temp[ind]  .=   Ttop

    # ==========================================================================
    # SETUP GEOMETRY
    # ==========================================================================
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
    ind         =   findall((Z .>= z_ocp) .& (Z .<= z_oc) .& (X .<= x2 .+ (Z .- z2) ./ tan(theta)) .& (Temp .<= 1200) );
    Phase[ind]  .=   OverridingPlate;

    # Weak zone
    ind         =   findall( (Z .<= z_oc) .& (X .>= x2 .+ (Z .- z2) ./ tan(theta)) .& (X .<= x1 .+ (Z .- z1) .* (x3 .- x1) ./ (z3 .- z1)));
    Phase[ind]  .=   WeakZone;

    # Save julia setup 
    Model3D     =   CartData(Grid, (Phases=Phase,Temp=Temp))   # Create LaMEM model:
    Write_Paraview(Model3D,"LaMEM_ModelSetup", verbose=false)                  # Save model to paraview   (load with opening LaMEM_ModelSetup.vts in paraview)  

    # Save LaMEM markers
    if NumberCores==1
        # 1 core
        Save_LaMEMMarkersParallel(Model3D, directory="./markers_p$NumberCores", verbose=false)                      # Create LaMEM marker input on 1 core
    else
        #> 1 cores; create partitioning file first
        PartFile = CreatePartitioningFile_local(ParamFile, NumberCores; LaMEM_dir="../../bin/opt", mpiexec=mpiexec)
        Save_LaMEMMarkersParallel(Model3D, PartitioningFile=PartFile,  directory="./markers_p$NumberCores", verbose=false,is64bit=is64bit)     
    end

    cd(cur_dir)
end