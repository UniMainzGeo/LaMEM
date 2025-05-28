# Load package that contains LaMEM I/O routines
using GeophysicalModelGenerator, SpecialFunctions, LaMEM  

function CreateMarkers_Subduction(dir="./", ParamFile="test.dat"; NumberCores=1,  mpiexec="mpiexec", is64bit=false)

    cur_dir = pwd()
    cd(dir)

    # Load LaMEM particles grid
    #ParamFile_2 =   "Subduction_MATLAB_Particles.dat"
    Grid        =   read_LaMEM_inputfile(ParamFile)

    # Geometry- related parameters
    ThickCrust          =  10;        # or thickness crust
    ThickWL             =  30;
    ThermalAge_Myrs     =  30;
    ThickOP             =   65;
    ThicknessPlate      =   ThickOP + ThickWL + ThickCrust;

    # Slab parameters
    w_op                =   310;
    w_max_op            =   1750;
    w_min_op            =   -1000;

    # ==========================================================================
    # PHASES
    # ==========================================================================
    Phases      =   ones(Int64, size(Grid.X))*3;     # Rock numbers
    Temp        =   ones(Float64,size(Grid.X))*1350;     # Temperature in C    
    T_surface   =   20;
    

    # ==========================================================================
    # SETUP GEOMETRY
    # ==========================================================================

    # Inclined part of slab        
    add_box!(Phases,Temp,Grid,
            xlim=(w_min_op-w_op, w_min_op), 
            zlim=(-ThicknessPlate   , 0.0),
            Origin = (w_min_op, 0.0, 0.0),
            DipAngle = -34,
            phase=LithosphericPhases(Layers=[ThickCrust ThickOP ThickWL], Phases=[0 1 2 3]),
            T=HalfspaceCoolingTemp(Age=ThermalAge_Myrs, Tsurface=T_surface) );               

    # Create horizontal part of slab with crust & mantle lithosphere
    add_box!(Phases,Temp,Grid,
            xlim=(w_min_op, w_max_op), 
            zlim=(-ThicknessPlate   , 0.0),
            phase=LithosphericPhases(Layers=[ThickCrust ThickOP ThickWL], Phases=[0 1 2 3]),
            T=HalfspaceCoolingTemp(Age=ThermalAge_Myrs, Tsurface=T_surface) );               

    # Air
    ind         = findall(Grid.Z .> 0.0);
    Phases[ind] .= 4;
    Temp[ind]   .= 0;

    # Save julia setup 
    Model3D     =   CartData(Grid, (Phases=Phases,Temp=Temp))   # Create LaMEM model:
    write_paraview(Model3D,"LaMEM_ModelSetup", verbose=false)   # Save model to paraview   (load with opening LaMEM_ModelSetup.vts in paraview)  

    # Save LaMEM markers
    if NumberCores==1
        # 1 core
        save_LaMEM_markers_parallel(Model3D, directory="./markers", verbose=false)                      # Create LaMEM marker input on 1 core
    else
        #> 1 cores; create partitioning file first
        #PartFile = CreatePartitioningFile_local(ParamFile, NumberCores; LaMEM_dir="../../bin/", mpiexec=mpiexec)
        #PartFile = create_partitioning_file(ParamFile, NumberCores; LaMEM_dir="../../bin/", mpiexec=mpiexec)
        PartFile = run_lamem_save_grid(ParamFile, NumberCores, verbose=true)
        
        save_LaMEM_markers_parallel(Model3D, PartitioningFile=PartFile,  directory="./markers", verbose=false)     
    end

    cd(cur_dir)
end
