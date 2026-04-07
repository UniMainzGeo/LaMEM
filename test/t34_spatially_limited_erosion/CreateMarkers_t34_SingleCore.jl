# Load package that contains LaMEM I/O routines
using GeophysicalModelGenerator, SpecialFunctions  

function CreateMarkers_t34_SingleCore(dir="./", ParamFile="spatially_limited_erosion.dat"; NumberCores=1, mpiexec="mpiexec", is64bit=false)

    cur_dir = pwd()
    cd(dir)

    # Load LaMEM particles grid
    #ParamFile_2 =   "Subduction_MATLAB_Particles.dat"
    Grid        =   read_LaMEM_inputfile(ParamFile)
    Phase       =   ones(Int64, size(Grid.X));           # Rock numbers
    Temp        =   ones(Float64,size(Grid.X))*1350;     # Temperature in C    
    
    X, Y, Z = Grid.X, Grid.Y, Grid.Z


  # add an air phase
    add_box!(Phase, Temp, Grid,
                    xlim    = (-50.0, 50.0), 
                    ylim    = (-1.0, 1.0), 
                    zlim    = (10.0, 20.0),
                    phase   = ConstantPhase(0),  
                    T       = nothing
                    )    
                    
    # add a crust phase
    add_box!(Phase, Temp, Grid,
                    xlim    = (-50.0, 50.0), 
                    ylim    = (-1.0, 1.0), 
                    zlim    = (-50.0, 10.0),
                    phase   = ConstantPhase(1),  
                    T       = nothing
                    ) 

    # Save julia setup 
    Model3D     =   CartData(Grid, (Phases=Phase,Temp=Temp))   # Create LaMEM model:
    write_paraview(Model3D,"LaMEM_ModelSetup_SpatiallyLimitedErosion", verbose=false)              # Save model to paraview   (load with opening LaMEM_ModelSetup.vts in paraview)  

    # Save LaMEM markers
    if NumberCores==1
        # 1 core
        save_LaMEM_markers_parallel(Model3D, directory="./markers", verbose=false)                      # Create LaMEM marker input on 1 core
    else
        #> 1 cores; create partitioning file first
        PartFile = CreatePartitioningFile_local(ParamFile, NumberCores; LaMEM_dir="../../bin/", mpiexec=mpiexec)
        save_LaMEM_markers_parallel(Model3D, PartitioningFile=PartFile, directory="./markers", verbose=false, is64bit=is64bit)
    end

    cd(cur_dir)
end
