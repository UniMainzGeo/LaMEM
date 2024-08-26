
function CreateMarkers_Temperature(dir="./", ParamFile="t12_Temperature_diffusion.dat", dir_markers="./markers_pT1"; NumberCores=1, is64bit=false)
    cur_dir = pwd()
    cd(dir)

    # Load LaMEM particles grid
    Grid        =   read_LaMEM_inputfile(ParamFile)

    Phases      =   zeros(Int64, size(Grid.X));      # Rock numbers
    Temp        =   zeros(Float64,size(Grid.X));     # Temperature in C    
    
    # Set initial T
    Temp_max    =   400 ;    # Maximum temperature perturbation[Celsius Degree]
    sigma       =   (Grid.H/2*1e3);			# Domain half-width [m]
    Temp        =   Temp_max.*exp.(-((Grid.Z.*1e3).^2.0)./sigma^2);

    # Save julia setup 
    Model3D     =   CartData(Grid, (Phases=Phases,Temp=Temp))   # Create LaMEM model:

    # Save LaMEM markers
    if NumberCores==1
        # 1 core
        Save_LaMEMMarkersParallel(Model3D, directory=dir_markers, verbose=false)                      # Create LaMEM marker input on 1 core
    else
        #> 1 cores; create partitioning file first
        PartFile = CreatePartitioningFile(ParamFile,NumberCores, LaMEM_dir="../../bin/opt/");
        Save_LaMEMMarkersParallel(Model3D, PartitioningFile=PartFile,  directory=dir_markers, verbose=false, is64bit=is64bit)     
    end

    cd(cur_dir)

    return nothing
end


function Plot_Analytics_vs_Numerics(z,T_anal, T, dir, filename="Analytics_vs_LaMEM.png")

    # Open figure 
    f = Figure(resolution = (1500, 800))
    ax = Axis(f[1, 1],  xlabel = "depth [km]", ylabel = "T [C]")
    lines!(ax, z, T_anal,  label = "Analytical T") 
    scatter!(ax, z, T,  label = "LaMEM T") 
    axislegend()

    ax = Axis(f[2, 1],  xlabel = "depth [km]", ylabel = "Tanal - Tnum [C]")
    lines!(ax, z, T_anal-T,  label = "Analytical T") 
    
    save(joinpath(dir,filename), f) 

    return f 
end 



function Analytical_1D(Z, t)

    rho=3000                     #[kg/m^3]
    sigma=50*1e3                 #[m]
    k=3                          #[W/m2]
    Cp=1050                      #[J/kg]
    TMax=400                     #[C]
    kappa= k/(Cp*rho)     #[m2/s]
    t=t*365.25*24*60*60*1e6      #[s]
    T_anal  = TMax./sqrt.(1.0 .+ 4*t*kappa/sigma^2) .* exp.(-(Z*1e3).^2 ./ (sigma^2 .+ 4*t*kappa))

    return T_anal
end
