#using Glob, DelimitedFiles, WriteVTK, Interpolations, GeophysicalModelGenerator, LaMEM

export clean_directory #, changefolder, read_phase_diagram, project_onto_crosssection

""" 
    clean_directory(DirName)

Removes all LaMEM timesteps & `*.pvd` files from the directory `DirName`

"""
function clean_directory(DirName="./")
    
    CurDir = pwd();

    # change to directory
    cd(DirName)

    # pvd files
    for f in glob("*.pvd")
         rm(f)
    end

    # vts files
    for f in glob("*.vts")
        rm(f)
    end

    #timestep directories
    for f in glob("Timestep*")
        rm(f, recursive=true, force=true)
    end


    cd(CurDir)
end

#=
"""
    changefolder()

Starts a GUI on Windowss or Mac machines, which allows you to change our working directory
"""
function changefolder()
    if Sys.iswindows() 
        command = """
        Function Get-Folder(\$initialDirectory) {
            [System.Reflection.Assembly]::LoadWithPartialName("System.windows.forms")|Out-Null

            \$foldername = New-Object System.Windows.Forms.FolderBrowserDialog
            \$foldername.Description = "Select a folder"
            \$foldername.rootfolder = "MyComputer"

            if(\$foldername.ShowDialog() -eq "OK")
            {
                \$folder += \$foldername.SelectedPath
            }
            return \$folder
        }

        Get-Folder
        """
        cd(chomp(read(`powershell -Command $command`, String)))
        println(pwd())
    elseif Sys.isapple() 
        command = """
        try
            set af to (choose folder with prompt "Folder?")
            set result to POSIX path of af
        on error
            beep
            set result to "$(pwd())"
        end
        result
        """
        cd(chomp(read(`osascript -e $command`, String)))
        println(pwd())
    else
        error("This only works on windows and mac")
    end
end

"""
    out = read_phase_diagram(name::String)

Reads a phase diagram from a file `name` and returns a NamedTuple with temperature `T`, pressure `P`, melt density `ρ_melt`, solid density `ρ_solid`, density `ρ` and melt fraction `ϕ`
"""
function read_phase_diagram(name::String)

    f = open(name)

    # Read dimensions
    for i = 1:49; readline(f); end
    minT = parse(Float64,   readline(f))
    ΔT   = parse(Float64,   readline(f))
    nT   = parse(Int64,     readline(f))
    minP = parse(Float64,   readline(f))
    ΔP   = parse(Float64,   readline(f))
    nP   = parse(Int64,     readline(f))
    close(f)
    
    data = readdlm(name, skipstart=55);     # read numerical data

    # reshape
    ρ_melt  = reshape(data[:,1],(nT,nP));
    ϕ       = reshape(data[:,2],(nT,nP));
    ρ_solid = reshape(data[:,3],(nT,nP));
    T_K     = reshape(data[:,4],(nT,nP));
    P_bar   = reshape(data[:,5],(nT,nP));
    ρ       = ρ_melt.*ϕ .+ ρ_solid.*(1.0 .- ϕ);

    
    return (;T_K,P_bar,ρ_melt,ρ_solid,ϕ, ρ)
end


"""
    data_projected = project_onto_crosssection(data::CartData, Cross::CartData)

This function is helpful if you used the `cross_section` routine of the `GeophysicalModelGenerator` package to create a 2D cross-section through a 3D model, which "flattens" the cross-section to a 2D model.    
If you later want to visualize the LaMEM model results in the original 3D context (along with topography, for example), you need to project the 2D model back onto the 3D model. This function does that.
"""
function project_onto_crosssection(data::CartData, Cross::CartData)

    data_projected = deepcopy(data)
    
    # Create interpolation object from x-y data in cross-section
    x_vec        = [Cross.fields.FlatCrossSection[1], Cross.fields.FlatCrossSection[end]]

    interp_linear_x = linear_interpolation(x_vec, [Cross.x.val[1], Cross.x.val[end]], extrapolation_bc=Line())
    interp_linear_y = linear_interpolation(x_vec, [Cross.y.val[1], Cross.y.val[end]], extrapolation_bc=Line())

    # interpolate coordinates
    data_projected.x.val .= interp_linear_x.(data.x.val)
    data_projected.y.val .= interp_linear_y.(data.x.val)

    return data_projected
end


"""
    project_onto_crosssection(simulation_name::String, Cross::CartData)

Reads the output of a LaMEM simulation and projects it onto a 2D cross-section `Cross`
"""
function project_onto_crosssection(simulation_name::String, Cross::CartData)

    # read LaMEM simulation
    Timestep, FileNames, _   = read_LaMEM_simulation(simulation_name)

    pvd_filename = simulation_name*"_project.pvd"

    pvd = paraview_collection(pvd_filename)
    for (i,it) in enumerate(Timestep)
        data, t = read_LaMEM_timestep(simulation_name, it)
        data_p  = project_onto_crosssection(data, Cross)
        write_paraview(data_p, FileNames[i][1:end-5], pvd=pvd, time=t[1])
    end
    close(pvd)
    println("Created new file $(pvd_filename) with projected 2D data")

    return nothing
end
=#