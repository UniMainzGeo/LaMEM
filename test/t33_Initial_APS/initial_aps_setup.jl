using Statistics
using GeophysicalModelGenerator.Printf

function compare_APS(dir, param_file, dir_markers_APS0, dir_markers_APS1)
    cur_dir = pwd();
    cd(dir)
    outfile0 = "test_markers_no_APS"
    outfile1 = "test_markers_APS"

    # Tests must be run on 1 core as markers were created with 1 core
    out0 = run_lamem_local_test(param_file, 1, "-mark_load_file $dir_markers_APS0/mdb -out_pvd 1 -out_file_name $outfile0 -out_ptr 1 -out_ptr_APS 1", opt=true)
    data0, t0 = read_LaMEM_timestep(outfile0, 2, pwd(), fields=("plast_strain [ ]",), surf=false, last=true) 
    
    out1 = run_lamem_local_test(param_file, 1, "-mark_load_file $dir_markers_APS1/mdb -out_pvd 1 -out_file_name $outfile1 -out_ptr 1 -out_ptr_APS 1", opt=true) 
    data1, t1 = read_LaMEM_timestep(outfile1, 2, pwd(), fields=("plast_strain [ ]",), surf=false, last=true) 
    
    # get the mean value of plast_strain from each model
    APS0 = data0.fields.plast_strain
    APS1 = data1.fields.plast_strain
    cd(cur_dir)
    return mean(APS0), mean(APS1)
end


# Modified from LaMEM.jl to allow writing markers with or without APS
function CreateMarkers_t33(dir="./", ParamFile="t33_setup.dat", dir_markers="./markers"; APS=nothing, NumberCores=1, is64bit=false)
    cur_dir = pwd()
    cd(dir)

    # Load LaMEM particles grid
    Grid =   read_LaMEM_inputfile(ParamFile)

    Phases =   zeros(Int64, size(Grid.X));      # Rock numbers
    Temp =   ones(Float64,size(Grid.X))*200;     # Temperature in C
    data = (Phases=Phases,Temp=Temp)
    
    if APS !== nothing
        APS = ones(Float64,size(Grid.X))*APS;     # Accumulated Plastic Strain
        data = (Phases=Phases,Temp=Temp,APS=APS)
    end
    
    # Save julia setup 
    Model3D = CartData(Grid, data)   # Create LaMEM model:

    # Save LaMEM markers
    if NumberCores==1
        # 1 core
        save_LaMEM_markers_parallel_t33(Model3D, directory=dir_markers, verbose=false) # Create LaMEM marker input on 1 core
    else
        #> 1 cores; create partitioning file first
        PartFile = CreatePartitioningFile(ParamFile,NumberCores, LaMEM_dir="../../bin/opt/");
        save_LaMEM_markers_parallel_t33(Model3D, PartitioningFile=PartFile,  directory=dir_markers, verbose=false, is64bit=is64bit)     
    end

    cd(cur_dir)

    return nothing
end

# Modified from GeophysicalModelGenerator to allow setting the header of the binary file
function save_LaMEM_markers_parallel_t33(Grid::CartData; PartitioningFile = empty, directory = "./markers", verbose = true, is64bit = false)

    x = ustrip.(Grid.x.val[:, 1, 1])
    y = ustrip.(Grid.y.val[1, :, 1])
    z = ustrip.(Grid.z.val[1, 1, :])

    if haskey(Grid.fields, :Phases)
        Phases = Grid.fields[:Phases]
    else
        error("You must provide the field :Phases in the structure")
    end

    if haskey(Grid.fields, :Temp)
        Temp = Grid.fields[:Temp]
    else
        if verbose
            println("Field :Temp is not provided; setting it to zero")
        end
        Temp = zeros(size(Phases))
    end

    if haskey(Grid.fields, :APS)
        APS = Grid.fields[:APS]
    else
        APS = nothing
    end

    if PartitioningFile == empty
        # in case we run this on 1 processor only
        Nprocx = 1
        Nprocy = 1
        Nprocz = 1
        xc, yc, zc = x, y, z
    else
        Nprocx, Nprocy, Nprocz,
            xc, yc, zc,
            nNodeX, nNodeY, nNodeZ = get_processor_partitioning(PartitioningFile, is64bit = is64bit)
        if verbose
            @show  Nprocx, Nprocy, Nprocz, xc, yc, zc, nNodeX, nNodeY, nNodeZ
        end
    end

    Nproc = Nprocx * Nprocy * Nprocz
    num, num_i, num_j, num_k = get_numscheme(Nprocx, Nprocy, Nprocz)

    xi, ix_start, ix_end = get_ind(x, xc, Nprocx)
    yi, iy_start, iy_end = get_ind(y, yc, Nprocy)
    zi, iz_start, iz_end = get_ind(z, zc, Nprocz)

    x_start = ix_start[num_i[:]]
    y_start = iy_start[num_j[:]]
    z_start = iz_start[num_k[:]]
    x_end = ix_end[num_i[:]]
    y_end = iy_end[num_j[:]]
    z_end = iz_end[num_k[:]]

    # Loop over all processors partition
    for n in 1:Nproc
        # Extract coordinates for current processor

        part_x = ustrip.(Grid.x.val[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]])
        part_y = ustrip.(Grid.y.val[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]])
        part_z = ustrip.(Grid.z.val[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]])
        part_phs = Phases[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]]
        part_T = Temp[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]]
        if(APS !== nothing)
            part_APS = APS[x_start[n]:x_end[n], y_start[n]:y_end[n], z_start[n]:z_end[n]]
        else
            part_APS = nothing
        end
        num_particles = size(part_x, 1) * size(part_x, 2) * size(part_x, 3)

        # Information vector per processor
        num_prop = APS===nothing ? 5 : 6       # number of properties we save [x/y/z/phase/T]
        lvec_info = num_particles

        lvec_prtcls = zeros(Float64, num_prop * num_particles)

        lvec_prtcls[1:num_prop:end] = part_x[:]
        lvec_prtcls[2:num_prop:end] = part_y[:]
        lvec_prtcls[3:num_prop:end] = part_z[:]
        lvec_prtcls[4:num_prop:end] = part_phs[:]
        lvec_prtcls[5:num_prop:end] = part_T[:]
        if(part_APS !== nothing)
            lvec_prtcls[6:num_prop:end] = part_APS[:]
        end


        # Write output files
        if ~isdir(directory)
            mkdir(directory)
        end         # Create dir if not existent
        fname = @sprintf "%s/mdb.%1.8d.dat"  directory (n - 1)    # Name
        if verbose
            println("Writing LaMEM marker file -> $fname")                   # print info
        end
        lvec_output = [lvec_info; lvec_prtcls]           # one vec with info about length

        if APS===nothing
            header = 1211214
        else
            header = 1211215
        end
        PetscBinaryWrite_Vec_t33(fname, lvec_output, header)            # Write PETSc vector as binary file

    end
    return
end

function PetscBinaryWrite_Vec_t33(filename, A, header)

    # Note: use "hton" to transfer to Big Endian type, which is what PETScBinaryRead expects
    return open(filename, "w+") do f
        n = length(A)
        nummark = A[1]            # number of markers

        write(f, hton(Float64(header)))     # header (not actually used)
        write(f, hton(Float64(nummark)))     # info about # of markers written

        for i in 2:n
            write(f, hton(Float64(A[i])))   # Write data itself
        end

    end
end

# Copied from LaMEM.jl
function get_numscheme(Nprocx, Nprocy, Nprocz)
    n = zeros(Int64, Nprocx * Nprocy * Nprocz)
    nix = zeros(Int64, Nprocx * Nprocy * Nprocz)
    njy = zeros(Int64, Nprocx * Nprocy * Nprocz)
    nkz = zeros(Int64, Nprocx * Nprocy * Nprocz)

    num = 0
    for k in 1:Nprocz
        for j in 1:Nprocy
            for i in 1:Nprocx
                num = num + 1
                n[num] = num
                nix[num] = i
                njy[num] = j
                nkz[num] = k
            end
        end
    end

    return n, nix, njy, nkz
end

# Copied from LaMEM.jl
function get_ind(x, xc, Nprocx)
    if Nprocx == 1
        xi = length(x)
        ix_start = [1]
        ix_end = [length(x)]
    else

        xi = zeros(Int64, Nprocx)
        for k in 1:Nprocx
            if k == 1
                xi[k] = length(x[(x .>= xc[k]) .& (x .<= xc[k + 1])])
            else
                xi[k] = length(x[(x .> xc[k]) .& (x .<= xc[k + 1])])
            end
        end
        ix_start = cumsum([0; xi[1:(end - 1)]]) .+ 1
        ix_end = cumsum(xi[1:end])
    end


    return xi, ix_start, ix_end
end