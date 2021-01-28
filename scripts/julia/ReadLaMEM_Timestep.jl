module ReadLaMEM_Timestep
# The routines specified in this module read a LaMEM timestep and return 2D or 3D
# Scalar or Vector Data to julia, for further quantitative analysis
# 
#         

# Note that for this to work we call the VTK toolkit from python
# Therefore you need to make sure that you installed:
#  1) A python version that includes the vtk toolkit on your system
#      Test this with: 
#           $ python
#           $ python> import vtk
#     If that gives an error message, install vtk with "$ pip install vtk"
#
#  2) The PyCall package in Julia that is linked to the correct version of python
#      a) First tell Julia the correct python path with 
#           ENV["PYTHON"] = "/Users/kausb/anaconda3/bin/python"
#      b) Compile PyCall
#           Pkg.add("PyCall")
                  
# Make these routines easily available outside the module:
export ReadField_3D_pVTR, Read_VTR_File, Read_VTU_File, ReadField_2D_pVTR;
export ReadField_VTU, WriteNewPoints_VTU, AddNewField_VTU, WritePVD_File, TransferVTK2MAT;

# Import VTK package from python; requires it to be installed in your python distribution
using PyCall, Printf, MAT


function __init__()
    py"""
    import vtk as vtk;
    import vtk.numpy_interface.dataset_adapter as dsa;

    def vtkXMLPUnstructuredGridReader(FileName):
        reader  = vtk.vtkXMLPUnstructuredGridReader()
        reader.SetFileName(FileName)
        reader.Update()
        data    = reader.GetOutput()
        return data
    
    def vtkXMLPRectilinearGridReader(FileName):
        reader1  = vtk.vtkXMLPRectilinearGridReader()
        reader1.SetFileName(FileName)
        reader1.Update()
        data    = reader1.GetOutput()
        return data

    def WriteNewPoints_VTU(FileName, points, New):
        mesh     = vtk.vtkUnstructuredGrid()
        mesh.SetPoints(points)              

        # write file back to disk that includes the previous data
        writer   = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(FileName)
        writer.SetInputData(New.VTKObject)
        writer.Update()

    def  WrapDataObject(data): 
        New  =   dsa.WrapDataObject(data);          # create new field if required
        return New
    
    def AppendPointData(New, dataField, dataName):
        New.PointData.append(dataField,  dataName) 
        return New
    """
end


#vtk     =   pyimport("vtk")
#dsa     =   pyimport("vtk.numpy_interface.dataset_adapter");
   
" Extracts a 3D data field from a pVTR data structure
Usage:
    output = ReadField_3D_pVTR(data, FieldName, Type)

Input:
-   `data`:       Data structure obtained with Read_VTR_File
-   `FieldName`:  Exact name of the field as specified in the *.vtr file
-   `Type`:       either `Scalar` or `Vector`, specifying the type of the field
    
Output:
    In case of: 
- `Scalar`:
    - `data_field`, `x`,`y`,`z`:   3D field with data, as well as 1D vectors with x,y,z coordinates
    
- `Vector`:
    - `Vx`,`Vy`,`Vz`,`x`,`y`,`z`:    3D fields with vector components in x,y,z direction
                                as well as 1D vectors with x,y,z coordinates    
"
function ReadField_3D_pVTR(data, FieldName, Type)
    
    x               =   data.GetXCoordinates();
    y               =   data.GetYCoordinates();
    z               =   data.GetZCoordinates();
   
    nx              =   length(x);
    ny              =   length(y);
    nz              =   length(z);
    
    data_Field      =   data.GetPointData().GetArray(FieldName);
    if      Type=="Scalar"
        data_Field      =   reshape(data_Field     ,(nx,ny,nz));
        return data_Field,x,y,z

    elseif  Type=="Vector"    
        Vx              =   reshape(data_Field[:,1],(nx,ny,nz));
        Vy              =   reshape(data_Field[:,2],(nx,ny,nz));
        Vz              =   reshape(data_Field[:,3],(nx,ny,nz));
        return Vx, Vy, Vz, x, y, z
    end

end

"Reads a 3D LaMEM timestep (from pVTR file)
Usage:

data = Read_VTR_File(DirName, FileName)
    
Input:
- `DirName` :    Name of timestep directory (e.g., `Timestep_00000001_1.10000000e-01`)
- `FileName`:   Filename (e.g., `Subduction2D_direct.pvtr`)    

Output:
- `data`    :   data structure containing the full content of the VTR file
"
function Read_VTR_File(DirName, FileName)
    CurDir = pwd();

    # change to directory
    cd(DirName)

    # read data from parallel rectilinear grid
    data = py"vtkXMLPRectilinearGridReader"(FileName);      # this is how it should be done within modules
    cd(CurDir)

    return data
end

"Reads a 3D point data from LaMEM timestep (from pVTU file)
Usage:

data = Read_VTU_File(DirName, FileName)
    
Input:
- `DirName` :    Name of timestep directory (e.g., `Timestep_00000001_1.10000000e-01`)
- `FileName`:   Filename (e.g., `Subduction2D_direct_passive_tracers.pvtu`)    

Output:
- `data`    :   data structure containing the full content of the VTU file
"
function Read_VTU_File(DirName, FileName)
    CurDir = pwd();

    # change to directory
    cd(DirName)

    # read data from parallel rectilinear grid
    data = py"vtkXMLPUnstructuredGridReader"(FileName);      # this is how it should be done within modules
    #reader  = vtk.vtkXMLPUnstructuredGridReader()
    #reader.SetFileName(FileName)
    #reader.Update()
    #data    = reader.GetOutput()
    cd(CurDir)

    return data
end

"
Reads a 2D slice from a LaMEM timestep (from pVTR file)
Usage:
    ReadField_2D_pVTR(data, FieldName, Type)

Input:
-   `data`:       Data structure obtained with Read_VTR_File
-   `FieldName`:  Exact name of the field as specified in the *.vtr file
-   `Type`:       either `Scalar` or `Vector`, specifying the type of the field
    
Output:
    In case of: 
- `Scalar`:
    - `data_field`, `x`,`z`:   2D array with data, as well as 1D vectors with x,z coordinates
    
- `Vector`:
    - `Vx`,`Vz`,`x`,`z`:    2D arrays with vector components in x,z direction
                                as well as 1D vectors with x,z coordinates
"
function ReadField_2D_pVTR(data, FieldName, Type)
    
  if      Type=="Scalar"
        data_Field3D,x,y,z      =   ReadField_3D_pVTR(data, FieldName, Type);

        data_Field  = data_Field3D[:,1,:]';

        return data_Field,x,z

    elseif  Type=="Vector"    
        Vx3D,Vy3D,Vz3D,x,y,z    =   ReadField_3D_pVTR(data, FieldName, Type);
        
        Vx = Vx3D[:,1,:]';
        Vz = Vz3D[:,1,:]';
        
        return Vx,Vz,x,z
    end

end

"
Reads data from a VTU data structure
Usage:
    data_field = ReadField_VTU(data, FieldName)
Input:
-   `data`:       Data structure obtained with Read_VTR_File
-   `FieldName`:  Exact name of the field as specified in the *.vtr file
                `coords` gives the coordinates

Output:
    In case of: 
- `Scalar`:
    - `data_field`   1D array with data [3D in case of coordinates]
    
"
function ReadField_VTU(data, FieldName)
    if      FieldName=="coords"
        data_Field      =   data.GetPoints().GetData();   
    else  
        data_Field      =   data.GetPointData().GetArray(FieldName); 
    end
    
    return data_Field
end

"
Writes new (point) to a VTU file, while also includes new data files
Usage:

Input:
- `DirName` :   Name of timestep directory (e.g., `Timestep_00000001_1.10000000e-01`)
- `FileName`:   Filename (e.g., `Subduction2D_direct_passive_tracers.vtu`)    
- `data`:       Data structure obtained with Read_VTR_File
- `NewData`:    Exact name of the field as specified in the *.vtr file
                `coords` gives the coordinates
    
"
function WriteNewPoints_VTU(DirName, FileName, data, New)

   CurDir   =   pwd();  
   points   =   data.GetPoints();

   # change to directory
   cd(DirName)

   # Set coordinates of points
   py"WriteNewPoints_VTU"(FileName, points, New);

   #mesh     = vtk.vtkUnstructuredGrid()
   #mesh.SetPoints(points)              

   # write file back to disk that includes the previous data
   #writer   = vtk.vtkXMLUnstructuredGridWriter()
   #writer.SetFileName(FileName)
   #writer.SetInputData(New.VTKObject)
   #writer.Update()

   cd(CurDir)
end



"
Adds data to a VTU (Point) Dataset
Usage:

Input:
-   `data`:       Data structure obtained with Read_VTR_File
-   `FieldName`:  Exact name of the field as specified in the *.vtr file
                `coords` gives the coordinates

- `New`:
    - `data_field`   1D array with data [3D in case of coordinates]
    
"
function AddNewField_VTU(New, data, dataField, dataName)
    
    if sizeof(New)==0
        #New  =   dsa.WrapDataObject(data);          # create new field if required
        New = py"WrapDataObject"(data);
    end
    
    New = py"AppendPointData"(New, dataField, dataName);

   # New.PointData.append(dataField,  dataName);     # add data

    return New
end


"
Writes a PVD file
Usage:

Input:
-   `PVDFileName`:      PVD FileName
-   `OutNames`:         Array with output Output filenames
-   `Time_vec`:         Array with output times    
"
function WritePVD_File(PVDFileName, OutNames, Time_vec)
    
  

    # Generate PVD file from the directories & filenames ----------------------------------
    io = open(PVDFileName, "w")
    println(io, "<?xml version=\"1.0\"?>")
    println(io, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">")
    println(io, "<Collection>")

    for i=1:length(OutNames)
        local Time, File;

        Time = @sprintf("%1.6e",Time_vec[i])
        File = OutNames[i];
        println(io, "    <DataSet timestep=\"$Time\" file=\"$File\"/>")
    end
    println(io, "</Collection>")
    println(io, "</VTKFile>")
    close(io)

    # Print    
    print("Wrote PVD file: $PVDFileName \n")

end





"
    This processes all LaMEM Timestep directories within the current directory, reads the 
    VTR files within those directories, reads in some of the most common data fields
    into julia and saves those as matlab files
"
function TransferVTK2MAT(FileName)
    

    # define a macro that finds all directories with a certain pattern
    searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
    
    if length(FileName)==0
        error("No VTR filename given on command-line; please specify one! \n")
    else
        
        FileName    = "$FileName.pvtr";
    end


    print("Postprocessing vtrfiles $FileName \n")

    # Go over all Timestep subdirectories in the current directory
    Directories         =   searchdir(pwd(),"Timestep_")

    nMax                =   length(Directories)
    for iStep=1:nMax
        global Vx,Vy,Vz,x,y,z, T_C, Eta_tot, Phase, Rho, E2nd, T2nd_MPa, P_MPa, Eta_creep, P_lithos, APS, MeltFrac;

        DirName     =   Directories[iStep];
        print("Processing directory $DirName \n")

        # Extract the timestep info from the directory directory name 
        id      =   findlast("_",DirName);
        Time    =   parse(Float64,DirName[id[1]+1:length(DirName)]);

        # ------------------------------------------------------------------------------------- 
        # Read data from timestep
        data            =   Read_VTR_File(DirName, FileName);

        Vx,Vy,Vz,x,y,z      =   ReadField_3D_pVTR(data, "velocity [cm/yr]",     "Vector");  # always read

        # Read fields if available
        try 
            P_MPa,x,y,z     =   ReadField_3D_pVTR(data, "pressure [MPa]",       "Scalar");
        catch
            P_MPa           =   [];
        end
        try
            T_C,x,y,z       =   ReadField_3D_pVTR(data, "temperature [C]",      "Scalar");
        catch 
            T_C             =   [];
        end
        try 
            Eta_tot,x,y,z   =   ReadField_3D_pVTR(data, "visc_total [Pa*s]",    "Scalar");
        catch 
            Eta_tot         =   [];
        end
        try
            Phase,x,y,z     =   ReadField_3D_pVTR(data, "phase [ ]",            "Scalar");
        catch
            Phase           =   [];
        end
        try
            Rho,x,y,z       =   ReadField_3D_pVTR(data, "density [kg/m^3]",     "Scalar");
        catch
            Rho             =   [];
        end
        try
            E2nd,x,y,z      =   ReadField_3D_pVTR(data, "j2_strain_rate [1/s]", "Scalar");
        catch
            E2nd            =   [];
        end
        try
            T2nd_MPa,x,y,z  =   ReadField_3D_pVTR(data, "j2_dev_stress [MPa]1",  "Scalar");
        catch
            T2nd_MPa        =   [];
        end
        try 
            Eta_creep,x,y,z =   ReadField_3D_pVTR(data, "visc_creep [Pa*s]",    "Scalar");
        catch 
            Eta_creep       =   [];
        end
        try 
            P_lithos,x,y,z  =   ReadField_3D_pVTR(data, "litho_press [MPa]",    "Scalar");
        catch 
            P_lithos        =   [];
        end
        try 
            APS,x,y,z       =   ReadField_3D_pVTR(data, "plast_strain [ ]",     "Scalar");
        catch 
            APS             =   [];
        end
        try 
            MeltFrac,x,y,z  =   ReadField_3D_pVTR(data, "melt_fraction [ ]",     "Scalar");
        catch 
            MeltFrac        =   [];
        end

        


        # ------------------------------------------------------------------------------------- 
        # Save data in matlab format
        CurDir   =   pwd();  
        cd(DirName)
     
        # Generate MATLAB file that saves all the info on the particles
        file = matopen("LaMEM_Data.mat", "w")
        write(file, "Time_Myrs",    Time);
        # We are making the assumption that velocity is always written to file
        write(file, "x_km",         x);
        write(file, "y_km",         y);
        write(file, "z_km",         z);
        write(file, "Vx_cmYr",      Vx);
        write(file, "Vy_cmYr",      Vy);
        write(file, "Vz_cmYr",      Vz);

        # The following variables may or may not be present
        if length(P_MPa)>0
            write(file, "P_MPa",            P_MPa);
        end
        if length(T_C)>0
            write(file, "T_C",              T_C);
        end
        if length(Eta_tot)>0
            write(file, "Eta_tot",          Eta_tot);
        end
        if length(Phase)>0
            write(file, "Phase",            Phase);
        end
        if length(Rho)>0
            write(file, "Rho_kgm3",         Rho );
        end
        if length(T2nd_MPa)>0
            write(file, "T2nd_MPa",         T2nd_MPa);
        end
        if length(Eta_creep)>0
            write(file, "Eta_creep",        Eta_creep);
        end
        if length(P_lithos)>0
            write(file, "Plithos_MPa",      P_lithos);
        end
        if length(APS)>0
            write(file, "PlasticStrain",    APS);
        end
        if length(MeltFrac)>0
            write(file, "MeltFraction",     MeltFrac);
        end
        close(file)

        cd(CurDir)
    end
    print("Saved data to matlab file: LaMEM_Data.mat in all directories above \n") 


end

end # end of module definition