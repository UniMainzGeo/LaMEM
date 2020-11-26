# The routines specified here read a LaMEM timestep and return 2D or 3D
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
                  

# Import VTK package from python; requires it to be installed in your python distribution
using PyCall
vtk     =   pyimport("vtk")
dsa     =   pyimport("vtk.numpy_interface.dataset_adapter");
   
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
    reader  = vtk.vtkXMLPRectilinearGridReader()
    reader.SetFileName(FileName)
    reader.Update()
    data    = reader.GetOutput()
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
    reader  = vtk.vtkXMLPUnstructuredGridReader()
    reader.SetFileName(FileName)
    reader.Update()
    data    = reader.GetOutput()
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
   mesh     = vtk.vtkUnstructuredGrid()
   mesh.SetPoints(points)              

   # write file back to disk that includes the previous data
   writer   = vtk.vtkXMLUnstructuredGridWriter()
   writer.SetFileName(FileName)
   writer.SetInputData(New.VTKObject)
   writer.Update()

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
        New  =   dsa.WrapDataObject(data);          # create new field if required
    end
    
    New.PointData.append(dataField,  dataName);     # add data

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