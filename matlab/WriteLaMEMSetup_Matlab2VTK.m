% Function WriteLaMEMSetup_Matlab2VTK using ASCII or BINARY format
% syntax:
%         WriteLaMEMSetup_Matlab2VTK(A,'BINARY')
%         WriteLaMEMSetup_Matlab2VTK(A,'ASCII')
% where A - structure built with ParallelMatlab_CreatePhases.m and contais:
%            W - width of domain in X-dir
%            L - length of domain in Y-dir
%            H - height of domain in Z-dir
%            nump_x - no. of particles in X-dir
%            nump_y - no. of particles in Y-dir
%            nump_z - no. of particles in Z-dir
%            Phase  - phase information of the particles
%            Temp  - Temperature information of the particles
%            x  - vector containing the x coordinates of particles
%            y  - vector containing the y coordinates of particles
%            z  - vector containing the y coordinates of particles
%            npart_x - no. of particles/cell in X-dir
%            npart_y - no. of particles/cell in Y-dir
%            npart_z - no. of particles/cell in Z-dir


function a = WriteLaMEMSetup_Matlab2VTK(A,format)
% Main function

if strcmp(format,'ASCII')
    %print in ASCII format
    WriteVTKModelSetup_ASCII(A);
elseif strcmp(format,'BINARY')
    WriteVTKModelSetup_BINARY(A);
else
    disp(['Required INCORRECT output format. Correct format: ASCII or BINARY '])
end

end


%% =======================================================================
function a = WriteVTKModelSetup_ASCII(A)
% Write rectilinear file in ASCII format

fname_vtk   = 'VTK_ModelSetup_paraview_ascii.vtr';
fid         = fopen(fname_vtk,'w','b');           % 'b': BigEndian

fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="RectilinearGrid" version="0.1" byte_order="BigEndian" >\n');

fprintf(fid,'\t<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n',1,A.nump_x,1,A.nump_y,1,A.nump_z);
fprintf(fid,'\t\t<Piece Extent=\"%d %d %d %d %d %d\">\n',1,A.nump_x,1,A.nump_y,1,A.nump_z);

fprintf(fid,'\t\t\t<CellData></CellData>\n');

disp(['Writing coordinates for Paraview... '])
fprintf(fid,'\t\t\t<Coordinates>\n');

fprintf(fid,'\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_X\" NumberOfComponents=\"1\" format=\"ascii\"/>\n');
for i=1:length(A.x)
    fprintf(fid,'\t\t\t\t\t %1.8g\n', A.x(i));
end

fprintf(fid,'\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_Y\" NumberOfComponents=\"1\" format=\"ascii\"/> \n');
for i=1:length(A.y)
    fprintf(fid,'\t\t\t\t\t %1.8g\n', A.y(i));
end

fprintf(fid,'\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_Z\" NumberOfComponents=\"1\" format=\"ascii\"/>\n');
for i=1:length(A.z)
    fprintf(fid,'\t\t\t\t\t %1.8g\n', A.z(i));
end

fprintf(fid,'\t\t\t</Coordinates>\n');
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'\t\t\t<PointData Scalars=\"Dimensional\">\n');
fprintf(fid,'\n');

disp(['Writing Phase information for Paraview... '])
fprintf(fid,'\t\t\t\t<DataArray type=\"Int32\" Name=\"Phase\" NumberOfComponents=\"1\" format=\"ascii\">\n');
p = 0;
for k=1:size(A.Phase,3)
    % Display progress
    if round(100*k/size(A.Phase,3))> p
        p= round(100*k/size(A.Phase,3));
        disp(['progress ' num2str(p) '%']);
    end
    for j=1:size(A.Phase,2)
        for i=1:size(A.Phase,1)
            fprintf(fid,'\t\t\t\t\t %ld\n', uint32(A.Phase(i,j,k)));
        end
    end
end

disp(['Writing Temp information... '])
fprintf(fid,'\t\t\t\t<DataArray type=\"Float64\" Name=\"Temp\" NumberOfComponents=\"1\" format=\"ascii\">\n');
for k=1:size(A.Phase,3)
    for j=1:size(A.Phase,2)
        for i=1:size(A.Phase,1)
            fprintf(fid,'\t\t\t\t\t %d\n', A.Temp(i,j,k));
        end
    end
end

fprintf(fid,'\t\t\t\t</DataArray>\n');
            
fprintf(fid,'\t\t\t</PointData>\n');
fprintf(fid,'\t\t</Piece>\n');
fprintf(fid,'\t</RectilinearGrid>\n');
fprintf(fid,'</VTKFile>\n');
end


%% =======================================================================
function a = WriteVTKModelSetup_BINARY(A)
% Write rectilinear file in BINARY format

% Definitions and initialization
sizeof_Float32  =   4;      
sizeof_Float64  =   4;     
sizeof_UInt64   =   8; 
sizeof_UInt32   =   4; 
sizeof_UInt16   =   2; 
sizeof_UInt8    =   1; 

offset = 0;

fname_vtk   = 'VTK_ModelSetup_paraview_binary.vtr';
fid         = fopen(fname_vtk,'w','b');           % 'b': BigEndian

fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="RectilinearGrid" version="0.1" byte_order="BigEndian" >\n');

fprintf(fid,'\t<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n',1,A.nump_x,1,A.nump_y,1,A.nump_z);
fprintf(fid,'\t\t<Piece Extent=\"%d %d %d %d %d %d\">\n',1,A.nump_x,1,A.nump_y,1,A.nump_z);

fprintf(fid,'\t\t\t<CellData></CellData>\n');

fprintf(fid,'\t\t\t<Coordinates>\n');

fprintf(fid,'\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_X\" NumberOfComponents=\"1\" format=\"appended\"  offset=\"%ld\"/>\n',int64(offset));
offset = offset + 1*sizeof_UInt32 + sizeof_Float32*length(A.x); 

fprintf(fid,'\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_Y\" NumberOfComponents=\"1\" format=\"appended\"  offset=\"%ld\"/>\n',int64(offset));
offset = offset + 1*sizeof_UInt32 + sizeof_Float32*length(A.y);

fprintf(fid,'\t\t\t\t<DataArray type=\"Float32\" Name=\"Coordinates_Z\" NumberOfComponents=\"1\" format=\"appended\"  offset=\"%ld\"/>\n',int64(offset));
offset = offset + 1*sizeof_UInt32 + sizeof_Float32*length(A.z);

fprintf(fid,'\t\t\t</Coordinates>\n');
fprintf(fid,'\n');

fprintf(fid,'\n');
fprintf(fid,'\t\t\t<PointData Scalars=\"Dimensional\">\n');
fprintf(fid,'\n');
% 
if  max(A.Phase(:))<256
    fprintf(fid,'\t\t\t\t<DataArray type=\"Int8\" Name=\"Phase\" NumberOfComponents=\"1\" format=\"appended\"  offset=\"%ld\"/>\n',int64(offset));
    offset = offset + 1*sizeof_UInt32 + sizeof_UInt8*(size(A.Phase,1)*size(A.Phase,2)*size(A.Phase,3));

else
    
    fprintf(fid,'\t\t\t\t<DataArray type=\"Int16\" Name=\"Phase\" NumberOfComponents=\"1\" format=\"appended\"  offset=\"%ld\"/>\n',int64(offset));
    offset = offset + 1*sizeof_UInt32 + sizeof_UInt16*(size(A.Phase,1)*size(A.Phase,2)*size(A.Phase,3));
end
fprintf(fid,'\t\t\t\t<DataArray type=\"Float32\" Name=\"Temp\" NumberOfComponents=\"1\" format=\"appended\"  offset=\"%ld\"/>\n',int64(offset));
offset = offset + 1*sizeof_UInt32 + sizeof_Float32*(size(A.Phase,1)*size(A.Phase,2)*size(A.Phase,3));
 

fprintf(fid,'\t\t\t</PointData>\n');
fprintf(fid,'\t\t</Piece>\n');
fprintf(fid,'\t</RectilinearGrid>\n');

%%% ---------------- Appended ---------------- data %%%
fprintf(fid,'  <AppendedData encoding=\"raw\">\n');
fprintf(fid,'_');

disp(['Writing Appended data for binary output in Paraview... '])
% X coord
fwrite(fid,int32(length(A.x)*sizeof_Float32),'int32');
fwrite(fid,A.x(:),'float32');

% Y coord
fwrite(fid,int32(length(A.y)*sizeof_Float32),'int32');
fwrite(fid,A.y(:),'float32');

% Z coord
fwrite(fid,int32(length(A.z)*sizeof_Float32),'int32');
fwrite(fid,A.z(:),'float32');

%% Properties
% Phase information - save as integer
% 
if  max(A.Phase(:))<256
    % we can safely save our data as 8 bit integers (implying that we can
    % have up to 256 phases)
    if size(A.Phase,1)*size(A.Phase,2)*size(A.Phase,3)*sizeof_UInt8>2^32
        warning('Problem with writing data to VTK file, as the number of entries is larger than 2^32')
    end
    
    fwrite(fid,int32(size(A.Phase,1)*size(A.Phase,2)*size(A.Phase,3)*sizeof_UInt8),'int32');
    fwrite(fid,int8(A.Phase(:)),'int8');
    
else
    % Use 16 bit, which would allow us to have 65535 phases, which seems
    % sufficient for most model setups
    fwrite(fid,int32(size(A.Phase,1)*size(A.Phase,2)*size(A.Phase,3)*sizeof_UInt16),'int32');
    fwrite(fid,int16(A.Phase(:)),'int16');

end

fwrite(fid,int32(size(A.Phase,1)*size(A.Phase,2)*size(A.Phase,3)*sizeof_Float32),'float32');
fwrite(fid,int8(A.Temp(:)),'float32');


fprintf(fid,'</VTKFile>\n');
end