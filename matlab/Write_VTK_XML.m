%Write_VTK_XML
%
% writes the grid into a VTK-XML file
%


disp(['Creating VTK files'])

%clear;
%$load test

%fname = 'test';

% Create the names
VTK_Grid = ['Grid',fname(1:end-4)];
VTK_Poly = ['Poly',fname(1:end-4)];



nnode = 0;
NumNodes = zeros(size(X));
for iz=1:size(X,3)
    for ix=1:size(X,2)
     for iy=1:size(X,1)
            NumNodes(iy,ix,iz) = nnode;
            nnode=nnode+1;
        end
    end
end
num_nodes = prod(size(X));


%==========================================================================
% Prepare the data for writing to VTK file
%==========================================================================
nel     =   nel_x*nel_y*nel_z;
nnel    =   8;
%Points  =   zeros(nel*nnel,3);
num_el  =   1;
Cells   =   zeros(nel,nnel+1);
Vx_vec  =   zeros(nel*8,1);

for iel_x=1:nel_x;
    for iel_y=1:nel_y;
        for iel_z=1:nel_z;
            
            % Retrieve indices
            ind_local = sub2ind(size(X),...
                [iel_x iel_x+1 iel_x+1 iel_x   iel_x   iel_x+1 iel_x+1 iel_x  ],...        
                [iel_y iel_y   iel_y+1 iel_y+1 iel_y   iel_y   iel_y+1 iel_y+1],...
                [iel_z iel_z   iel_z   iel_z   iel_z+1 iel_z+1 iel_z+1 iel_z+1]);
          
            %Points( ((num_el-1)*nnel+1):num_el*nnel,1)  = X(ind_local);
            %Points( ((num_el-1)*nnel+1):num_el*nnel,2)  = Y(ind_local);
            %Points( ((num_el-1)*nnel+1):num_el*nnel,3)  = Z(ind_local);
            
            
            %Cells(num_el,:) = [8, [((num_el-1)*nnel+1):num_el*nnel]-1 ];
            Cells(num_el,:) = [8, NumNodes(ind_local)];
            
            % Store velocity
            %Velocity_vec( ((num_el-1)*nnel+1):num_el*nnel,1)  = Vx(ind_local);
            %Velocity_vec( ((num_el-1)*nnel+1):num_el*nnel,2)  = Vy(ind_local);
            %Velocity_vec( ((num_el-1)*nnel+1):num_el*nnel,3)  = Vz(ind_local);
            
            % Density, Viscosity
            %Density_vec(((num_el-1)*nnel+1):num_el*nnel,1)     =    Rho(iel_x,iel_y,iel_z);
            %Viscosity_vec(((num_el-1)*nnel+1):num_el*nnel,1)   =    Mu(iel_x,iel_y,iel_z);
            Density_vec( num_el,1)       =    Rho(iel_x,iel_y,iel_z);
            Viscosity_vec( num_el,1)     =    Mu(iel_x,iel_y,iel_z);
            
            
            num_el  =   num_el+1;
        end
    end
end

Velocity_Vec = [];
Velocity_Vec(:,1) = Vx(:);
Velocity_Vec(:,2) = Vy(:);
Velocity_Vec(:,3) = Vz(:);

Points = [];
Points(:,1)         =   X(:);
Points(:,2)         =   Y(:);
Points(:,3)         =   Z(:);



%nel = 1;
%==========================================================================
% End of preparing the data
%==========================================================================



%==========================================================================
%
% START WRITING VTU FILE
%
%==========================================================================

% header
fid = fopen([VTK_Grid,'.vtu'],'w');
fprintf(fid,'<?xml version="1.0"?> \n');
%fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n');
fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" >\n');
fprintf(fid,'  <UnstructuredGrid> \n');
fprintf(fid,'    <Piece NumberOfPoints="%i"  NumberOfCells="%i" >\n',num_nodes, nel);


% Point-wise data==========================================================
fprintf(fid,'      <PointData Vectors="Velocity" Scalars="Vx">\n');
fprintf(fid,'           <DataArray type="Float32"  Name="Velocity" NumberOfComponents="3"  format="ascii"> \n');
fprintf(fid,'               ');
offset = nnel;
for i=1:num_nodes
    for k=1:3
        fprintf(fid,' %g',Velocity_Vec(i,k));
    end
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');

fprintf(fid,'           <DataArray type="Float32"  Name="Vx"  format="ascii"> \n');
fprintf(fid,'               ');
offset = nnel;
for i=1:num_nodes
        fprintf(fid,' %g',Velocity_Vec(i,1));
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');


fprintf(fid,'           <DataArray type="Float32"  Name="Vy"  format="ascii"> \n');
fprintf(fid,'               ');
offset = nnel;
for i=1:num_nodes
        fprintf(fid,' %g',Velocity_Vec(i,2));
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');

fprintf(fid,'           <DataArray type="Float32"  Name="Vz"  format="ascii"> \n');
fprintf(fid,'               ');
offset = nnel;
for i=1:num_nodes
        fprintf(fid,' %g',Velocity_Vec(i,3));
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');

fprintf(fid,'      </PointData>\n');
%==========================================================================

% Cell-wise data===========================================================
fprintf(fid,'      <CellData Scalar="Density">\n');

fprintf(fid,'           <DataArray type="Float32"  Name="Density" format="ascii"> \n');
fprintf(fid,'               ');
offset = nnel;
for i=1:nel
        fprintf(fid,' %g',Density_vec(i));
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');

fprintf(fid,'           <DataArray type="Float32"  Name="Viscosity" format="ascii"> \n');
fprintf(fid,'               ');
offset = nnel;
for i=1:nel
        fprintf(fid,' %g',Viscosity_vec(i));
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');


fprintf(fid,'      </CellData>\n');
%==========================================================================

% Nodal coordinates =======================================================
fprintf(fid,'        <Points>\n');
fprintf(fid,'           <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
%fprintf(fid,'           <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="0">\n');
fprintf(fid,'               ');
for ipoint=1:num_nodes
   for k=1:3
       fprintf(fid,' %g',Points(ipoint,k));
   end
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');
fprintf(fid,'        </Points> \n');
%==========================================================================

% describe the cells=======================================================
fprintf(fid,'        <Cells> \n');

fprintf(fid,'           <DataArray type="Int32" Name="connectivity" format="ascii"> \n');
fprintf(fid,'               ');
for iel=1:nel
    for inode=2:nnel+1
        fprintf(fid,' %i',Cells(iel,inode));
    end
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');

fprintf(fid,'           <DataArray type="Int32" Name="offsets" format="ascii"> \n');
fprintf(fid,'               ');
offset = nnel;
for iel=1:nel
        fprintf(fid,' %i',offset);
        offset = offset+nnel;
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');

fprintf(fid,'           <DataArray type="UInt8" Name="types" format="ascii"> \n');
fprintf(fid,'               ');
types = 12;
for iel=1:nel
        fprintf(fid,' %i',types);
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');
fprintf(fid,'        </Cells> \n');
%==========================================================================


fprintf(fid,'     </Piece>\n');    
fprintf(fid,'  </UnstructuredGrid> \n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);


%==========================================================================
%
% FINISHED WRITING VTU FILE
%
%==========================================================================





%==========================================================================
%
% START WRITING VTP FILE
%
%==========================================================================

% header
fid = fopen([VTK_Poly,'.vtp'],'w');
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" >\n');
fprintf(fid,'  <PolyData> \n');
fprintf(fid,'    <Piece NumberOfPoints="%i"  NumberOfVerts="%i" NumberOfLines="%i" \n',size(PassiveParticles,1), 0 ,0);
fprintf(fid,'           NumberOfStrips="%i"  NumberOfPolys="%i">  \n',0,size(TRI,1));

% Point-wise data==========================================================
fprintf(fid,'      <PointData>\n');
fprintf(fid,'      </PointData>\n');
%==========================================================================

% Cell-wise data===========================================================
fprintf(fid,'      <CellData>\n');
fprintf(fid,'      </CellData>\n');
%==========================================================================

% Nodal coordinates =======================================================
fprintf(fid,'        <Points>\n');
fprintf(fid,'           <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n');
fprintf(fid,'               ');
for ipoint=1:size(PassiveParticles,1)
   for k=1:3
       fprintf(fid,' %g',PassiveParticles(ipoint,k));
   end
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');
fprintf(fid,'        </Points> \n');
%==========================================================================

% Polys =======================================================
fprintf(fid,'        <Polys>\n');

% Connectivity
fprintf(fid,'           <DataArray type="Int32" Name="connectivity" format="ascii">\n');
fprintf(fid,'               ');
for itri=1:size(TRI,1)
   for k=1:3
       fprintf(fid,' %i',TRI(itri,k)-1);
   end
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');

fprintf(fid,'           <DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(fid,'               ');
offset = 3;
for itri=1:size(TRI,1)
   fprintf(fid,' %i',offset);
   offset = offset+3;
end
fprintf(fid,'\n');
fprintf(fid,'           </DataArray> \n');


fprintf(fid,'        </Polys> \n');
%==========================================================================

fprintf(fid,'     </Piece>\n');    
fprintf(fid,'  </PolyData> \n');
fprintf(fid,'</VTKFile>\n');
fclose(fid);

%==========================================================================
%
% FINISHED WRITING VTP FILE
%
%==========================================================================

disp(['Finished writing VTK files: ', [VTK_Grid,'.vtu'],' and ' [VTK_Poly,'.vtp']])



% Write the Poly file in binary form
% Write points
fid = fopen([fname,'VTK_Poly.vtk'],'w','b');                                    % note the 'b' for big endian!
fprintf(fid,'# vtk DataFile Version 3.0 \n');                                  % VTK format
fprintf(fid,'Passive tracer surface \n');                                  
fprintf(fid,'BINARY \n'); 
fprintf(fid,'DATASET POLYDATA\n');

% Write points
PassiveParticles = single(full(PassiveParticles));                             % ensure that it has correct dataformat
fprintf(fid,'POINTS %i float \n',size(PassiveParticles,1));                    % # of points
fwrite(fid,PassiveParticles.','single');                                       % Write binary data

% Write triangle-edgepoints
Data_TRI = [ones(size(TRI,1),1,'int32')*3, int32(TRI-1)];                      % should be integer*32
fprintf(fid,'\nTRIANGLE_STRIPS %i %i \n',size(TRI,1),size(TRI,1)*4 );          % description of triangles
fwrite(fid,Data_TRI.','int')                                                   % Write Binary data

fclose(fid)                  


