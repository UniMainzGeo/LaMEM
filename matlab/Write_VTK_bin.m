%Write_VTK
%
% writes the grid into a VTK file, readable with e.g. paraview
%

%clear;
%load test

%fname = 'test';

% Create the names
VTK_Grid = ['Grid_',fname(1:end-4)];

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


%===================================================
% Write passsive tracer surface to file
%===================================================

% Write the Poly file in binary form
fid = fopen([VTK_Poly,'.vtk'],'w','b');                                    % note the 'b' for big endian!
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
fwrite(fid,Data_TRI.','int');                                                   % Write Binary data
fclose(fid);                  
%===================================================
% Finished writing passsive tracer surface to file
%===================================================


%===================================================
% Write grid info to file
%===================================================

% Write the Grid file in binary form
fid = fopen([VTK_Grid,'.vtk'],'w','b');                                    % note the 'b' for big endian!
fprintf(fid,'# vtk DataFile Version 3.0 \n');                                  % VTK format
fprintf(fid,'Unstructured Grid Data \n');                                  
fprintf(fid,'BINARY \n'); 
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');

% Write the points to the file
fprintf(fid,'POINTS %i float \n',size(Points,1));
Points = single(full(Points));
fwrite(fid,Points.','single');                                       % Write binary data

% Describe the cells
fprintf(fid,'\nCELLS %i %i \n',nel,nel*(nnel+1));
Cells = int32(Cells);
fwrite(fid,Cells.','int');                                                   % Write Binary data

% Describe celltypes
Celltypes = ones(1,nel,'int32')*12;  % for 3D brick elements
fprintf(fid,'CELL_TYPES %i \n',nel);
fwrite(fid,Celltypes,'int');                                                   % Write Binary data


% Write Cell data
fprintf(fid,'\n');
fprintf(fid,'CELL_DATA %i \n',nel);

% Describe density
fprintf(fid,'SCALARS Density float\n');
fprintf(fid,'LOOKUP_TABLE default \n');
Density_vec = single(Density_vec(:));
fwrite(fid,Density_vec.','single');

% Viscosity
fprintf(fid,'SCALARS Viscosity float\n');
fprintf(fid,'LOOKUP_TABLE default \n');
Viscosity_vec = single(Viscosity_vec(:));
fwrite(fid,Viscosity_vec.','single');

% Write scalar data
fprintf(fid,'\n');
fprintf(fid,'POINT_DATA %i \n',size(Points,1));

% Describe Scalar data
fprintf(fid,'SCALARS Vx float\n');
fprintf(fid,'LOOKUP_TABLE default \n');
Vx_Vec = single(Velocity_Vec(:,1));
fwrite(fid,Vx_Vec.','single'); 

% Describe Scalar data
fprintf(fid,'SCALARS Vy float\n');
fprintf(fid,'LOOKUP_TABLE default \n');
Vy_Vec = single(Velocity_Vec(:,2));
fwrite(fid,Vy_Vec.','single'); 

% Describe Scalar data
fprintf(fid,'SCALARS Vz float\n');
fprintf(fid,'LOOKUP_TABLE default \n');
Vz_Vec = single(Velocity_Vec(:,3));
fwrite(fid,Vz_Vec.','single'); 

% Describe Vector data
fprintf(fid,'VECTORS velocity float\n');
Velocity_Vec = single(Velocity_Vec);
fwrite(fid,Velocity_Vec.','single'); 


fclose(fid);                  
%===================================================
% Finished writing passsive tracer surface to file
%===================================================


disp(['Finished writing VTK files: ', [VTK_Grid,'.vtu'],' and ' [VTK_Poly,'.vtp']])

