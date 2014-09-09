function [] = WriteInitialMesh(X,Y,Z,info,InitialMesh_Filename)
% WriteInitialMesh
%
% This routine writes a FEM mesh to file such that ir can be read in by
% LaMEM
% 
% $Id$ 


addpath('./LaMEM_MATLAB');  % add I/O routines


% Initialize---------------------------------------------------------------
nnode_x = size(X,1);
nnode_y = size(X,2);
nnode_z = size(X,3);

coord_x_vec = zeros(nnode_x*nnode_y*nnode_z,1);
coord_y_vec = zeros(nnode_x*nnode_y*nnode_z,1);
coord_z_vec = zeros(nnode_x*nnode_y*nnode_z,1);

num = 1;
for ix=1:nnode_x
    for iy=1:nnode_y
        for iz=1:nnode_z
            coord_x_vec(num) = X(ix,iy,iz);
            coord_y_vec(num) = Y(ix,iy,iz);
            coord_z_vec(num) = Z(ix,iy,iz);
            num=num+1;
        end
    end
end
%--------------------------------------------------------------------------

% Save mesh ---------------------------------------------------------------
PetscBinaryWrite(InitialMesh_Filename,info, coord_x_vec, coord_y_vec, coord_z_vec);
%--------------------------------------------------------------------------


disp(['Wrote initial mesh file to ',InitialMesh_Filename])
disp(['Ensure that this is present in your LaMEM input file!'])
disp(['InitialMeshFromFile 	=   1'])
disp(['InitialMeshFileName 	=   ',InitialMesh_Filename])
