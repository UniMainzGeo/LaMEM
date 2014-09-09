function [coord_x_3D,coord_y_3D,coord_z_3D,info] = ReadInitialMesh(Directory)
% ReadInitialMesh
%
% This routine reads the initial mesh that's created by a LaMEM simulation
% 
% Note: the initial mesh MUST be located in the directory:
%           ./InitialMesh
%
%
% $Id$ Boris Kaus

if nargin==0
    Directory   =   'InitialMesh';
end



addpath('./LaMEM_MATLAB');  % add I/O routines



[info]   =   PetscBinaryRead(['./',Directory,'/InitialMesh.0.out']);
num_proc =   info(4);  nnode_x  =   info(1);    nnode_y  =   info(2); nnode_z  =   info(3);

coord_x_3D = zeros(nnode_x,nnode_y,nnode_z);    coord_y_3D = zeros(nnode_x,nnode_y,nnode_z);   coord_z_3D = zeros(nnode_x,nnode_y,nnode_z);
for i=1:num_proc  % read files from all CPU's

    %newest version
    [info,coord_x,coord_y,coord_z]  =   PetscBinaryRead(['./',Directory,'/InitialMesh.',num2str(i-1),'.out']);

    % (1) Deal with nodal based data
    xs = info(5);    ys = info(6);    zs = info(7);    xm = info(8);    ym = info(9);    zm = info(10);
    char.km = info(12); char.Length	=info(11);
    
    % Shape local data into 3D arrays
    coord_x_3D_loc  = reshape(coord_x,[xm ym zm]);   coord_y_3D_loc = reshape(coord_y,[xm ym zm]);   coord_z_3D_loc = reshape(coord_z,  [xm ym zm]);

    % Put local 3D arrays in global 3D arrays
    coord_x_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_x_3D_loc;
    coord_y_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_y_3D_loc;
    coord_z_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_z_3D_loc;
   
end

disp(['Read the initial mesh from directory ./',Directory,' that was generated on ',num2str(num_proc),' procs'])

X = coord_x_3D;
Y = coord_y_3D;
Z = coord_z_3D;





