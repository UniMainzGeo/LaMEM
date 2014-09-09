%Matlab2VTK
%
% Writes a VTK file from 3D data that has been collected with matlab (from
% several CPUs) for easier 3D visualization with Paraview



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
nnodes          =   nnode_x*nnode_y*nnode_z;
Velocity_vec    =   zeros(nnodes,3);
Points          =   zeros(nnodes,3);

Points          = [X(:),        Y(:),       Z(:)];
Velocity_vec    = [Vx_3D(:),    Vy_3D(:),   Vz_3D(:)];
Velocity_vec    = Velocity_vec*characteristic.cmYear;

Temperature_vec = T_3D(:);
CPU_vec         = CPU_3D(:);




%==========================================================================
% Write the Mesh-based data in ascii
%==========================================================================
fid = fopen([fname,'.',num2str(time_step+100000),'.MeshData.vtk'],'w');

fprintf(fid,'# vtk DataFile Version 2.0 \n');
fprintf(fid,'Structured Grid Example \n');
fprintf(fid,'ASCII \n \n');
fprintf(fid,'DATASET STRUCTURED_GRID\n');
fprintf(fid,['DIMENSIONS ',num2str(size(X,1)),' ',num2str(size(X,2)),' ',num2str(size(X,3)),' \n' ] );


% Write the points to the file
num_points = nnode_x*nnode_y*nnode_z;
fprintf(fid,'POINTS %i double \n',num_points);
for i=1:num_points
    fprintf(fid,'%g %g %g \n',Points(i,1),Points(i,2),Points(i,3));
end

% Write the some pointwise data to the file
fprintf(fid,['POINT_DATA ',num2str(num_points),' \n']);

fprintf(fid,['SCALARS Temperature[C] float 1 \n']);
fprintf(fid,['LOOKUP_TABLE default \n']);
fprintf(fid,'%g \n',Temperature_vec(:)-273);    

fprintf(fid,['SCALARS Processor int 1 \n']);
fprintf(fid,['LOOKUP_TABLE default \n']);
fprintf(fid,'%g \n',int32(CPU_3D(:)));


% Write the velocity to the file
fprintf(fid,'VECTORS velocity[cm/yr] double \n');
% fprintf(fid,['LOOKUP_TABLE default \n']);
for i=1:num_points
    fprintf(fid,'%g %g %g \n',Velocity_vec(i,1),Velocity_vec(i,2),Velocity_vec(i,3));    
end

fclose(fid)



