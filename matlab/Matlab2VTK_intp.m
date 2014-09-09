%Matlab2VTK_intp
%
% Writes a VTK file from 3D data that has been collected with matlab (from
% several CPUs) for easier 3D visualization with Paraview



nnode = 0;
NumNodes = zeros(size(X));
for iz=1:size(intpx_3D,3)
    for ix=1:size(intpx_3D,2)
     for iy=1:size(intpx_3D,1)
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
nnodes          =   prod(size(intpx_3D));
Points          =   zeros(nnodes,3);

% intpx_3D=shiftdim(intpx_3D,2);
% intpx_3D=shiftdim(intpy_3D,2);
% intpz_3D=shiftdim(intpz_3D,2);


% Prepare input for tensor
Tensor_row1             =   [Txx_3D(:), Txy_3D(:), Txz_3D(:)];
Tensor_row2             =   [Txy_3D(:), Tyy_3D(:), Tyz_3D(:)];
Tensor_row3             =   [Txz_3D(:), Tyz_3D(:), Tzz_3D(:)];
Tensor                  =   zeros(size(Tensor_row1,1)*3,3);
Tensor(1:3:end-2,:)     =   Tensor_row1;
Tensor(2:3:end-1,:)     =   Tensor_row2;
Tensor(3:3:end  ,:)     =   Tensor_row3;
Tensor                  =   single(Tensor*characteristic.Stress/characteristic.MPa);        % single precision, in MPa 






Points          = [intpx_3D(:),        intpy_3D(:),       intpz_3D(:)];




%==========================================================================
% Write the Mesh-based data in ascii
%==========================================================================
fid = fopen([fname,'.',num2str(time_step+100000),'.MeshData_intp.vtk'],'w');

fprintf(fid,'# vtk DataFile Version 2.0 \n');
fprintf(fid,'Structured Grid Example \n');
fprintf(fid,'ASCII \n \n');
fprintf(fid,'DATASET STRUCTURED_GRID\n');
fprintf(fid,['DIMENSIONS ',num2str(size(intpx_3D,1)),' ',num2str(size(intpx_3D,2)),' ',num2str(size(intpx_3D,3)),' \n' ] );


% Write the points to the file
num_points = prod(size(intpx_3D));
fprintf(fid,'POINTS %i float \n',num_points);
for i=1:num_points
    fprintf(fid,'%g %g %g \n',Points(i,1),Points(i,2),Points(i,3));
end

% Write the some pointwise data to the file
fprintf(fid,['POINT_DATA ',num2str(num_points),' \n']);


 fprintf(fid,['SCALARS Density float 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',rho_3D(:));

 fprintf(fid,['SCALARS Viscosity float 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',mu_3D(:)); 
 
 fprintf(fid,['SCALARS Phases int 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',int32(Phases_3D(:)));
 
 fprintf(fid,['SCALARS T2nd[MPa] float 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',Tau2nd_3D(:)*characteristic.Stress/characteristic.MPa);
 
 fprintf(fid,['SCALARS log10(E2nd)[1/s] float 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',log10(E2nd_3D(:)*(1/characteristic.Time)));
 
 
 fprintf(fid,['SCALARS Pressure[MPa] float 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',Pressure_3D(:)*characteristic.Stress/characteristic.MPa);
 
 
 fprintf(fid,['SCALARS Strain float 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',single(Strain_3D(:)));

 
 fprintf(fid,['SCALARS PlasticStrain float 1 \n']);
 fprintf(fid,['LOOKUP_TABLE default \n']);
 fprintf(fid,'%g \n',single(PlasticStrain_3D(:)));
 
 
 fprintf(fid,['TENSORS DeviatoricStressTensor[MPa] float\n']);
 for i=1:size(Tensor,1)
    fprintf(fid,'%g %g %g \n',Tensor(i,:));    
 end

 
fclose(fid)



