function []=plot_mesh_3D(X3d,Y3d,Z3d,GridType, x,y,z, color)
% Functions plots the 3D mesh
%

if length(x)>0
    if length(x)>size(X3d,2)
        error('x dimension too large!')
    end
    x2d = squeeze(X3d(:,x,:));
    y2d = squeeze(Y3d(:,x,:));
    z2d = squeeze(Z3d(:,x,:));
end
if length(y)>0
    if length(y)>size(X3d,1)
        error('x dimension too large!')
    end
    x2d = squeeze(X3d(y,:,:));
    y2d = squeeze(Y3d(y,:,:));
    z2d = squeeze(Z3d(y,:,:));
end
if length(z)>0
     if length(z)>size(X3d,3)
        error('z dimension too large!')
    end
    x2d = squeeze(X3d(:,:,z));
    y2d = squeeze(Y3d(:,:,z));
    z2d = squeeze(Z3d(:,:,z));
end


switch GridType
 case '3D_Lin_rectangular'
  plot3(x2d,y2d,z2d,color,  x2d.',y2d.',z2d.',color), hold on, view(3)
 case '3D_Quad_rectangular'
  
  step = 2;
  plot3(x2d(1:step:end,1:step:end)  ,y2d(1:step:end,1:step:end)  ,z2d(1:step:end,1:step:end),color,...
        x2d(1:step:end,1:step:end).',y2d(1:step:end,1:step:end).', ...
        z2d(1:step:end,1:step:end).',color) 
  hold on, view(3)
 % plot3(x2d,y2d,z2d,'k.')
end


