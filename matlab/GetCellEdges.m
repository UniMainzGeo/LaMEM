% --------------------------------------
function [xcoor,ycoor,zcoor] = GetCellEdges(test, Is64BIT)

% Read Processor Partitioning file
fid=PetscOpenFile(test);

if Is64BIT
    % In case file was written  by 64 BIT compiled PETSC version
    Precision_INT       = 'int64';
    Precision_SCALAR    = 'float64';
else
    Precision_INT       = 'int32';
    Precision_SCALAR    = 'double';
end

Nprocx=read(fid,1,Precision_INT);
Nprocy=read(fid,1,Precision_INT);
Nprocz=read(fid,1,Precision_INT);

nnodx=read(fid,1,Precision_INT);
nnody=read(fid,1,Precision_INT);
nnodz=read(fid,1,Precision_INT);

read(fid,Nprocx+1,Precision_INT);
read(fid,Nprocy+1,Precision_INT);
read(fid,Nprocz+1,Precision_INT);

CharLength=read(fid,1,Precision_SCALAR);

xcoor=read(fid,nnodx,Precision_SCALAR);
ycoor=read(fid,nnody,Precision_SCALAR);
zcoor=read(fid,nnodz,Precision_SCALAR);

close(fid);

% Dimensionalize
xcoor = xcoor*CharLength;
ycoor = ycoor*CharLength;
zcoor = zcoor*CharLength;

end