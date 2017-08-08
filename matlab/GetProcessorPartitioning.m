function [P] = GetProcessorPartitioning(filename, Is64BIT)

%========================================================
% Read processor partitioning into data structure
%
%   P - structure contains:
%
%      Nprocx - number of processors
%      Nprocy
%      Nprocz
%      nnodx  - number of nodes
%      nnody
%      nnodz
%      ix     - first node index on each processor (last entry - number of nodes)
%      iy
%      iz
%      xcoor  - nodal coordinates
%      ycoor
%      zcoor
%      xc     - processor coordinate bounds
%      yc
%      zc
%========================================================

% open file
fid = PetscOpenFile(filename);

if Is64BIT
    % In case file was written by 64 BIT compiled PETSC version
    Precision_INT       = 'int64';
    Precision_SCALAR    = 'float64';
else
    Precision_INT       = 'int32';
    Precision_SCALAR    = 'double';
end

% get number of processors
P.Nprocx = read(fid, 1, Precision_INT);
P.Nprocy = read(fid, 1, Precision_INT);
P.Nprocz = read(fid, 1, Precision_INT);

% get number of nodes
P.nnodx = read(fid, 1, Precision_INT);
P.nnody = read(fid, 1, Precision_INT);
P.nnodz = read(fid, 1, Precision_INT);

% get first node indices
P.ix = read(fid, P.Nprocx + 1, Precision_INT);
P.iy = read(fid, P.Nprocy + 1, Precision_INT);
P.iz = read(fid, P.Nprocz + 1, Precision_INT);

% get coordinate scaling
CharLength = read(fid, 1, Precision_SCALAR);

% get nodal coordinates
P.xcoor = read(fid, P.nnodx, Precision_SCALAR);
P.ycoor = read(fid, P.nnody, Precision_SCALAR);
P.zcoor = read(fid, P.nnodz, Precision_SCALAR);

close(fid);

% shift indices (no difference for last entry, it's a last node C-index)
P.ix = P.ix + 1;
P.iy = P.iy + 1;
P.iz = P.iz + 1;

% scale coordinates
P.xcoor = P.xcoor*CharLength;
P.ycoor = P.ycoor*CharLength;
P.zcoor = P.zcoor*CharLength;

% processor coordinate bounds
P.xc = P.xcoor(P.ix);
P.yc = P.ycoor(P.iy);
P.zc = P.zcoor(P.iz);

end

