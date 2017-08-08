%% Generate parallel constraint cells flags from 3D raw image data

clc, clear all , close all

%% Add important paths

addpath('/home/anton/PROG/LaMEM/matlab')

%% Defining Parameters
mesh_file      = 'ProcessorPartitioning_16cpu_2.2.4.bin';
pixel_file     = 'voxel3D-512x512x512.raw';
npx            = 512;        % x-pixel grid size
npy            = 512;        % y-pixel grid size
npz            = 512;        % z-pixel grid size
rock_threshold = 0.5;        % upscaling threshold
SizeOfVoxel    = 3.528e-6;   % pixel physical resolution
Is64BIT        = 0;

%==========================================================================

%% LOAD DATA

% Read raw iamge data file
fp  = fopen(pixel_file, 'r');
P3D = fread(fp); 
P3D = reshape(P3D, npx, npy, npz);

% Create pixel coordinates
lx = SizeOfVoxel*(npx-1)/2;
ly = SizeOfVoxel*(npy-1)/2;
lz = SizeOfVoxel*(npz-1)/2;
X  = linspace(-lx, lx, npx);
Y  = linspace(-ly, ly, npy);
Z  = linspace(-lz, lz, npz);

% Read processor partitioning
[P] = GetProcessorPartitioning(mesh_file, Is64BIT);

% get grid data
nx     = P.nnodx-1;
ny     = P.nnody-1;
nz     = P.nnodz-1;
xcoor  = P.xcoor;
ycoor  = P.ycoor;
zcoor  = P.zcoor;
ncell  = nx*ny*nz;

% check whether resolution is too high (silly)
if(nx > npx || ny > npy || nz > npz)
    error('Grid resolution cannot exceed pixel resolution');
end

% check whether upscaling is necessary (most likely)
if(nx ~= npx || ny ~= npy || nz ~= npz)
   
    %% PERFORM UPSCALING
    
    % map pixels on computational grid
    NX = histcounts(X, xcoor);
    NY = histcounts(Y, ycoor);
    NZ = histcounts(Z, zcoor);

    % compute number of pixels in each cell
    N3D = repmat (NX'*NY, [1, 1, numel(NZ)]);
    NZR = reshape(NZ,     [1, 1, numel(NZ)]);
    N3D = N3D.*NZR;
    
    % check for empty cell (just in case)
    if(min(N3D(:)) == 0)
        error('Empty cell');    
    end
        
    % create matrix of host cell ID for each pixel
    cellID = reshape(1:ncell, nx, ny, nz);
    cellID = repelem(cellID,  NX, NY, NZ);

    % find all rock pixels
    pind = find(P3D == 1);
    
    % get cells IDs containing rock pixels 
    cind = cellID(pind);

    % count rock pixels in each cell
    P3D = histcounts(cind, 0.5:1:(ncell + 0.5));
    P3D = reshape(P3D, nx, ny, nz);
    
    % get rock phase ratio
    P3D = P3D./N3D;
    
    % round-off
    P3D(P3D <  rock_threshold) = 0;
    P3D(P3D >= rock_threshold) = 1;
       
end

clearvars -except P3D P

%% STORE FILES

if ~isdir('bc')
    mkdir bc
end

% get partitioning data
Nprocx = P.Nprocx;
Nprocy = P.Nprocy;
Nprocz = P.Nprocz;
ix     = P.ix;
iy     = P.iy;
iz     = P.iz;

% initialize subdomain ID
num = 0;

for k = 1:Nprocz
    for j = 1:Nprocy
        for i = 1:Nprocx
            
            % get cell subset of current subdomain
            P3DS = P3D(ix(i):ix(i+1)-1, iy(j):iy(j+1)-1, iz(k):iz(k+1)-1);
            
            % open file
            fname = sprintf('./bc/cdb.%1.8d.dat', num);
            
            disp(['Writing file -> ', fname])

            fp = fopen(fname, 'w');
            
            % dump data on disk
            fwrite(fp, P3DS, 'uint8');

            fclose(fp);

            clear P3DS fname fp
     
            % update processor ID
            num = num + 1;
            
        end
    end
end


