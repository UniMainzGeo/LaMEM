% 3D 
% =========================================================================

% clean up
clear
close all

% settings
inputFile           = ['Simple.EW.svg'];
opt                 = geomIO_Options();
opt.outputDir       = ['./Output'];
opt.inputFileName   = inputFile;
opt.LaMEMinputFileName ='Simple.dat';
opt.readLaMEM       = true;
opt.writeParaview   = true;
opt.writePolygons   = true;
opt.interp = true;
opt.zi = [0:5:200];
opt.getPolygons= true;
opt.gravity.lenUnit = 'm';

% Density assignment
paths = {
    'Air', 0, 0
    'Sediments', 2500, 1
    'Crust', 2700, 2
    'Magma', 2400, 3
     };
opt.pathNames       = {paths{:,1}}; 
opt.gravity.drho    = [paths{:,2}]; 
opt.phase           = [paths{:,3}];    

% Run geomIO
[PathCoord,Volumes,opt] = run_geomIO(opt,'default');

%% Setup temperature
% get bounds of LaMEM box
x_box       = opt.setup.coord_x;
y_box       = opt.setup.coord_y;
z_box       = opt.setup.coord_z;

% fit marker distribution
nx  = length(opt.setup.marker_x);
ny  = length(opt.setup.marker_y);
nz  = length(opt.setup.marker_z);

% initialize matrix
T3D = zeros(nx,ny,nz);

% set up coordinate grid (for assigning specific temperatures based on coordinates)
dx      = (x_box(2)-x_box(1))/nx;
dy      = (y_box(2)-y_box(1))/ny; 
dz      = (z_box(2)-z_box(1))/nz;
x_vec   = (x_box(1)+dx/2:dx:x_box(2)-dx/2);
y_vec   = (y_box(1)+dy/2:dy:y_box(2)-dy/2);
z_vec   = (z_box(1)+dz/2:dz:z_box(2)-dz/2);
% DO NOT USE MESHGRID HERE!
[X,Y,Z] = ndgrid(x_vec,y_vec,z_vec);

% setup linear Temperature profile starting at the surface
z_surf      = 3;
gradient    = 30;      % [K/km]
zeros(nz,1);
ind_surf    = find(z_vec < z_surf);
T_vec       = zeros(nz,1);
T_vec(ind_surf)     = (z_surf - z_vec(ind_surf)) * gradient;

for i = 1 : nz
    T3D(:,:,i) = T_vec(i);
end

% assign Temperature to volume
vol     = Volumes.Magma;
T_vol   = 1180;

% run through all depths
for k = 1 : nz
    depth   = z_vec(k);
    Slices  = {};
    % find slices at that depth
    for s = 1 : length(vol.Polygons)
        if vol.Polygons{s}(1,3) <= depth + dz/10 && vol.Polygons{s}(1,3) >= depth - dz/10
            Slices = [Slices; vol.Polygons{s}];
        end
    end
    
    % if slices of the body exist at that depth
    if ~isempty(Slices)
        allInds = false(nx,ny);
        % add all markers that lie inside the body
        for s = 1 : length(Slices)
            inds = inpolygon(X(:,:,k),Y(:,:,k),Slices{s}(:,1),Slices{s}(:,2));
            allInds = allInds | inds;
        end
        % assign the markers a temperature
        T_Slice             = squeeze(T3D(:,:,k));
        T_Slice(allInds)    = T_vol; 
        T3D(:,:,k) = T_Slice;
    end
    
end

% write binary to be read by LaMEM
%(Filename,[# of points in x-dir, # of points in y-dir, # of points in z-dir, Temp-Matrix in vector form])
petscBinaryWrite('T3D.dat',[nx;ny;nz;T3D(:)]);

%% Setup topography

% import your topography data
% Create example data:
Topo        = peaks(301);
Topo        = Topo/4;
Easting     = (0:1:300);
Northing    = (0:1:300);

% For orientation
% Topo(1,1):     SW
% Topo(1,end):   SE
% Topo(end,1):   NW
% Topo(end,end): NE

% compute grid spacing
dx = (max(Easting) - min(Easting))/(length(Easting)-1);
dy = (max(Northing) - min(Northing))/(length(Northing)-1);

% possibly add offset in data
Easting     = Easting - 200;
Northing    = Northing -100;

% transpose Topo to be read correctly by LaMEM
Topo    = Topo';

% write binary to be read by LaMEM
% (FILENAME,[# of points in x-dir, # of points in y-dir, x-coord of SW corner, y_coord of SW corner,
% grid spacing in x-dir, grid spacing in y-dir, TopoMatrix in vector form])
petscBinaryWrite('Topo.dat', [size(Topo,1); size(Topo,2); min(Easting);min(Northing); dx; dy; Topo(:)]);