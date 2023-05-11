# 3.6 How to use a topography input file

A topography file can be used to introduce an **initial topography** into a LaMEM model. This can be set up as a 2D matrix in Matlab. The matrix has to cover the entire LaMEM box, or can be bigger. The number of points in the matrix does not matter but it's dimensions have to be included in the file as well as the coordinate of the SW corner and the gridspacing. Make sure that the units fit the LaMEM length unit (km for geo units, m for SI units).

```
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
```

![TopoFromFile](./Pictures/TopoFromFile.png)