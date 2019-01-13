% Create a 2D subduction setup with particles & temperature, 


clear, clf, close all


% Tell the code where the LaMEM matlab routines are 
addpath ../../matlab


% Define parallel partition file
Parallel_partition     = 'ProcessorPartitioning_4cpu_4.1.1.bin'


% Number of markers in a grid cell [Note that this needs to be changed if 
% this is changed in the *.dat file!]
npart_x = 3;
npart_y = 3;
npart_z = 3;

%==========================================================================
% OUTPUT OPTIONS (standard)
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        =    1;

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output  =    1;

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh               =    1;

% random noise of particles
RandomNoise             =   logical(0);

Is64BIT                 =   logical(0); % only used for some 

% Load grid from parallel partitioning file
[X,Y,Z,xcoor,ycoor,zcoor,Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition,RandomNoise,Is64BIT);

% Update other variables
nump_x  = size(X,2);
nump_y  = size(X,1);
nump_z  = size(X,3);



% Domain parameters
W       =   xcoor(end)-xcoor(1);    % x-dir
L       =   ycoor(end)-ycoor(1);    % y-dir
H       =   zcoor(end)-zcoor(1);    % z-dir

%==========================================================================


%==========================================================================
% SPECIFY PARAMETERS OF THE SLAB
%==========================================================================
Trench_x_location   = -500;     % trench location
Length_Subduct_Slab =  300;     % length of subducted slab
Length_Horiz_Slab   =  1500;    % length of overriding plate of slab
Width_Slab          =  1000;    % Width of slab (in case we run a 3D model)         

SubductionAngle     =   34;     % Subduction angle
ThermalAge_Myrs     =   50;     % Thermal age of the slab in Myrs
ThicknessCrust      =   10;



z_surface           =   0;      % initial free surface


%==========================================================================
% DEFINE SLAB
%==========================================================================

%% 1) Create the top (and bottom) of the slab

dx              =   (Length_Horiz_Slab+Length_Subduct_Slab)/100;

Slab_Top        =   ndgrid(Trench_x_location-Length_Subduct_Slab:dx:Trench_x_location+Length_Horiz_Slab); 
Slab_Top        =   unique([Slab_Top; Trench_x_location; Trench_x_location+Length_Horiz_Slab])'; % ensure that the Trench is included
Slab_Top(2,:)   =   ones(size(Slab_Top))*z_surface;                           % top of slab
ind_T           =   find(Slab_Top(1,:)==Trench_x_location);                   % trench location

Slab_Bot        =   Slab_Top;
Slab_Bot(2,:)   =   -500;

% Rotate inclined piece of slab
alpha           =   -SubductionAngle*pi/180;
R               =   [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

% Top & Bottom of slab
SlabInclined    =   [Slab_Top(:,1:ind_T), Slab_Bot(:,1:ind_T)];
SlabInclined(1,:)    =   SlabInclined(1,:)-Trench_x_location;
SlabInclined    =   R*SlabInclined;
SlabInclined(1,:)    =   SlabInclined(1,:)+Trench_x_location;

Slab_Top(:,1:ind_T) = SlabInclined(:,1:ind_T);
Slab_Bot        =   [SlabInclined(:,ind_T+1:end), [Slab_Top(1,end); SlabInclined(2,end)]];


%% 2) Compute the perpendicular distance of all "slab" points to top of slab
X2d             =   squeeze(X(1,:,:));
Z2d             =   squeeze(Z(1,:,:));

Distance        =   ones(size(X))*NaN;      % 3D matrix with distance to 
Dist_2D         =   ones(size(X2d))*NaN;    % in x-z plane

[d_min]         =   p_poly_dist(X2d(:), Z2d(:), Slab_Top(1,:), Slab_Top(2,:), false);
Dist_2D(find(Dist_2D)) = d_min;

SlabPolygon     =   [Slab_Top, Slab_Bot(:,end:-1:1)];
in              =   inpolygon(X2d,Z2d,SlabPolygon(1,:),SlabPolygon(2,:));
Dist_2D(~in)    =   NaN;

for iy=1:size(X,1);
    Distance(iy,:,:) = Dist_2D;
end


%% Set phases and temperature based on distance of top of slab
T_surface   =   20;
T_mantle    =   1350;
Phase       =   ones(size(X)); % initialize to have mantle phase
Temp        =   T_mantle*ones(size(Phase));

% Set air
ind        =    find(Z>z_surface);
Phase(ind) =    0;
Temp(ind)  =    0;

% Set Crust
ind        =    find(Distance<ThicknessCrust);
Phase(ind) =    2;

% Set 3D temperature based on halfspace cooling & distance to top of slab
kappa       =   1e-6;
ThermalAge  =   ThermalAge_Myrs*1e6*(365*24*3600);

% halfspace cooling
ind         =   ~isnan(Distance);
Temp(ind)   =   (T_mantle -T_surface) * erf(abs(Distance(ind))*1000/2/sqrt(kappa*ThermalAge)) + T_surface;

% Set Mantle Lithosphere for mantle points that have temperatures < 1200 Celcius
ind        =    find(Temp<1200 & Phase==1);
Phase(ind) =    3;


% Limit lateral size of slab (in 3D cases)
ind         =   find(Y>Width_Slab & Phase>0);
Phase(ind)  =   1;
Temp(ind)   =   T_mantle;



%==========================================================================
% PREPARE DATA FOR VISUALIZATION/OUTPUT (no need to change this)
%==========================================================================

% Prepare data for visualization/output
A = struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);

Phase    = permute(Phase,[2 1 3]);
Temp     = permute(Temp, [2 1 3]);

% Linear vectors containing coords
x        =  X(1,:,1);
y        =  Y(:,1,1);
z        =  Z(1,1,:);
X        =  permute(X,[2 1 3]);
Y        =  permute(Y,[2 1 3]);
Z        =  permute(Z,[2 1 3]);

A.W      =  W;
A.L      =  L;
A.H      =  H;
A.nump_x =  nump_x;
A.nump_y =  nump_y;
A.nump_z =  nump_z;
A.Phase  =  Phase;
A.Temp   =  Temp;
A.x      =  x(:);
A.y      =  y(:);
A.z      =  z(:);
A.Xpart  =  X;
A.Ypart  =  Y;
A.Zpart  =  Z;
A.npart_x=  npart_x;
A.npart_y=  npart_y;
A.npart_z=  npart_z;

% PARAVIEW VISUALIZATION
if (Paraview_output == 1)
    FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option
end

% SAVE PARALLEL DATA (parallel)
if (LaMEM_Parallel_output == 1)
    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition,Is64BIT);
end




