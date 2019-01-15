% Create a 2D subduction setup with particles & temperature, 


clear


% Tell the code where the LaMEM matlab routines are 
addpath ../../matlab


% Define parallel partition file
Parallel_partition     = 'ProcessorPartitioning_2cpu_2.1.1.bin'


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
Length_Subduct_Slab =  200;     % length of subducted slab
Length_Horiz_Slab   =  1500;    % length of overriding plate of slab
Width_Slab          =  750;     % Width of slab (in case we run a 3D model)         


SubductionAngle     =   34;     % Subduction angle
ThicknessCrust      =   10;
ThicknessSlab       =   75;    % Thickness of mantle lithosphere


%==========================================================================
% DEFINE SLAB
%==========================================================================

%% 1) Create the top (and bottom) of the slab

% Create polygon of crust and mantle lithoshere
x_s             =   Trench_x_location-Length_Subduct_Slab;  % start
x_t             =   Trench_x_location;                      % trench
x_e             =   Trench_x_location+Length_Horiz_Slab;    % end

% First create a horizontal crust/slab
Slab_Crust(1,:) =   [x_s x_t x_e x_e x_t x_s];
Slab_Crust(2,:) =   [0   0     0 -ThicknessCrust -ThicknessCrust -ThicknessCrust];

Slab_ML         =   Slab_Crust;
Slab_ML(2,:)    =   [0   0     0 -ThicknessSlab -ThicknessSlab -ThicknessSlab];

% Next rotate the inclined portions
alpha           =   -SubductionAngle*pi/180;
R               =   [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

% Crust
Inclined        =   Slab_Crust(:,[1:2 end-1:end]);
Inclined(1,:)   =   Inclined(1,:)-Trench_x_location;
Inclined        =   R*Inclined;
Inclined(1,:)   =   Inclined(1,:)+Trench_x_location;
Inclined(2,end-1) =   -ThicknessCrust;

Slab_Crust(:,1:2) = Inclined(:,1:2);
Slab_Crust(:,end-1:end) = Inclined(:,end-1:end);

% rotate slab
Inclined        =   Slab_ML(:,[1:2 end-1:end]);
Inclined(1,:)   =   Inclined(1,:)-Trench_x_location;
Inclined        =   R*Inclined;
Inclined(1,:)   =   Inclined(1,:)+Trench_x_location;

Slab_ML(:,1:2) = Inclined(:,1:2);
Slab_ML(:,end-1:end) = Inclined(:,end-1:end);
Slab_ML(2,end-2) = Slab_ML(2,end-1);                    % NOTE: the final thickness is no longer the specified slab thicknes!!

% Now we have a polygon that described the crust and one that describes the
% mantle lithosphere.
% You can plot them with
% fill(Slab_ML(1,:),Slab_ML(2,:),'bo-',Slab_Crust(1,:),Slab_Crust(2,:),'ro-'), axis equal

%% Set phases based on whether we are in the crust or in the mantle lithosphere
Phase           =   zeros(size(X)); % initialize to have mantle phase

% Note that the mantle lithosphere should be set first in the way we define
% the polygons
in              =   inpolygon(X,Z,Slab_ML(1,:),Slab_ML(2,:));       % ML
Phase(in)       =   2;

in              =   inpolygon(X,Z,Slab_Crust(1,:),Slab_Crust(2,:)); % crust
Phase(in)       =   1;

% For 3D setups, take the width of the slab in Y-direction into account
ind             =  find(Y>Width_Slab);
Phase(ind)      =   0;

% Temperature is not used in the setup, so set it to a constant value
T_mantle    =   1350;
Temp        =   T_mantle*ones(size(Phase));



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




