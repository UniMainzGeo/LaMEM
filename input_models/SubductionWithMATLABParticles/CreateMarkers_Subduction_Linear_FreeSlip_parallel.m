% Create a 2D subduction setup with particles & temperature, 


clear


% Tell the code where the LaMEM matlab routines are 
addpath('../../matlab')


%==========================================================================
% OUTPUT OPTIONS (standard)
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        =    1;

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output  =    0;

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh               =    1;

% random noise of particles
RandomNoise             =   logical(0);

Is64BIT                 =   logical(0); % only used for some 

if LaMEM_Parallel_output
    % We perform a paralell simulation; or this a 'ProcessorPartitioning'
    % file shoule be created first by running LaMEM on the desired # of
    % processors as:
    %   mpiexec -n 2 ../../bin/opt/LaMEM -ParamFile Subduction2D_FreeSlip_MATLABParticles_Linear_DirectSolver.dat -mode save_grid
     
    % Define parallel partition file
    Parallel_partition     = 'ProcessorPartitioning_2cpu_2.1.1.bin'
    
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
    

else
    
    % In the other case, we create a setup for 1 processor and defined the
    % parameters here. 
    % Important: the resolution you use here should be identical to what
    % is specified in then *.dat file!
    disp('Creating setup for 1 processor')
    
    % Domain parameters
    W       =   3000;       % x-dir
    L       =   10;         % y-dir
    H       =   660;        % z-dir
    
    x_left  =   1500;      % coord of the left margin
    y_front =   0;      % coord of the front margin
    z_bot   =   -660;      % coord of the bottom margin

    % Element resolution
    nel_x   =   256;
    nel_y   =   2;
    nel_z   =   128;
    
    % Number of markers in a grid cell
    npart_x = 3;
    npart_y = 3;
    npart_z = 3;
    
    %-----------------------------------------------------------------
    % not to be changed
    % Number of markers
    nump_x  =   nel_x*npart_x;
    nump_y  =   nel_y*npart_y;
    nump_z  =   nel_z*npart_z;
    
    dx      =   W/nump_x;
    dy      =   L/nump_y;
    dz      =   H/nump_z;
    x       =   [x_left  + dx*0.5 : dx : x_left+W  - dx*0.5 ];
    y       =   [y_front + dy*0.5 : dy : y_front+L - dy*0.5 ];
    z       =   [z_bot   + dz*0.5 : dz : z_bot+H   - dz*0.5 ];
    [X,Y,Z] =   meshgrid(x,y,z);
    %
    %----------------------------------------------------------------------
    
end

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




