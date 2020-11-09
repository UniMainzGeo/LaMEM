% Create a 2D subduction setup with particlesm but with no temperature structure
clear


% Tell the code where the LaMEM matlab routines are, relative to the current directory 
addpath('../../matlab')

%==========================================================================
% OUTPUT OPTIONS (standard)
%==========================================================================
LaMEM_Parallel_output   =    0;

RandomNoise             =   logical(0); % add random noise to particles?

LaMEM_input_file        =   'Subduction2D_FreeSlip_Particles_Linear_DirectSolver.dat';

%% Compute 3D grid, depending on whether we are on 1 or >1 processors
if ~LaMEM_Parallel_output 
    % In the other case, we create a setup for 1 processor and defined the
    % parameters here. 
    % Important: the resolution you use here should be identical to what
    % is specified in then *.dat file!
    disp(['Creating setup for 1 processor using LaMEM file: ', LaMEM_input_file])

    % Read info from LaMEM input file & create 3D grid    
    [Grid,X,Y,Z,npart_x,npart_y,npart_z,W,L,H] =   LaMEM_ParseInputFile(LaMEM_input_file);
    
    nump_x              =   Grid.nel_x*npart_x;
    nump_y              =   Grid.nel_y*npart_y;
    nump_z              =   Grid.nel_z*npart_z;
    Parallel_partition  =   [];   % since we run this on one code
    
else
    % We perform a paralel simulation; or this a 'ProcessorPartitioning'
    % file shoule be created first by running LaMEM on the desired # of
    % processors as:
    %   mpiexec -n 2 ../../bin/opt/LaMEM -ParamFile Subduction2D_FreeSlip_MATLABParticles_Linear_DirectSolver.dat -mode save_grid
     
    % Define parallel partition file
    Parallel_partition     = 'ProcessorPartitioning_2cpu_2.1.1.bin'
    
    % Load grid from parallel partitioning file
    [Grid,X,Y,Z,npart_x,npart_y,npart_z,W,L,H] =   LaMEM_ParseInputFile(LaMEM_input_file);
    [X,Y,Z,xcoor,ycoor,zcoor,Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition,RandomNoise);
    
    % Update other variables
    nump_x  = size(X,2);
    nump_y  = size(X,1);
    nump_z  = size(X,3);
    
    % Domain parameters
    W       =   xcoor(end)-xcoor(1);    % x-dir
    L       =   ycoor(end)-ycoor(1);    % y-dir
    H       =   zcoor(end)-zcoor(1);    % z-dir
end
%==========================================================================
%%


%% 
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


%% 
%==========================================================================
% DEFINE SLAB
%==========================================================================

% 1) Create the top (and bottom) of the slab

% Create polygon of crust and mantle lithosphere
x_s             =   Trench_x_location-Length_Subduct_Slab;  % start
x_t             =   Trench_x_location;                      % trench
x_e             =   Trench_x_location+Length_Horiz_Slab;    % end

% 1a)First create a horizontal crust/slab
Slab_Crust(1,:) =   [x_s x_t x_e x_e x_t x_s];
Slab_Crust(2,:) =   [0   0     0 -ThicknessCrust -ThicknessCrust -ThicknessCrust];

Slab_ML         =   Slab_Crust;
Slab_ML(2,:)    =   [0   0     0 -ThicknessSlab -ThicknessSlab -ThicknessSlab];

% 1b) Next rotate the inclined portions
alpha           =   -SubductionAngle*pi/180;
R               =   [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

% 2) do the same for the Crust
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

Slab_ML(:,1:2)   = Inclined(:,1:2);
Slab_ML(:,end-1:end) = Inclined(:,end-1:end);
Slab_ML(2,end-2) = Slab_ML(2,end-1);                    

% Now we have a polygon that described the crust and one that describes the mantle lithosphere.
% You can plot them with:
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
T_mantle        =   1350;
Temp            =   T_mantle*ones(size(Phase));

%%
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
FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option

% SAVE DATA (parallel)
FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition);



        
