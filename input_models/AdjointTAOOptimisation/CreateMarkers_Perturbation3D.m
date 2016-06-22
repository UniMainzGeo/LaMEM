% --------------------------------------------------------%
%%%%% 3D shear localization model with random noise%%%%%
% --------------------------------------------------------%

% This script creates LaMEM input files (parallel and/or sequential) for markers 
% Files contain: marker coordinates, phase and temperature distributions
% WARNING: The model setup should be dimensional! Non-dimensionalization is done internally in LaMEM!

clear
% add matlab file to create phases
addpath /home/greuber/Codes/lamem/matlab

%==========================================================================
% OUTPUT OPTIONS
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;          

% Output a single file containing particles information for LaMEM (msetup = redundant)
LaMEM_Redundant_output = 0; %0;        

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output  = 1;     

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh               = 1;  

% random noise of particles
RandomNoise             = logical(0);

% Parallel partition file
%Parallel_partition     = 'ProcessorPartitioning_2048cpu_16.16.8.bin';
Parallel_partition     = 'ProcessorPartitioning_2cpu_1.1.2.bin'

Is64BIT                 =   logical(0);

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Domain parameters
W       =   1.5;    % x-dir
L       =   0.001;   % y-dir
H       =   1;    % z-dir

% Number of markers in a grid cell
npart_x = 5;
npart_y = 5;
npart_z = 5;

% Element resolution
nel_y   =   2;
nel_x   =   6;
nel_z   =   20;

% Number of markers
nump_x  =   nel_x*npart_x;  
nump_y  =   nel_y*npart_y;
nump_z  =   nel_z*npart_z;


% Model specific parameters
dx  =   W/nump_x;
dy  =   L/nump_y;
dz  =   H/nump_z;
x_left  =   0;      % coord of the left margin
y_front =   0;    % coord of the front margin
z_bot   =   0;      % coord of the bottom margin


% Geometry related parameters
Hi			=	0.5;		% height of air
A0          =   0.01;



%==========================================================================
% MESH GRID
%==========================================================================

% Create new uniform grid
if LoadMesh == 0
    x = [x_left  + dx*0.5 : dx : x_left+W  - dx*0.5 ];
    y = [y_front + dy*0.5 : dy : y_front+L - dy*0.5 ];
    z = [z_bot   + dz*0.5 : dz : z_bot+H   - dz*0.5 ];
    [X,Y,Z] =   meshgrid(x,y,z);
end

% Load grid from parallel partitioning file 
if LoadMesh == 1    
    [~,~,~,xcoor,ycoor,zcoor,Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition,RandomNoise,Is64BIT);
    X = Xpart;
    Y = Ypart;
    Z = Zpart;
    % Update other variables
    nump_x = size(X,2);
    nump_y = size(X,1);
    nump_z = size(X,3);
end

%==========================================================================
% PHASES
%==========================================================================
Freeair   		= 0;
Overburden      = 1;
Salt     		= 2;

Phase   =   zeros(size(X)); % initialize phases

%==========================================================================
% SETUP GEOMETRY
%==========================================================================
% Add AIR
ind         =   find( Z>(1) );
Phase(ind)  =   Freeair;

ind         =   find( Z>=(0.5) & Z<=(1));
Phase(ind)  =   Overburden;

cd = linspace(0,pi/2,100);
ce = sin(cd)./(1/A0);  % to have A0 as maximum

ud = linspace(0,1.5,100);


x_Pol = [0     0      ud        1.5   1.5];
z_Pol = [0     0.5    0.5+ce    0.5    0];

ind = inpolygon(X(:),Z(:),x_Pol,z_Pol);

Phase(ind)  =   Salt;

% ind         =   find( Z>=(0) & Z<(0.5));
% Phase(ind)  =   Salt;

% %==========================================================================
% % ADJIOINT - Find indices
% %==========================================================================
% length_vel = (nel_y*nel_x*nel_z)*3;
% length_pres = (nel_y*nel_x*nel_z);
%  
% length_tot = length_vel + length_pres;
% 
% dx  =   W/nel_x;
% dy  =   L/nel_y;
% dz  =   H/nel_z;
% x = [x_left  + dx*0.5 : dx : x_left+W  - dx*0.5 ];
% y = [y_front + dy*0.5 : dy : y_front+L - dy*0.5 ];
% z = [z_bot   + dz*0.5 : dz : z_bot+H   - dz*0.5 ];
% [X_mesh,Y_mesh,Z_mesh] =   meshgrid(x,y,z);
% 
% [val,ind] = max(z(z<=0.5+A0));
% 
% ind = find(Z_mesh==val);
% 
% % since I'm looking for a y velocity it should be the 2. entry in the
% % velocity entries -1 because they start with 0 in C
% ind_adj = 2+(4*ind(15))-1
% 

%==========================================================================
% TEMPERATURE - in Celcius
%==========================================================================

Temp    =   ones(size(Z));% .*273.15; % initialize temperature

%==========================================================================
% PREPARE DATA FOR VISUALIZATION/OUTPUT
%==========================================================================

% Prepare data for visualization/output
A = struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);

Phase       = permute(Phase,[2 1 3]);
Temp        = permute(Temp, [2 1 3]);

% Linear vectors containing coords
x = X(1,:,1);
y = Y(:,1,1);
z = Z(1,1,:);
X        = permute(X,[2 1 3]);
Y        = permute(Y,[2 1 3]);
Z        = permute(Z,[2 1 3]);

A.W      = W;
A.L      = L;
A.H      = H;
A.nump_x = nump_x;
A.nump_y = nump_y;
A.nump_z = nump_z;
A.Phase  = Phase;
A.Temp   = Temp;
A.x      = x(:);
A.y      = y(:);
A.z      = z(:);
A.Xpart  = X;
A.Ypart  = Y;
A.Zpart  = Z;
A.npart_x= npart_x;
A.npart_y= npart_y;
A.npart_z= npart_z;



% SAVE DATA IN 1 FILE (redundant)
if (LaMEM_Redundant_output == 1)
    PhaseVec(1) = nump_z;
    PhaseVec(2) = nump_y;
    PhaseVec(3) = nump_x;
    PhaseVec    = [PhaseVec(:); X(:); Y(:); Z(:); Phase(:); Temp(:)];
    
    % Save data to file
    ParticleOutput  =   'MarkersInput3D.dat';
    
    PetscBinaryWrite(ParticleOutput, PhaseVec);
    
end

% Clearing up some memory for parallel partitioning
clearvars -except A Paraview_output LaMEM_Parallel_output Parallel_partition Is64BIT

% PARAVIEW VISUALIZATION
if (Paraview_output == 1)
    
    FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option
    %FDSTAGWriteMatlab2VTK(A,'ASCII'); % for debugging only (slow) 
end

% SAVE PARALLEL DATA (parallel)
if (LaMEM_Parallel_output == 1)
    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition,Is64BIT);
end




