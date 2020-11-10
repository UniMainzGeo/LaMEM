% --------------------------------------------------------%
%%%%% 3D shear localization model with random noise%%%%%
% --------------------------------------------------------%

% This script creates LaMEM input files (parallel and/or sequential) for markers 
% Files contain: marker coordinates, phase and temperature distributions
% WARNING: The model setup should be dimensional! Non-dimensionalization is done internally in LaMEM!

clear
% add matlab file to create phases
% addpath /data/ljeff/software/LaMEM/matlab % sith
% addpath /local/home/ljeff/software/LaMEM/matlab%gaia
addpath('../../matlab')

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
RandomNoise             = logical(1);

% Parallel partition file
%Parallel_partition     = 'ProcessorPartitioning_2048cpu_16.16.8.bin';
Parallel_partition     = 'ProcessorPartitioning_2cpu_2.1.1.bin'

Is64BIT                 =   logical(0);

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Domain parameters
W       =   100;    % x-dir
L       =   0.7;    % y-dir
H       =   5;      % z-dir

% Number of markers in a grid cell
npart_x = 3;
npart_y = 3;
npart_z = 3;

% Element resolution
nel_y   =   2;
nel_x   =   256;
nel_z   =   128;

% Number of markers
nump_x  =   nel_x*npart_x;  
nump_y  =   nel_y*npart_y;
nump_z  =   nel_z*npart_z;


% Model specific parameters
dx  =   W/nump_x;
dy  =   L/nump_y;
dz  =   H/nump_z;
x_left  =   -W/2;      % coord of the left margin
y_front =   -L/2;      % coord of the front margin
z_bot   =   -25;      % coord of the bottom margin

%Slab - related parameters
ThickAir            =   5;
ThickUpperCrust1    =   5;
ThickUpperCrust2    =   5;
ThickUpperCrust3    =   5;
ThickUpperCrust4    =   5;

z_air               =   H-ThickAir;
z_UC1               =   z_air-ThickUpperCrust1;
z_UC2               =   z_UC1-ThickUpperCrust2;
z_UC3               =   z_UC2-ThickUpperCrust3;
z_UC4               =   z_UC3-ThickUpperCrust4;


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
Air             =   0;
UC1             =   1;
UC2             =   2;
UC3             =   3;
UC4             =   4;

Phase   =   zeros(size(X)); % initialize phases

%==========================================================================
% SETUP GEOMETRY
%==========================================================================
ind         =   find(Z>=z_UC1 & Z<z_air);
Phase(ind)  =   1;
ind         =   find(Z>=z_UC2 & Z<z_UC1);
Phase(ind)  =   2;
ind         =   find(Z>=z_UC3 & Z<z_UC2);
Phase(ind)  =   3;
% ind         =   find(Z>=z_UC4 & Z<z_UC3);
ind         =   find( Z<z_UC3);
Phase(ind)  =   4;

% add notch in center
ind         =   find( Z<(z_UC3+0.1) & abs(X)<1 & abs(Y)<1);
Phase(ind)  =   4;

% Add pattern at the base 
Wpattern = 2;

x1 = x_left;
for i=1:100
    x2 = x1+Wpattern;
    ind = find(X>x1 & X<x2 & Phase==4);
    Phase(ind) = 5;
    
    x1 = x2 + 2*Wpattern;
    
    if x1>(x_left+W)
        break
    end
    
end

    
% % Add low density anomaly at base (graite intrusion)
% ind         =   find( Z<(z_bot+10) & abs(X)<10 & abs(Y)<10);
% Phase(ind)  =   6;
%     



%==========================================================================
% TEMPERATURE - in Celcius
%==========================================================================
Tbottom     =   275;
Temp        =   Tbottom.*(Z/z_bot);     % linear temp gradient between top & bottom
ind         =   find(Z>z_air);          % air has zero
Temp(ind)   =   0;
                

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




