% ----------------------------%
%%%%% 3D Subduction Model %%%%%
% ----------------------------%

% This script creates LaMEM input files (parallel and/or sequential) for markers 
% Files contain: marker coordinates, phase and temperature distributions
% WARNING: The model setup should be dimensional! Non-dimensionalization is done internally in LaMEM!

clear
addpath ../../../matlab

%==========================================================================
% OUTPUT OPTIONS
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;          

% Output a single file containing particles information for LaMEM (msetup = redundant)
LaMEM_Redundant_output = 0;        

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output  = 0;     

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh               = 0;  

% Parallel partition file
Parallel_partition     = 'ProcessorPartitioning_1cpu_1.1.1.bin';

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Domain parameters
W       =   300e3;  % x-dir
L       =   100e3;  % y-dir
H       =   200e3;  % z-dir

% Number of markers
nump_x  =   180;  
nump_y  =   120;
nump_z  =   60;

% Number of markers in a grid cell
npart_x = 3;
npart_y = 3;
npart_z = 3;

% Model specific parameters
dx  =   W/nump_x;
dy  =   L/nump_y;
dz  =   H/nump_z;
x_left  =   -150e3; % coord of the left margin
y_front =   0;      % coord of the front margin
z_bot   =   0;      % coord of the bottom margin

%Slab - related parameters
margin              =   10e3; % attached margin=0 or unattached to the boundaries    
H_slab              =   H/9;
Depth_slab          =   0.75*H;
ThicknessAir        =   0.05*H;

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
    [X,Y,Z,x,y,z] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition);
    
    % Update other variables
    nump_x = size(X,2);
    nump_y = size(X,1);
    nump_z = size(X,3);
end

%==========================================================================
% PHASES
%==========================================================================
mantle   = 0;
air      = 1;
slab     = 2;

Phase   =   zeros(size(X)); % initialize phases

%==========================================================================
% SETUP GEOMETRY
%==========================================================================
% SLAB
ind         =   find( X>(x_left+margin) & X<0 & Z>(H-H_slab) & Y>margin &Y<(L-margin) ); 
Phase(ind)  =   slab;

ind         =   find( (Z<(H-(X-0/2))) &  (Z>((H-H_slab)-(X-0/2))) & Z>Depth_slab & Y>margin &Y<(L-margin)  ); 
Phase(ind)  =   slab;

% Add AIR
ind         =   find( Z>(H-ThicknessAir) );
Phase(ind)  =   air;

%==========================================================================
% TEMPERATURE - in Celcius
%==========================================================================
%Temp = (H-Z)./H*0 + 0.5 + (rand(size(Z))-0.5)*0.05 + 0*1000;
Temp    =   zeros(size(Z)); % initialize temperature                 

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
A.npart_x= npart_x;
A.npart_y= npart_y;
A.npart_z= npart_z;

X        = permute(X,[2 1 3]);
Y        = permute(Y,[2 1 3]);
Z        = permute(Z,[2 1 3]);

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
clearvars -except A Paraview_output LaMEM_Parallel_output Parallel_partition

% PARAVIEW VISUALIZATION
if (Paraview_output == 1)
    
    FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option
    %FDSTAGWriteMatlab2VTK(A,'ASCII'); % for debugging only (slow) 
end

% SAVE PARALLEL DATA (parallel)
if (LaMEM_Parallel_output == 1)
    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition);
end




