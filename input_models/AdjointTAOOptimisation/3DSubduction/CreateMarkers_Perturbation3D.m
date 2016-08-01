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
Parallel_partition     = 'ProcessorPartitioning_4cpu_2.1.2.bin'

Is64BIT                 =   logical(0);

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Domain parameters
W       =   1000;    % x-dir
L       =   200;   % y-dir
H       =   400;    % z-dir

% Number of markers in a grid cell
npart_x = 3;
npart_y = 3;
npart_z = 3;

% Element resolution
nel_y   =   16;  % 128
nel_x   =   128;  % 512 
nel_z   =   64;  % 256

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
Freeair  = 0;
Mantle   = 1;
UCL      = 2;
LCL      = 3;
Litho    = 4;
UCR      = 5;
LCR      = 6;

Phase   =   zeros(size(X)); % initialize phases

%==========================================================================
% SETUP GEOMETRY
%==========================================================================
% Add AIR
ind         =   find( Z>395 );
Phase(ind)  =   Freeair;

% Mantle
ind         =   find( Z<=260 );
Phase(ind)  =   Mantle;


% UCL
ind                 = find(Z<=395 & Z>380 & X>=0 & X<=500);
Phase(ind)  =   UCL;

% LCL
ind                 = find(Z<=380 & Z>360 & X>=0 & X<=500);
Phase(ind)  =   LCL;

% Litho
ind                 = find(Z<=360 & Z>260 & X>=0 & X<=W);
Phase(ind)  =   Litho;

% UCR
ind                 = find(Z<=395 & Z>380 & X>500 & X<=W);
Phase(ind)  =   UCR;

% LCR
ind                 = find(Z<=380 & Z>360 & X>500 & X<=W);
Phase(ind)  =   LCR;


for i = 1:nump_y/2
    X2D = X(i,:,:);
    Z2D = Z(i,:,:);
    Phase2D = Phase(i,:,:);
    
    % UCR
    Pol_x = [W 500 400 425 510 W];
    Pol_z = [395 395 200 200 380 380];
    ind = find(inpolygon(X2D,Z2D,Pol_x,Pol_z));
    Phase2D(ind)  =   UCR;
    
    % LCR
    Pol_x = [W 510 425 450 520 W];
    Pol_z = [380 380 200 200 360 360];
    ind = find(inpolygon(X2D,Z2D,Pol_x,Pol_z));
    Phase2D(ind)  =   LCR;
    
    Phase(i,:,:) = Phase2D;
    
end



%==========================================================================
% TEMPERATURE - in Celcius
%==========================================================================

Temp    =   ones(size(Z)) .*1; % initialize temperature
grad1 = 11; % K/km
grad2 = 910/100; % K/km
grad3 = 0.5; %K/km

Z_temp = abs(Z-H);

ind = find(Z >= 360);
Temp(ind) = Temp(ind) + grad1*Z_temp(ind);
Tcur = 1;

ind = find(Z >= 260 & Z<360);
Temp(ind) = Temp(ind) + grad2*Z_temp(ind) + Tcur;
Tcur = max(max(max(Temp)));

ind = find(Z < 260);
Temp(ind) = Temp(ind) + grad3*Z_temp(ind) + Tcur;


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




