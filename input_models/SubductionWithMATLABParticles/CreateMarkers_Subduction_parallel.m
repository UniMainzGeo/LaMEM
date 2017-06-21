% --------------------------------------------------------%
%%%%% 3D shear localization model with random noise%%%%%
% --------------------------------------------------------%

% This script creates LaMEM input files (parallel and/or sequential) for markers 
% Files contain: marker coordinates, phase and temperature distributions
% WARNING: The model setup should be dimensional! Non-dimensionalization is done internally in LaMEM!

clear

addpath ../../matlab


%==========================================================================
% OUTPUT OPTIONS
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;          

% Output a single file containing particles information for LaMEM (msetup = redundant)
LaMEM_Redundant_output = 1;%0;        

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output  = 0;     

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh               = 1;  

% random noise of particles
RandomNoise             = logical(1);

% Parallel partition file
Parallel_partition     = 'ProcessorPartitioning_1cpu_1.1.1.bin';

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Domain parameters
W           =   500;  % x-dir
L           =   1;    % y-dir
H           =   130;  % z-dir

% Number of markers in a grid cell
npart_x     =   3;
npart_y     =   3;
npart_z     =   3;
% Element resolution
nel_x       =   64;
nel_y       =   2;
nel_z       =   32;
% Number of markers
nump_x      =   nel_x*npart_x;  
nump_y      =   nel_y*npart_y;
nump_z      =   nel_z*npart_z;

ThickAir    =   10;
% Model specific parameters
dx          =   W/nump_x;
dy          =   L/nump_y;
dz          =   H/nump_z;
x_left      =   -W/2;      % coord of the left margin
y_front     =   0;      % coord of the front margin
z_bot       =   ThickAir-H;      % coord of the bottom margin

ThickSlab   =   40;
XSlab       =   -50;
z_air       =   0;
z_slab      =   z_air - ThickSlab;





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
    [~,~,~,xcoor,ycoor,zcoor,Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition,RandomNoise,logical(0));
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
Slab            =   1;
Mantle          =   2;


Phase   =   2*ones(size(X)); % initialize phases

%==========================================================================
% SETUP GEOMETRY
%==========================================================================
% Air
ind         =   find(Z>z_air);
Phase(ind)  =   Air;
% slab
ind         =   find(Z>=z_slab & Z<=z_air & X>=XSlab);
Phase(ind)  =   Slab;



%==========================================================================
% TEMPERATURE - in Celcius
%==========================================================================
Ttop        =   0;
Tbottom     =   1300;
Temp        =   1200 + (Tbottom - 1200)/(z_bot-z_air).*(Z-z_air);
ind         =   find(Z>=z_air);
Temp(ind)   =   Ttop;
ind         =   find(Phase==Slab);
kappa       =   1e-6;
Tage        =   50e6*3600*24*365.25;
Temp(ind)   =   Ttop + (Tbottom-Ttop).*erf(abs(Z(ind)-z_air)*1000./2/sqrt(kappa*Tage));

                

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
clearvars -except A Paraview_output LaMEM_Parallel_output Parallel_partition

% PARAVIEW VISUALIZATION
if (Paraview_output == 1)
    
    FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option
    %FDSTAGWriteMatlab2VTK(A,'ASCII'); % for debugging only (slow) 
end

% SAVE PARALLEL DATA (parallel)
if (LaMEM_Parallel_output == 1)
    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition,logical(0));
end




