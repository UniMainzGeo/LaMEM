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
Paraview_output        = 1 ;

% Output a single file containing particles information for LaMEM (msetup = redundant)
LaMEM_Redundant_output = 0;%0;

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output  = 1;

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh               = 1;

% random noise of particles
RandomNoise             = logical(0);

% Parallel partition file

Parallel_partition     = 'ProcessorPartitioning_1cpu_1.1.1.bin';

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Domain parameters
W       =   2                       ; % x-dir
L       =   2                        ; % y-dir
H       =   100                    ; % z-dir

% Number of markers in a grid cell
npart_x =   3                      ;
npart_y =   3                      ;
npart_z =   3                      ;
% Element resolution
nel_x   =  2                      ;
nel_y   =  2                      ;
nel_z   =  64                    ; 
% Number of markers
nump_x  =   nel_x*npart_x          ;
nump_y  =   nel_y*npart_y          ;
nump_z  =   nel_z*npart_z          ;

% Model specific parameters
dx      =   W/nump_x             ;
dy      =   L/nump_y              ;
dz      =   H/nump_z              ;
% Temperature distribution input data
Temp_max   =   400 ;    % Maximum temperature perturbation[Celsius Degree]
sigma=(H/2*1e3);			% Domain half-width [m]



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
%=========================================================================
Phase   =   zeros(size(X));
Temp    =  ones(size(X))*200;

%==========================================================================
% SETUP Temperature
%==========================================================================
Temp=Temp_max.*exp(-(Z.*1e3).^2./sigma^2);
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

end

% SAVE PARALLEL DATA (parallel)
if (LaMEM_Parallel_output == 1)

    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition,logical(0));

	!rm -rf ../markers_pT1
	!mv -f markers ../markers_pT1

end




