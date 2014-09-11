% Create input files (parallel or sequential) for markers with Phase and Temperature distributions
% Non-dimensionalization of the setup is done internally in LaMEM!

clear
addpath ../../../matlab

% See model setup in Paraview 1-YES; 0-NO
Paraview_output       = 1;          

% Output parallel files for LaMEM, using a processor distribution file (Setup.Extern = 1)
LaMEM_Parallel_output = 0;          
Parallel_partition    = 'ProcessorPartitioning_1cpu_1.1.1.bin';

% Output a single file containing particles information for LaMEM (Setup.Extern = 2)
LaMEM_OLD_WAY_output  = 1;          

% Domain parameters
W       =   300e3;
L       =   100e3;
H       =   200e3;

% Number of markers
nump_x  =   2*90;  
nump_y  =   2*60;
nump_z  =   60;

% Number of markers in a grid cell
npart_x = 3;
npart_y = 3;
npart_z = 3;

% Model specific parameters
dx_reg  =   W/(nump_x-1);
dy_reg  =   L/(nump_y-1);
dz_reg  =   H/(nump_z-1);
x_left  =   -150e3;
y_front =   0;
z_bot   =   0;

% Create grid and initialize phases and temperature distribution
[X,Y,Z] =   meshgrid([x_left:dx_reg:x_left+W], [y_front:dy_reg:y_front+L], [z_bot:dz_reg:z_bot+H]);

Phase   =   zeros(size(X));     %   Contains phases
Temp    =   zeros(size(X));     %   Contains temperatures

% Create slab-like initial distribution====================================
%Slab
margin              =   10e3;          % attached margin=0 or unattached to the boundaries    
H_slab              =   H/9;
Depth_slab          =   0.75*H;
ThicknessAir        =   0.05*H;

%==========================================================================
%PHASES
mantle   = 0;
air      = 1;
slab     = 2;

%SLAB
ind         =   find( X>(x_left+margin) & X<0 & Z>(H-H_slab) & Y>margin &Y<(L-margin) ); 
Phase(ind)  =   slab;

ind         =   find( (Z<(H-(X-0/2))) &  (Z>((H-H_slab)-(X-0/2))) & Z>Depth_slab & Y>margin &Y<(L-margin)  ); 
Phase(ind)  =   slab;

% Add AIR
ind         =   find( Z>(H-ThicknessAir) );
Phase(ind)  =   air;

%==========================================================================
% TEMPERATURE - in Celcius
%Temp = (H-Z)./H*0 + 0.5 + (rand(size(Z))-0.5)*0.05 + 0*1000;

%==========================================================================
% Prepare data for visualization/output
A = struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);
Phase       = permute(Phase,[2 1 3]);
Temp        = permute(Temp, [2 1 3]);

% Linear vectors containing coords
x = X(1,:,1)-x_left;
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

% SAVE DATA IN 1 FILE (sequential)
if (LaMEM_OLD_WAY_output == 1)
    PhaseOrig   = Phase;
    PhaseVec(1) = nump_z;
    PhaseVec(2) = nump_y;
    PhaseVec(3) = nump_x;
    PhaseVec    = [PhaseVec(:); Phase(:); Temp(:)];
    
    % Save data to file
    ParticleOutput  =   'MarkersInput3D.dat';
    PassiveOutput   =   'PassiveInput.dat';
    
    PetscBinaryWrite(ParticleOutput, PhaseVec);
    
end

clearvars -except A Paraview_output LaMEM_Parallel_output Parallel_partition

% PARAVIEW VISUALIZATION
if (Paraview_output == 1)
    
    WriteLaMEMSetup_Matlab2VTK(A,'BINARY'); % default option
    %WriteLaMEMSetup_Matlab2VTK(A,'ASCII'); % for debugging only (slow) 
end

% SAVE PARALLEL DATA (parallel)
if (LaMEM_Parallel_output == 1)
    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition);
end




