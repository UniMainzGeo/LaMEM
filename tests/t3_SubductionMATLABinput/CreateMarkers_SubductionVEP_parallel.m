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
Paraview_output        = 0;

% Output a single file containing particles information for LaMEM (msetup = redundant)
LaMEM_Redundant_output = 0;%0;

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output  = 1;

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh               = 1;

% random noise of particles
RandomNoise             = logical(1);

% Parallel partition file

Parallel_partition     = '../ProcessorPartitioning_8cpu_2.1.4.bin';

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Domain parameters
W       =   2000; % x-dir
L       =   1;    % y-dir
H       =   680;  % z-dir

% Number of markers in a grid cell
npart_x = 3;
npart_y = 3;
npart_z = 3;
% Element resolution
nel_x   =   256;
nel_y   =   2;
nel_z   =   128;
% Number of markers
nump_x  =   nel_x*npart_x;
nump_y  =   nel_y*npart_y;
nump_z  =   nel_z*npart_z;

ThickAir=   20;
% Model specific parameters
dx  =   W/nump_x;
dy  =   L/nump_y;
dz  =   H/nump_z;
x_left  =   -W/2;      % coord of the left margin
y_front =   0;      % coord of the front margin
z_bot   =   ThickAir-H;      % coord of the bottom margin

ThickOC     =   8;  % thickness of oceanic crust [km]
ThickSP     =   80; % thickness of subducting plate [km]
ThickOP     =   80; % thickness of overriding plate [km]
ThickSML    =   ThickSP - ThickOC; % thickness of subducting mantle lithosphere
ThickOML    =   ThickOP - ThickOC; % thickness of overriding mantle lithosphere

z_air       =   0;
z_oc        =   z_air - ThickOC;
z_sp        =   z_air - ThickSP;
z_op        =   z_air - ThickOP;

xtrench     =   0; % initial trench position
theta       =   30 * pi/180; % initial subducting angle [rad]





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
Mantle          =   1;
WeakZone        =   2;
OceanicCrust    =   3;
SubductPlate    =   4;
OverridingPlate =   5;



Phase   =   ones(size(X)); % initialize phases


% Weak zone geometry
Hweak = ThickSP + 20;
x1 = xtrench; z1 = z_oc;
x2 = x1-50;   z2 = z_oc;
x3 = x2-Hweak/tan(theta); z3 = z_oc - Hweak;






%==========================================================================
% TEMPERATURE - in Celcius
%==========================================================================
Ttop        =   0;
Tmantle     =   1280; % mantle potential temperature
dTdz        =   0.3;  % adiabatic gradient [oC/km]
Tbottom     =   Tmantle + (H-ThickAir)*dTdz; % bottom temperature

SecYear     =   3600*24*365.25; % 1 year in second
kappa       =   1e-6; % thermal diffusivity
Tage_SP     =   80e6*SecYear; % Theraml age of subducting plate
Tage_OP     =   40e6*SecYear; % Thermal age of overriding plate
% thermal structure in the mantle
Temp        =   Tmantle + abs(Z-z_air)*dTdz;
% thermal structure in subducting plate
TSP         =   Tmantle + abs(z_sp-z_air)*dTdz;
dT          =   (TSP - Ttop)*(1 - erf(abs(z_sp-z_air)*1000/2/sqrt(Tage_SP*kappa)) );
ind         =   find(Z>=z_sp & Z<=z_air);
Temp(ind)   =   Ttop + (TSP+dT-Ttop)*erf(abs(Z(ind)-z_air).*1000/2/sqrt(Tage_SP*kappa));
% thermal structure in overriding plate
TOP         =   Tmantle + abs(z_op-z_air)*dTdz;
dT          =   (TOP - Ttop)*(1 - erf(abs(z_op-z_air)*1000/2/sqrt(Tage_OP*kappa)) );
ind         =   find(Z>=z_op & Z<=z_air & X<=x2+(Z-z2)./tan(theta) );
Temp(ind)   =   Ttop + (TOP+dT-Ttop)*erf(abs(Z(ind)-z_air).*1000/2/sqrt(Tage_OP*kappa));
% constrain air to Ttop
ind         =   find(Z>=z_air);
Temp(ind)   =   Ttop;



%==========================================================================
% SETUP GEOMETRY
%==========================================================================
% Air
ind         =   find(Z>z_air);
Phase(ind)  =   Air;
% Oceanic crust
ind         =   find(Z>=z_oc & Z<=z_air);
Phase(ind)  =   OceanicCrust;
% Subducting plate
ind         =   find(Z>=z_sp & Z<=z_oc & X>=x1+(Z-z1).*(x3-x1)/(z3-z1) & Temp<=1200);
Phase(ind)  =   SubductPlate;
% Overriding plate
ind         =   find(Z>=z_op & Z<=z_oc & X<=x2+(Z-z2)./tan(theta) & Temp<=1200);
Phase(ind)  =   OverridingPlate;
% weak zone
ind         =   find(Z<=z_oc & X>=x2+(Z-z2)./tan(theta) & X<=x1+(Z-z1).*(x3-x1)/(z3-z1));
Phase(ind)  =   WeakZone;


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

	!rm -rf ../markers_p8
	!mv -f markers ../markers_p8

end




