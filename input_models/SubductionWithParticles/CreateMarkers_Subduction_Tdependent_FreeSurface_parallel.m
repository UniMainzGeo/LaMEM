% Create a 2D subduction setup with particles & temperature
clear

% Tell the code where the LaMEM matlab routines are 
addpath ../../matlab

LaMEM_Parallel_output   =    0;

RandomNoise             =   logical(0); % add random noise to particles?

LaMEM_input_file        =   'Subduction2D_FreeSurface_Particles_Nonlinear_DirectSolver.dat';

%% Compute 3D grid, depending on whether we are on 1 or >1 processors
if ~LaMEM_Parallel_output 
    % In the other case, we create a setup for 1 processor and defined the
    % parameters here. 
    % Important: the resolution you use here should be identical to what
    % is specified in then *.dat file!
    disp(['Creating setup for 1 processor using LaMEM file: ', LaMEM_input_file])

    % Read info from LaMEM input file & create 3D grid    
    [npart_x,npart_y,npart_z,Grid,X,Y,Z,W,L,H] =   LaMEM_ParseInputFile(LaMEM_input_file);
    
    nump_x              =   Grid.nel_x*npart_x;
    nump_y              =   Grid.nel_y*npart_y;
    nump_z              =   Grid.nel_z*npart_z;
    Parallel_partition  =   [];   % since we run this on one code
    
else
    % We perform a paralel simulation; or this a 'ProcessorPartitioning'
    % file shoule be created first by running LaMEM on the desired # of
    % processors as:
    %   mpiexec -n 2 ../../bin/opt/LaMEM -ParamFile Subduction2D_FreeSlip_MATLABParticles_Linear_DirectSolver.dat -mode save_grid
    disp(['Creating setup in parallel using LaMEM file: ', LaMEM_input_file])
      
    % Define parallel partition file
    Parallel_partition                          =   'ProcessorPartitioning_4cpu_4.1.1.bin'
    
    % Load grid from parallel partitioning file
    [npart_x,npart_y,npart_z]  =   LaMEM_ParseInputFile(LaMEM_input_file);
    [X,Y,Z,xcoor,ycoor,zcoor,Xpart,Ypart,Zpart] =   FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition,RandomNoise);
    
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
%==========================================================================
% SPECIFY PARAMETERS OF THE SLAB
%==========================================================================
Trench_x_location   = -500;     % trench location
Length_Subduct_Slab =  200;     % length of subducted slab
Length_Horiz_Slab   =  1500;    % length of overriding plate of slab
Width_Slab          =  1000;    % Width of slab (in case we run a 3D model)         

SubductionAngle     =   20;     % Subduction angle
ThermalAge_Myrs     =   50;     % Thermal age of the slab in Myrs
ThicknessCrust      =   10;

ThicknessSlab       =   400;    % Thickness of slab box; 

z_surface           =   0;      % initial free surface

T_mantle            =   1350;
T_surf              =   20;     


%==========================================================================
% DEFINE SLAB
%==========================================================================
Phase               =   ones(size(X));                 % initialize to have mantle phase
Temp                =   T_mantle*ones(size(Phase));     

% Add horizontal part of slab 
BoxSides            =   [Trench_x_location (Trench_x_location+Length_Horiz_Slab) min(Y(:)) Width_Slab -ThicknessSlab 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 1,...
                                    'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);                         % Set slab to mantle lithosphere phase

BoxSides            =   [Trench_x_location (Trench_x_location+Length_Horiz_Slab) min(Y(:)) Width_Slab -ThicknessCrust 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 2);               % Add crust (will override the mantle lithosphere phase above)

% Add inclined part of slab
RotPt               =   [Trench_x_location, 0, 0];

BoxSides            =   [(Trench_x_location-Length_Subduct_Slab) Trench_x_location min(Y(:)) Width_Slab -ThicknessSlab 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 1, 'RotationPoint',RotPt, 'DipAngle', -SubductionAngle, ...
                                    'TempType','Halfspace', 'topTemp',T_surf,'botTemp', T_mantle);                          % mantle lithosphere

BoxSides            =   [(Trench_x_location-Length_Subduct_Slab) Trench_x_location min(Y(:)) Width_Slab -ThicknessCrust 0];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 2, 'RotationPoint',RotPt, 'DipAngle', -SubductionAngle);       	% crust (will override the slab phase above

% Set Mantle Lithosphere for mantle points that have temperatures < 1200 Celcius
ind                 =    find(Temp<1200 & Phase==1);
Phase(ind)          =    3;

% Add sticky air
BoxSides            =   [min(X(:)) max(X(:))  min(Y(:)) max(Y(:)) 0 max(Z(:))];  % [Left Right Front Back Bottom Top] of the box
[Phase,Temp]        =   AddBox(Phase,Temp,X,Y,Z,BoxSides, 0,  'TempType','Constant','cstTemp',T_surf);                      % sticky air with constant temperature 


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

% SAVE PARALLEL DATA (parallel)
FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition);




