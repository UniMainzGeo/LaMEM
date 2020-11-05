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
    disp(['Creating setup in parallel using LaMEM file: ', LaMEM_input_file])
      
    % Define parallel partition file
    Parallel_partition                          =   'ProcessorPartitioning_4cpu_4.1.1.bin'
    
    [Grid,X,Y,Z,npart_x,npart_y,npart_z,W,L,H]  =   LaMEM_ParseInputFile(LaMEM_input_file);
    
    % Load grid from parallel partitioning file
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
Length_Subduct_Slab =  300;     % length of subducted slab
Length_Horiz_Slab   =  1500;    % length of overriding plate of slab
Width_Slab          =  1000;    % Width of slab (in case we run a 3D model)         

SubductionAngle     =   34;     % Subduction angle
ThermalAge_Myrs     =   50;     % Thermal age of the slab in Myrs
ThicknessCrust      =   10;
z_surface           =   0;      % initial free surface


%==========================================================================
% DEFINE SLAB
%==========================================================================

%% 1) Create the top (and bottom) of the slab

dx              =   (Length_Horiz_Slab+Length_Subduct_Slab)/100;

Slab_Top        =   ndgrid(Trench_x_location-Length_Subduct_Slab:dx:Trench_x_location+Length_Horiz_Slab); 
Slab_Top        =   unique([Slab_Top; Trench_x_location; Trench_x_location+Length_Horiz_Slab])'; % ensure that the Trench is included
Slab_Top(2,:)   =   ones(size(Slab_Top))*z_surface;                           % top of slab
ind_T           =   find(Slab_Top(1,:)==Trench_x_location);                   % trench location

Slab_Bot        =   Slab_Top;
Slab_Bot(2,:)   =   -500;

% Rotate inclined piece of slab
alpha           =   -SubductionAngle*pi/180;
R               =   [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

% Top & Bottom of slab
SlabInclined        =   [Slab_Top(:,1:ind_T), Slab_Bot(:,1:ind_T)];
SlabInclined(1,:)	=   SlabInclined(1,:)-Trench_x_location;
SlabInclined        =   R*SlabInclined;
SlabInclined(1,:) 	=   SlabInclined(1,:)+Trench_x_location;

Slab_Top(:,1:ind_T) = SlabInclined(:,1:ind_T);
Slab_Bot        =   [SlabInclined(:,ind_T+1:end), [Slab_Top(1,end); SlabInclined(2,end)]];


%% 2) Compute the perpendicular distance of all "slab" points to top of slab
X2d             =   squeeze(X(1,:,:));
Z2d             =   squeeze(Z(1,:,:));

Distance        =   ones(size(X))*NaN;      % 3D matrix with distance to 
Dist_2D         =   ones(size(X2d))*NaN;    % in x-z plane

[d_min]         =   p_poly_dist(X2d(:), Z2d(:), Slab_Top(1,:), Slab_Top(2,:), false);
Dist_2D(find(Dist_2D)) = d_min;

SlabPolygon     =   [Slab_Top, Slab_Bot(:,end:-1:1)];
in              =   inpolygon(X2d,Z2d,SlabPolygon(1,:),SlabPolygon(2,:));
Dist_2D(~in)    =   NaN;

for iy=1:size(X,1);
    Distance(iy,:,:) = Dist_2D;
end


%% Set phases and temperature based on distance of top of slab
T_surface   =   20;
T_mantle    =   1350;
Phase       =   ones(size(X)); % initialize to have mantle phase
Temp        =   T_mantle*ones(size(Phase));

% Set air
ind        =    find(Z>z_surface);
Phase(ind) =    0;
Temp(ind)  =    0;

% Set Crust
ind        =    find(Distance<ThicknessCrust);
Phase(ind) =    2;

% Set 3D temperature based on halfspace cooling & distance to top of slab
kappa       =   1e-6;
ThermalAge  =   ThermalAge_Myrs*1e6*(365*24*3600);

% halfspace cooling
ind         =   ~isnan(Distance);
Temp(ind)   =   (T_mantle -T_surface) * erf(abs(Distance(ind))*1000/2/sqrt(kappa*ThermalAge)) + T_surface;

% Set Mantle Lithosphere for mantle points that have temperatures < 1200 Celcius
ind        =    find(Temp<1200 & Phase==1);
Phase(ind) =    3;


% Limit lateral size of slab (in 3D cases)
ind         =   find(Y>Width_Slab & Phase>0);
Phase(ind)  =   1;
Temp(ind)   =   T_mantle;



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




