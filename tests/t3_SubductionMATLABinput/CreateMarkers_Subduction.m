function CreateMarkers_Subduction(NumberCores)


if ismac
    addpath ../../../LaMEM/matlab
else
    addpath ~/LaMEM/matlab/
    
end
%==========================================================================
% OUTPUT OPTIONS
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;

% Output a single file containing particles information for LaMEM (msetup = redundant)
LaMEM_Redundant_output = 0;%0;

% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
% WARNING: Need a valid 'Parallel_partition' file!
LaMEM_Parallel_output   = 1;

% Mesh from file 1-YES (load uniform or variable mesh from file); 0-NO (create new uniform mesh)
% WARNING: Need a valid 'Parallel_partition' file!
LoadMesh                =   1;

% random noise of particles
RandomNoise             =   logical(0);

if      NumberCores==1
    Parallel_partition  =   'ProcessorPartitioning_1cpu_1.1.1.bin'
    
elseif  NumberCores==4
    Parallel_partition  =   'ProcessorPartitioning_4cpu_4.1.1.bin'
    
end

Is64BIT                 =   logical(0);
%==========================================================================

%==========================================================================
% DOMAIN PARAMETERS (DIMENSIONAL)
%==========================================================================
% Number of markers in a grid cell
npart_x = 3;
npart_y = 3;
npart_z = 3;

% Geometry- related parameters
ThickCrust          =  10;        % or thickness crust
ThickWL             =  30;
ThermalAge_Myrs     =  30;

% Slab parameters
w_op                =   310;
w_max_op            =   1750;
w_min_op            =   -1000;
Depth_slab          =   -250;

%==========================================================================
% MESH GRID
%==========================================================================
if LoadMesh == 1
    [~,~,~,xcoor,ycoor,zcoor,Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition,RandomNoise,Is64BIT);
    X       = Xpart;
    Y       = Ypart;
    Z       = Zpart;
    % Update other variables
    nump_x  = size(X,2);
    nump_y  = size(X,1);
    nump_z  = size(X,3);
    
    z_bot   =   min(Z(:));
    H       =   max(Z(:))-min(Z(:));
    W       =   max(X(:))-min(X(:));
    L       =   max(Y(:))-min(Y(:));
    
    
end

%==========================================================================
% PHASES
%==========================================================================
crust       	=   0;
UC1             =   1; % oceanic plate
UC2             =   2;
UC3             =   3;% mantle

Phase           =   zeros(size(X)); % initialize phases

T_surface       =   20;
kappa           =   1e-6;
TLAB            =   1300;
T_M             =   1280;

depth           =   abs(zcoor);
ThermalAge     	=   ThermalAge_Myrs*1e6*(365*24*3600);
T               =   1350* erf(depth*1000/2/sqrt(kappa*ThermalAge)) +T_surface;

ind             =   find(T<1300);
ThickOP         =   max(depth(ind));


z_surface       =   z_bot+H;
z_UC0           =   -ThickCrust;
z_UC1           =   z_UC0-ThickOP;
z_UC2           =   z_UC1-ThickWL;




%==========================================================================
% SETUP GEOMETRY
%==========================================================================

%%% Mantle
ind         =   find( Z< z_surface);
Phase(ind)  =   3;

% Weak Layer
ind         =   find( X>=w_min_op & X<=w_max_op & Z>=z_UC2 & Z<z_surface);
Phase(ind)  =   2;

alpha       =   -34*pi/180;
R           =   [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

coord_WL    =   [-w_op  0 0 -w_op -w_op; z_surface z_surface z_UC2  z_UC2 z_surface];
coord_WL    =   R*coord_WL;

Poly_x      =   coord_WL(1,:) + w_min_op;
Poly_z      =   coord_WL(2,:) ;
ind         =   find(inpolygon(X,Z,Poly_x,Poly_z));
Phase(ind)  =   2;

% Oceanic plate
ind         =   find( X>=w_min_op & X<=w_max_op & Z>=z_UC1 & Z<=z_surface );
Phase(ind)  =   1;

coord       =   [-w_op  0 0 -w_op -w_op; z_surface z_surface z_UC1  z_UC1 z_surface];
coord       =   R*coord;
Poly_x      =   coord(1,:)+ w_min_op;
Poly_z      =   coord(2,:);
ind       	=   find(inpolygon(X,Z,Poly_x,Poly_z)) ;
Phase(ind)  =   1;

% crust
ind         =   find( X>=w_min_op & X<=w_max_op & Z>=z_UC0 & Z<=z_surface);
Phase(ind)  =   0;

coord       =   [-w_op  0 0 -w_op -w_op; z_surface z_surface z_UC0  z_UC0 z_surface];
coord       =   R*coord;
Poly_x      =   coord(1,:)+ w_min_op;
Poly_z      =   coord(2,:);
ind         =   find(inpolygon(X,Z,Poly_x,Poly_z));
Phase(ind)  =   0;


%==========================================================================
% TEMPERATURE - in Celcius
%==========================================================================
% Temperature at the surface
Temp        =   1350*ones(size(Phase));

% Half cooling method for horizontal part of the slab
ind         =   find(Z>=z_UC2 & X>=w_min_op & X<=w_max_op);  %if you want to have the crust with different gradient just pust smth diff from z_surfaceace
Temp(ind)   =   1350* erf(abs(Z(ind)-z_surface)*1000/2/sqrt(kappa*ThermalAge)) +T_surface;

%Half space cooling for the subducting lithosphere
coord       =   [-w_op  0 0 -w_op -w_op; 0 0 z_UC2  z_UC2 0];
coord       =   R*coord;
Poly_x      =   coord(1,:)+ w_min_op;
Poly_z      =   coord(2,:);
ind         =   find(inpolygon(X,Z,Poly_x,Poly_z));
X_in        =   X(ind);
Z_in        =   Z(ind);
Xc_in       =   [X_in(:)'; Z_in(:)'];

R1          =   [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
Xc_tran     =   R1*Xc_in;       % rotate points back

dist        =   -(Xc_tran(2,:)-max(Xc_tran(2,:)));
Temp(ind)   =   1350.*erf((dist*1000./(2.*sqrt(kappa*ThermalAge))))+ T_surface;


%==========================================================================
% PREPARE DATA FOR VISUALIZATION/OUTPUT
%==========================================================================

% Prepare data for visualization/output
A = struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);

Phase       = permute(Phase,[2 1 3]);
Temp        = permute(Temp, [2 1 3]);

% Linear vectors containing coords
x        = X(1,:,1);
y        = Y(:,1,1);
z        = Z(1,1,:);
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
clearvars -except A Paraview_output LaMEM_Parallel_output Parallel_partition Is64BIT NumberCores

% PARAVIEW VISUALIZATION
if (Paraview_output == 1)
    FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option
end

% SAVE PARALLEL DATA (parallel)
if (LaMEM_Parallel_output == 1)
    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition,Is64BIT);
    
    
    if      NumberCores==1
        !mv -f MatlabInputParticles/ MatlabInputParticles_p1
        
    elseif  NumberCores==4
        !mv -f MatlabInputParticles/ MatlabInputParticles_p4
        
    end
end



