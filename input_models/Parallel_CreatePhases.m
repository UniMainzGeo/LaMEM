%% CreatePhases
clear

Paraview_output       = 1;          % See model setup in Paraview 1-YES; 0-NO

LaMEM_Parallel_output = 1;          % Output parallel files containing particles information for LaMEM, using a processor distribution file
Parallel_partition    = 'ProcessorPartitioning_8cpu_4.1.2.bin';

LaMEM_OLD_WAY_output  = 0;          % Output a single file containing particles information for LaMEM

% Domain
W       =   3000e3;                 % Width in x-dir
L       =   1000e3;                 % Length in y-dir
H       =   800e3;                  % Height in z-dir

% Number of particles
nump_x  =   450;  %for simulations
nump_y  =   150;
nump_z  =   120;

% No of particles in a grid cell
npart_x = 3;
npart_y = 3;
npart_z = 3;

margin  =   50e3;               %specify whether the plates are attached margin=0 or unattached to the boundaries
dx_reg  =   W/(nump_x-1);
dy_reg  =   L/(nump_y-1);
dz_reg  =   H/(nump_z-1);
x_left  =   -1500e3;
y_front =   0;
z_bot   =   0;

[X,Y,Z] =   meshgrid([x_left:dx_reg:x_left+W], [y_front:dy_reg:y_front+L], [z_bot:dz_reg:z_bot+H]);

Phase   =   zeros(size(X));     %   Contains phases
Temp    =   zeros(size(X));     %   Contains temperatures

if round(nump_x/npart_x)~=nump_x/npart_x
    disp(['watch it: probably mistake with number of particles in x-direction'])
end
if round(nump_y/npart_y)~=nump_y/npart_y
    disp(['watch it: probably mistake with number of particles in y-direction'])
end
if round(nump_z/npart_z)~=nump_z/npart_z
    disp(['watch it: probably mistake with number of particles in z-direction'])
end

%==========================================================================
% Set initial phase distribution
%
%Slab
H_slab              =   H/8;
Width_slab          =   L; %2*.9*L;
Depth_slab          =   0.75*H;
Length_slab         =   .35*W;
Depth_LowerMantle   =   670e3;
ThicknessAir        =   0.05*H;

%India block
H_India             =   1.2*H_slab;
Width_India         =   .15*W;
Length_India        =   .15*W;

%Lithosphere_Asia
H_Asia              =   H/8;

%Thickness of crust Asia
Crust_Asia          =   30e3;

%WEAK ZONE
Width_weak          =   0e3;

%PHASES
mantle     = 0;
air        = 1;
slab       = 2;
sed_layer  = 3;
india      = 4;
asia       = 5;
crust_asia = 6;


Thickness_sed = 10e3;
Lith_thickness= H_slab-ThicknessAir-Thickness_sed;

%SLAB
ind         =   find( X>(x_left+margin) & X<0 & Z>(H-H_slab)  );% & Y>margin & Y<(L-margin));
Phase(ind)  =   slab;

ind         =   find( (Z<(H-(X-0/2))) &  (Z>((H-H_slab)-(X-0/2))) & Z>Depth_slab );% & Y>margin & Y<(L-margin) );
Phase(ind)  =   slab;

%Sedimentary layer
ind         =   find( X>(x_left+margin) & X<Thickness_sed & Z>(H-ThicknessAir-Thickness_sed) );% & Y>margin & Y<(L-margin));
Phase(ind)  =   sed_layer;

ind         =   find((Z<(H-(X-0/2))) & Z>Depth_slab & (Z>((H-ThicknessAir)-(X-0/2))) );% & Y>margin & Y<(L-margin));
Phase(ind)  =   sed_layer;

%INDIA
%lithosphere India
ind         =   find( X>-(Length_slab/2+Width_India)/2 & X<-(Length_slab/2-Width_India)/2 & Z>(H-H_India) & Y>(L-Width_India)/2 & Y<(L+Width_India)/2); %& (abs(Y)+y_front)<=Width_slab/2 %%half-slab
Phase(ind)  =   slab;

%crust India
ind         =   find( X>-(Length_slab/2+Width_India)/2 & X<-(Length_slab/2-Width_India)/2 & Z>(H-H_India+Lith_thickness) & Y>(L-Width_India)/2 & Y<(L+Width_India)/2); %& (abs(Y)+y_front)<=Width_slab/2 %%half-slab
Phase(ind)  =   india;

%ASIA
ind         =   find((X>Width_weak) & (Z>(H-(X-Width_weak))) & Z>(H-H_Asia) & X>0 );% & Y>margin & Y<(L-margin) & X<(-x_left-margin));
Phase(ind)  =   asia;

%crust Asia
ind         =   find((X>Width_weak) & (Z>(H-(X-Width_weak))) & Z>(H-ThicknessAir-Crust_Asia) & X>0 );% & Y>margin & Y<(L-margin) & X<(-x_left-margin));
Phase(ind)  =   crust_asia;

% Add AIR
ind         =   find( Z>(H-ThicknessAir) );
Phase(ind)  =   air;


% Add lower mantle in a primitive way
ind         =   find(Z<(H-Depth_LowerMantle)) ;
Phase(ind)  =   mantle;

%==========================================================================

%==========================================================================
% Set initial temperature distribution - in Celcius
%
Temp = (H-Z)./H*0 + 0.5 + (rand(size(Z))-0.5)*0.05 + 0*1000;

%==========================================================================

% Prepare data for visualization/output
A = struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[],'CL',[]);
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
A.CL     = 1000e3;

% SAVE DATA IN 1 FILE (old way)
if (LaMEM_OLD_WAY_output == 0)
    PhaseOrig   = Phase;
    PhaseVec(1) = nump_z;
    PhaseVec(2) = nump_y;
    PhaseVec(3) = nump_x;
    PhaseVec    = [PhaseVec(:); Phase(:); Temp(:)];
    
    % Save data to file
    ParticleOutput  =   'ParticlesInput3D.dat';
    PassiveOutput   =   'PassiveInput.dat';
    PetscBinaryWrite(ParticleOutput, PhaseVec);
end

% % Plot in matlab - not recommended - takes a lot of memory and it's slow
% % Use Paraview visualization option below
% figure(1), clf, hold on
% marker = {'r','k','g','b','m','y'};
% for phase=0:max(Phase(:))
%     ind = find(Phase==phase);
%     ind = find(abs(Phase-phase)<1e-9);
%     ind_plot = phase+1;
%     if ind_plot>length(marker); ind_plot = ind_plot-fix(ind_plot/length(marker))*length(marker)+1; end
%     plot3(X(ind),Y(ind), Z(ind),[marker{ind_plot},'.']);
% end
% view(3)
% axis equal

clearvars -except A Paraview_output LaMEM_Parallel_output Parallel_partition

% PARAVIEW VISUALIZATION
if (Paraview_output == 1)
    
    WriteLaMEMSetup_Matlab2VTK(A,'BINARY'); % fast and should be the default option
    %WriteLaMEMSetup_Matlab2VTK(A,'ASCII'); % for debugging - takes a while
end

% SAVE PARALLEL DATA
if (LaMEM_Parallel_output == 1)
    SaveParticlesParallelFromMatlab(A,Parallel_partition);
end
