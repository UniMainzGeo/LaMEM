function [TotalError, AllowedError ] = LaMEM_Benchmark_MultilayerFolding(CreatePlot, PertType)

%Multilayer_Numerical_vs_Analytical.m
%Compares the results of the quasi-2D setup with the multilayer analytical
%solution of Boris (Biot_nlayer_powerlaw_time.
%Different length (all with 1 m amplitude cos perturbation) setups are run
%1 time step and the Vz maximum velocity at the salt/overburden interface
%is extracted.
%The cos perturbation can be applied, only to the salt/overburden interface
%(PertType = 0) or to all the interfaces (PertType=1)
% The results are compared with the Growth Rate - Lambda graph of the
%analytical solution. Computes the error at the calculation points

%Required files:
%   DetachmentFolding_Zagros_Cos_Amplitude1.dat
%   Multilayer_Analytical_Inputs.m (a modyfication of DetachmentFolding_7layersOverburden_InputForOptimization.m)
%   Biot_nlayer_powerlaw_time.m

clf, %clear



addpath('../../matlab');  % add I/O routines

if nargin==0
    CreatePlot = logical(0);
    if ismac
        CreatePlot = logical(1);
    end
    PertType = 0; % 0 = Perturbation only on the interface between salt & overburden; 1 = Perturbation in all the layers
    
end

MATLAB_Library=[];
if isunix & ~ismac % on some linux machines, this seems to be required
    MATLAB_Library = 'export LD_LIBRARY_PATH=/usr/lib64; ';
end


e_bg = 1e-15;

%Choose Perturbation Type

BCLower = 8; % 1 = Free Slip; 8 = No slip and imposed velocity (used when the bottom topography is irregular

% Define vector containin Lambda Value(s)
if PertType == 0
    lam_vec= 10e3:5e3:50e3;
    Ampl = [0 1 0 0 0 0 0 0 0];
    lam_vec_analytical	=   0.25:0.1:50;
elseif PertType == 1
    lam_vec=1e3:1e3:10e3;
    Ampl = [0 1 1 1 1 1 1 1 1];
    lam_vec_analytical	=   0.1:0.05:25;
end


for ilam=1:length(lam_vec);
    lambda = lam_vec(ilam);
    system('rm *.breakpoint');
    
    % with direct solver
    str_command = [MATLAB_Library, ' ../../bin/LaMEM -ParamFile ./input/Benchmark_Multilayer_Folding_Zagros_Cos_Amplitude1.dat   -BC.LowerBound ',num2str(BCLower), ' -Perturbation ', num2str(PertType) ,'  -W  ', num2str(lambda),' -BC.Exx ',num2str(e_bg)];
    
    system(str_command);
    
    %======================================================================
    % READ MATLAB OUTPUT
    %time_step=(itime-1)*5;
    Directory = (['./Timestep_000000']);
    fname='DetachmentFolding3D_Zagros_Output';
    
    %newest version
    [info,coord_x,coord_y,coord_z, Vx, Vy, Vz, mu, rho, n,G,...
        C,phi,k,Cp,Q,alpha,FK,intpx, intpy, intpz, charac, temp, Phases,...
        P,Tau_2nd, E2nd, Txx,Tyy,Tzz,Txy,Tyz,Txz, Exx,Eyy,Ezz,Exy,Eyz,Exz, TimeDependentData, NumParticles, ...
        Strain, PlasticStrain]  =   PetscBinaryRead([Directory,'/',fname,'.0.1000000.out']);
    
    % (1) Deal with nodal based data
    xs = info(19);    ys = info(20);    zs = info(21);    xm = info(22);    ym = info(23);    zm = info(24);
    % Shape local data into 3D arrays
    coord_x_3D_loc  = reshape(coord_x,[xm ym zm]);   coord_y_3D_loc = reshape(coord_y,[xm ym zm]);   coord_z_3D_loc = reshape(coord_z,  [xm ym zm]);
    Vx_loc          = reshape(Vx,     [xm ym zm]);   Vy_loc         = reshape(Vy,     [xm ym zm]);   Vz_loc         = reshape(Vz,       [xm ym zm]);
    T_loc           = reshape(temp,   [xm ym zm]);
    % Put local 3D arrays in global 3D arrays
    coord_x_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_x_3D_loc;
    coord_y_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_y_3D_loc;
    coord_z_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_z_3D_loc;
    Vx_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)         = Vx_loc;
    Vy_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)         = Vy_loc;
    Vz_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)         = Vz_loc;
    T_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)          = T_loc;
    CPU_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)        = T_loc*0 + (i-1);
    
    
    % Find the deformed interface in an automated manner
    for iz=1: size(coord_z_3D,3)
        z_level = squeeze(coord_z_3D(:,:,iz));
        z_level_perturbation(iz) = (max(z_level(:))-min(z_level(:)));
        
    end
    
    [Amplitude,ind_interface] = max(z_level_perturbation/2);
    Amplitude = (z_level_perturbation(13)/2);
    ind_interface = 13;
    
    %======================================================================
    
    % extract thge vertical velocity at this interface
    Vz_int= Vz_3D(:,:,ind_interface);
    x_int = coord_x_3D(:,:,ind_interface);
    y_int = coord_y_3D(:,:,ind_interface);
    z_int = coord_z_3D(:,:,ind_interface);
    
    Vz_kin = z_int*e_bg;        % kinematic velocity
    Vz_act = Vz_int - Vz_kin; % active velocity field
    
    % compute growthrate from numerical code
    nx=size(x_int,1);
    
    q_num(ilam) = (max(Vz_act( (nx-1)/2+1 ))/Amplitude)/e_bg;
    
    
end

% Bottom BC for the analytical solution
if BCLower == 1
    LowerBC = 1;
elseif BCLower == 8
    LowerBC = 2;
end

% compute analytical solution
[misfit, lam_vec_an, q_vec_an, l_c] = Multilayer_Analytical_Inputs([22 20 21 20 22 20 22], Ampl, lam_vec_analytical, LowerBC);


q_an_comp   =   interp1(lam_vec_an*l_c, q_vec_an, lam_vec);
TotalError  =   norm(q_an_comp - q_num);
if PertType==0
    AllowedError=   2.8;
elseif PertType==1
    AllowedError=   3.7;
end

clc
disp(['===================================================================='])
if TotalError<AllowedError
    disp([' Multilayer Folding benchmark with perturbation type ', num2str(PertType) ,' PASSED '] )
else
    disp([' Multilayer Folding benchmark with perturbation type ', num2str(PertType) ,' FAILED '] )
end
disp(['===================================================================='])



if CreatePlot
    figure(1), clf
    % Create plot:
    % Plot both analytical and numerical
    plot(lam_vec_an*l_c/1e3, q_vec_an,'k-',lam_vec/1e3,q_num,'ob')
    
    
    xlabel('Wavelength [km]')
    ylabel('q / e_{bg} [ ]')
    if PertType==0
        axis([0 55 0 60])
        title(['multilayer folding benchmark with perturbation @ salt-sediments interface; error=',num2str(TotalError), ' and an allowed error=',num2str(AllowedError)])
    elseif PertType==1
        axis([0 20 0 30])
        title(['multilayer folding benchmark with perturbation @ all interfaces; error=',num2str(TotalError), ' and an allowed error=',num2str(AllowedError)])
    end
    
    
    legend('Numerical Q2Pm1\_global','Analytical','Location','SouthEast')
    
    
    if PertType==0
        print('-djpeg',['Benchmark_MultilayerFolding_SaltOverburdenPerturbed_',date])
    elseif PertType==1
        print('-djpeg',['Benchmark_MultilayerFolding_AllLayersPerturbed_',date])
    end
end





