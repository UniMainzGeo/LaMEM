function [TotalError, AllowedError ] = LaMEM_Benchmark_2D_SingleLayer_Folding(CreatePlot, n_matrix, n_layer)
%
% 2D powerlaw single-layer powerlaw folding benchmark using LaMEM
%
%

% clear

close all

addpath('../../matlab')

if nargin==0
    CreatePlot = logical(0);
    if ismac
        CreatePlot = logical(1);
    end
    n_matrix = 1;
    n_layer  = 1;
    
end
MATLAB_Library=[];
if isunix & ~ismac % on some linux machines, this seems to be required
    MATLAB_Library = 'export LD_LIBRARY_PATH=/usr/lib64; ';
end



% n_matrix    =   1;
% n_layer     =   1;

R           =   100;        % viscosity contrast
ampl2D      =   1e-3;       % 2D amplitude

lam_vec     =   [5:5:30];

for ilam=1:length(lam_vec)
    
    % Wavelength
    lambda = lam_vec(ilam);
    
    
    
    % Run LaMEM with the correct parameters
    %
    % Note that solver options etc. are specified in growthrate_test_2D.dat
    %
    if isunix
        system([MATLAB_Library,'../../bin/LaMEM -ParamFile ./input/growthrate_test_2D.dat  -MuMeanMethod 0 -W ',num2str(lambda),...
            ' -FoldingBenchmark_R ',num2str(R),...
            ' -FoldingBenchmark_nL ',num2str(n_layer),...
            ' -FoldingBenchmark_nM ',num2str(n_matrix),...
            ' -ampl2D ',num2str(ampl2D),...
            ' -NonlinearIterationsAccuracy  1e-5 -restart 0 -save_breakpoints 0']);
    elseif ispc
        system([MATLAB_Library,'../../bin/LaMEM -ParamFile ./input/growthrate_test_2D.dat  -MuMeanMethod 0 -W ',num2str(lambda),...
            ' -FoldingBenchmark_R ',num2str(R),...
            ' -FoldingBenchmark_nL ',num2str(n_layer),...
            ' -FoldingBenchmark_nM ',num2str(n_matrix),...
            ' -ampl2D ',num2str(ampl2D),...
            ' -NonlinearIterationsAccuracy  1e-5 -restart 0 -save_breakpoints 0']);
    end
    
    
    %======================================================================
    % Read LaMEM output
    fname = './Timestep_000001/Growthrate_Output';
    time_step = 1;
    
    [info]   =   PetscBinaryRead([fname,'.',num2str(0),'.',num2str(time_step+1000000),'.out']);
    num_proc =   info(18);  nnode_x  =   info(1);    nnode_y  =   info(2); nnode_z  =   info(3);
    ElementType = info(10); nel_x    =   info(4);    nel_y    =   info(5); nel_z    =   info(6);
    ngp_vel     = info(11);
    ngp_1D      = ngp_vel^(1/3);
    
    coord_x_3D = zeros(nnode_x,nnode_y,nnode_z);    coord_y_3D = zeros(nnode_x,nnode_y,nnode_z);   coord_z_3D = zeros(nnode_x,nnode_y,nnode_z);
    Vx_3D      = zeros(nnode_x,nnode_y,nnode_z);    Vy_3D      = zeros(nnode_x,nnode_y,nnode_z);   Vz_3D      = zeros(nnode_x,nnode_y,nnode_z);
    
    
    for i=1:num_proc  % read files from all CPU's
        
        if (1==1)
            %newest version
            [info,coord_x,coord_y,coord_z, Vx, Vy, Vz, mu, rho, n,G,...
                C,phi,k,Cp,Q,alpha,FK,intpx, intpy, intpz, charac, temp, Phases,...
                P,Tau_2nd, E2nd, Txx,Tyy,Tzz,Txy,Tyz,Txz, Exx,Eyy,Ezz,Exy,Eyz,Exz, TimeDependentData, NumParticles]  =   PetscBinaryRead([fname,'.',num2str(i-1),'.',num2str(time_step+1000000),'.out']);
            
        end
        
        
        % (1) Deal with nodal based data
        xs = info(19);    ys = info(20);    zs = info(21);    xm = info(22);    ym = info(23);    zm = info(24);
        % Shape local data into 3D arrays
        coord_x_3D_loc  = reshape(coord_x,[xm ym zm]);   coord_y_3D_loc = reshape(coord_y,[xm ym zm]);   coord_z_3D_loc = reshape(coord_z,  [xm ym zm]);
        Vx_loc          = reshape(Vx,     [xm ym zm]);   Vy_loc         = reshape(Vy,     [xm ym zm]);   Vz_loc         = reshape(Vz,       [xm ym zm]);
        
        % Put local 3D arrays in global 3D arrays
        coord_x_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_x_3D_loc;
        coord_y_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_y_3D_loc;
        coord_z_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)    = coord_z_3D_loc;
        Vx_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)         = Vx_loc;
        Vy_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)         = Vy_loc;
        Vz_3D(xs+1:xs+xm,ys+1:ys+ym,zs+1:zs+zm)         = Vz_loc;
        
        
    end
    
    X       = coord_x_3D;
    Y       = coord_y_3D;
    Z       = coord_z_3D;
    Exx      = info(12);   Eyy    = info(13); Ezz = -(Exx+Eyy);
    
    % Compute active velocity fields:
    Vz_act  = Vz_3D + Z*Ezz;
    Vy_act  = Vy_3D + Y*Eyy;
    Vx_act  = Vx_3D + X*Exx;
    
    %======================================================================
    
    
    
    
    
    %======================================================================
    % Compute growthrate
    maxVz       =  max(Vz_act(:));
    
    Growthrate(ilam)  =   maxVz/ampl2D;
    %
    %======================================================================
    
end





% Compute analytical growthrate
lam_anal    =   .1:.5:1.1*lam_vec(end);
nL          =   n_layer  + 1e-6;
nM          =   n_matrix + 1e-6;
for ilam=1:length(lam_anal);
    lambda          = lam_anal(ilam);
    Gr_ana(ilam)    = 2.0*nL*(1.0-(1.0/R))/(-1.0+(1/(R*R)*nL/nM)+(sqrt((nL-1.0)))*((1.0+1.0/R*(sqrt((nL/nM))))*(1.0+1.0/R*(sqrt((nL/nM))))*(exp(2.0*3.14156/lambda*(sqrt(1.0/nL))))-(1.0-1.0/R*(sqrt((nL/nM))))*(1.0-1.0/R*(sqrt((nL/nM))))*(exp(-2.0*3.14156/lambda*(sqrt(1.0/nL)))))/sin(2.0*(sqrt((1.0-1.0/nL)))*3.14156/lambda)/2.0);
end

Gr_ana_num  =   interp1(lam_anal,Gr_ana,lam_vec);
TotalError  =   norm(Gr_ana_num-Growthrate);
AllowedError    =   2.7;

clc
disp(['==============================================='])
if TotalError<AllowedError
    disp(' Single layer Folding benchmark PASSED ' )
else
    disp(' Single layer Folding benchmark FAILED ' )
end
disp(['==============================================='])





if CreatePlot
    % Create plot:
    figure(1),clf
    plot(lam_vec,Growthrate,'+',lam_anal,Gr_ana);
    xlabel('\lambda/H')
    ylabel('q/\epsilon_{bg}')
    legend('Numerical Q2Pm1\_global','Analytical','Location','SouthEast')
    title(['Viscosity contrast =',num2str(R),' n_{layer}=',num2str(nL),' n_{matrix}=',num2str(nM),' Error=',num2str(TotalError)])
    
     print('-djpeg',['Benchmark_SingleLayerFolding_nMatrix',num2str(n_matrix),'_nLayer',num2str(n_layer),'_',date])
     
end



