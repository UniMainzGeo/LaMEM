% LaMEM_Benchmark_3D_Folding
%
% 3D newtonian single-layer folding benchmark using LaMEM.
% This benchmark takes a relatively long amount of time to perform and is
% therefore not included in the standard MATLAB benchmarks

clear

addpath('../../matlab')
CreatePlot = logical(0);
if ismac
    CreatePlot = logical(1);
end
MATLAB_Library=[];
if isunix & ~ismac % on some linux machines, this seems to be required
    MATLAB_Library = 'export LD_LIBRARY_PATH=/usr/lib64; ';
end

n_matrix    =   1;
n_layer     =   1;
R           =   100;        % viscosity contrast
ampl2D      =   0e-5;      % 2D amplitude
ampl3D      =   1e-3;      % 3D amplitude

lam_vec     =   [5:20:65];
% lam_vec     =   [5:1:100];

% lam_vec = 10;

for ilam=1:length(lam_vec)
    for jlam=1:length(lam_vec)
        
        % Wavelength
        lambda_x = lam_vec(ilam);
        lambda_y = lam_vec(jlam);
        
        % Run LaMEM with the correct parameters
        %
        % Note that solver options etc. are specified in growthrate_test_3D.dat
        %
        if isunix
            system([MATLAB_Library, ' ../../bin/LaMEM -ParamFile ./input/SingleLayerFolding_growthrate_test_3D.dat  -MuMeanMethod 0 -W ',num2str(lambda_x), ' -L ',num2str(lambda_y),...
                ' -FoldingBenchmark_R ',num2str(R),...
                ' -FoldingBenchmark_nL ',num2str(n_layer),...
                ' -FoldingBenchmark_nM ',num2str(n_matrix),...
                ' -ampl2D ',num2str(ampl2D),...
                ' -ampl3D ',num2str(ampl3D),...
                ' -NonlinearIterationsAccuracy  1e-5']);
        end
        
        if ispc
            system(['../../bin/LaMEM -ParamFile ./input/growthrate_test_3D.dat  -MuMeanMethod 0 -W ',num2str(lambda_x), ' -L ',num2str(lambda_y),...
                ' -FoldingBenchmark_R ',num2str(R),...
                ' -FoldingBenchmark_nL ',num2str(n_layer),...
                ' -FoldingBenchmark_nM ',num2str(n_matrix),...
                ' -ampl2D ',num2str(ampl2D),...
                ' -ampl3D ',num2str(ampl3D),...
                ' -NonlinearIterationsAccuracy  1e-5']);
        end
        
        
        %======================================================================
        % Read LaMEM output
        fname = './Timestep_000001/SingleLayerFolding_3D_Growthrate';
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
        
        Growthrate(ilam, jlam)  =   maxVz/ampl3D;
       
        Wavex_2D(ilam, jlam)  =    lambda_x;
        Wavey_2D(ilam, jlam)  =    lambda_y;
        
        
    end
end


save test



%--------------------------------------------------------------------------
% COMPUTE ANALYTICAL SOLUTION [reproduces fig 1a of Kaus and Schmalholz
% 2006, GRL]
lam_ana = [0.01:.5:70];
[Wavex_2Da,Wavey_2Da] = meshgrid (lam_ana, lam_ana);
H                   = 1.0;

kx                   = 2.0*3.14156./Wavex_2Da;
ky                   = 2.0*3.14156./Wavey_2Da;
k                   = sqrt(kx.^2 + ky.^2);
str_xx              = 0.75;
str_yy              = 0.25;
str_zz              = -(str_xx+str_yy);

% See Kaus & Schmalholz 2006 for analytical solution
inv_nu          = 1/R
q               = -4*(1-inv_nu)*k./[ 2*k*(1-inv_nu^2) - (1+inv_nu).^2.*exp(k) + (1-inv_nu)^2*exp(-k) ]*str_zz;

dAdt            = (q/2).*[ (ky.^2./k.^2)*(-str_xx/str_zz-1) - (kx.^2./k.^2)*-1*str_xx/str_zz - 1 ]*str_zz;
q_anal          = dAdt/str_zz;


%
%--------------------------------------------------------------------------





%--------------------------------------------------------------------------
%
% Compute error
Gr_ana_num = interp2(Wavex_2Da,Wavey_2Da,q_anal, Wavex_2D,Wavey_2D);

TotalError  =   norm(Gr_ana_num(:)-Growthrate(:));
AllowedError    =   2.3;

clc
disp(['==============================================='])
if TotalError<AllowedError
    disp(' Single layer 3D Folding benchmark PASSED ' )
else
    disp(' Single layer 3D Folding benchmark FAILED ' )
end
disp(['==============================================='])
%
%--------------------------------------------------------------------------




%--------------------------------------------------------------------------
if CreatePlot
    figure(1), clf
    
    % analytics
%     [cc,h]=contour(Wavex_2Da,Wavey_2Da,q_anal,[0:2.5:25]);
    [cc,h]=contour(Wavex_2D,Wavey_2D,Gr_ana_num,[0:2.5:25]);
    
    
    set(h,'LineStyle','-','linecolor','k')
    clabel(cc,h)
    xlabel('\lambda_x/H')
    ylabel('\lambda_y/H')
    title(['Growthrate (q/\epsilon_{zz})   \epsilon_{xx}/\epsilon_{zz}=',num2str(str_xx/str_zz),' \epsilon_{yy}/\epsilon_{zz}=',...
        num2str(str_yy/str_zz),'  \mu_{layer}/\mu_{matrix}=',num2str(R),' error=',num2str(TotalError),' allowed error=',num2str(AllowedError)])
    hold on
    
    
    % numerics
    
    [cc,h]=contour(Wavex_2D, Wavey_2D, Growthrate,[0:2.5:25]);
    set(h,'LineStyle','--','linecolor','r')
    
    legend('Analytical','Numerical')
    
      print('-djpeg',['Benchmark_SingleLayerFolding_3D_',date])
    
end
%--------------------------------------------------------------------------


