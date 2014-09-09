function [TotalError, AllowedError] = LaMEM_Benchmark_2D_RayleighTaylor(CreatePlot);
%
% 2D Rayleigh-Taylor benchmark using LaMEM, comparing vs analytical
% solution for both Q2Pm1_global and Q1P0 elements
% 

addpath('../../matlab')


if nargin==0
    CreatePlot = logical(0);
    if ismac
        CreatePlot = logical(1);
    end
end
MATLAB_Library=[];
if isunix & ~ismac % on some linux machines, this seems to be required
    MATLAB_Library = 'export LD_LIBRARY_PATH=/usr/lib64; ';
end

TotalError = 100;

R           =   1e2;              % Viscosity contrast
ampl2D      =   1e-3;           % 2D amplitude
Hi          =   0.5;            % Height of interface
lam_vec     =   [.1:.75:3.5];


for Element=1:2
    
    if      Element==1
        vpt_element = 'Q1P0';
    elseif  Element==2
        vpt_element = 'Q2Pm1_global';
    end
    
    
    for ilam=1:length(lam_vec)
        
        % Wavelength
        lambda = lam_vec(ilam);
        
        
        
        % Run LaMEM with the correct parameters
        %
        % Note that solver options etc. are specified in growthrate_test_2D.dat
        %
        if isunix
            system([MATLAB_Library, '../../bin/LaMEM -ParamFile ./input/RayleighTaylor_test_2D.dat  -MuMeanMethod 0 -W ',num2str(lambda),...
                ' -FoldingBenchmark_R ',num2str(R),...
                ' -ampl2D ',num2str(ampl2D),...
                ' -Hinterface ',num2str(Hi),' -restart 0 -VelocitySolver 1 -ph_tol 1e-12 -vpt_element ',vpt_element]);
        else
            error('code not yet tested on PC')
        end
        
        %======================================================================
        % Read LaMEM output
        fname = './Timestep_000000/RT_Benchmark_Output';
        time_step = 0;
        
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
        
        
        z1D = squeeze(Z(1,1,:));
        [dum,indz] = min(abs(z1D-Hi));
        
        
        %======================================================================
        
        
        
        
        
        %======================================================================
        % Compute growthrate
        maxVz       =  max(max(Vz_3D(:,:,indz)));
        
        Growthrate(ilam,Element)  =   maxVz/ampl2D;
        %
        %======================================================================
        
    end
end



b           = 1/2.0;
mu2          = R;
rho1       = 2;
rho2       = 1;
g           = 1;
hi          =   Hi;


aa          = 2.0*3.14159265*b;
lam_anal    =   .1:.1:1.1*lam_vec(end);
for ilam=1:length(lam_anal);
    lambda          = lam_anal(ilam);
    %     Gr_ana(ilam)    = 4.0*mu/((rho_1-rho_2)*g*b)*((lambda/aa+1/(sinh(aa/lambda)*cosh(aa/lambda)))/((lambda/aa)*(lambda/aa)*tanh(aa/lambda)-1/(sinh(aa/lambda)*cosh(aa/lambda))));
    
    w = 2.*pi./lam_anal(ilam);
    str=0;
    
    
    Gr_ana(ilam)   = 1./2.*(-2.*w.^3.*str.*mu2-w.^3.*hi.*rho2.*g-4.*w.^5.*hi.^3.*mu2.^2.*str+4.*w.^3.*hi.*str.*mu2-4.*w.^5.*hi.^2.*str.*mu2+2.*w.^5.*hi.^4.*mu2.^2.*str-  ...
        4.*w.^5.*hi.^4.*str.*mu2+8.*w.^5.*hi.^3.*str.*mu2+2.*w.^5.*hi.^2.*str-4.*w.^3.*hi.*mu2.^2.*str+2.*w.^5.*hi.^2.*mu2.^2.*str-4.*w.^5.*hi.^3.*str+   ...
        2.*w.^5.*hi.^4.*str-w.*hi.*mu2.*rho2.*g.*cosh(w.*hi).^2+w.*mu2.*rho2.*g.*cosh(w.*hi).^2+w.*hi.*mu2.*rho2.*g+w.*hi.*rho2.*g.*cosh(w.*(-1+hi)).^2-  ...
        4.*w.^3.*hi.^2.*str.*mu2+2.*w.^3.*hi.^2.*mu2.^2.*str-w.^3.*hi.^3.*rho2.*g+rho2.*g.*sinh(w.*hi).*cosh(w.*hi)-2.*w.^3.*mu2.^2.*str.*cosh(w.*hi).^2+ ...
        2.*w.^3.*str.*mu2.*cosh(w.*hi).^2+w.^3.*hi.^3.*mu2.*rho2.*g+2.*w.^3.*hi.^2.*mu2.*cosh(w.*(-1+hi)).^2.*str-...
        rho2.*g.*sinh(w.*hi).*cosh(w.*hi).*cosh(w.*(-1+hi)).^2-w.^3.*hi.^2.*mu2.*rho2.*g+w.^2.*hi.^2.*rho2.*g.*sinh(w.*hi).*cosh(w.*hi)-...
        w.^2.*hi.^2.*sinh(w.*(-1+hi)).*mu2.*rho2.*g.*cosh(w.*(-1+hi))-4.*w.^3.*hi.*str.*mu2.*cosh(w.*hi).^2-2.*w.^3.*hi.^2.*mu2.^2.*str.*cosh(w.*hi).^2+...
        4.*w.^3.*hi.*mu2.^2.*str.*cosh(w.*hi).^2-sinh(w.*(-1+hi)).*mu2.*rho2.*g.*cosh(w.*(-1+hi))+sinh(w.*(-1+hi)).*mu2.*rho2.*g.*...
        cosh(w.*(-1+hi)).*cosh(w.*hi).^2+2.*w.^3.*hi.^2.*rho2.*g-2.*w.^2.*hi.*rho2.*g.*sinh(w.*hi).*cosh(w.*hi)+2.*w.^3.*str.*hi.^2.*cosh(w.*hi).^2.*mu2+...
        w.^2.*rho2.*g.*sinh(w.*hi).*cosh(w.*hi)+2.*w.^3.*mu2.^2.*str-w.*mu2.*rho2.*g+2.*w.^3.*hi.^2.*str-2.*w.^3.*hi.^2.*str.*cosh(w.*(-1+hi)).^2-...
        w.*hi.*rho2.*g)./(w.*(cosh(w.*hi).^2.*w.^2.*hi.^2-cosh(w.*(-1+hi)).^2.*mu2.^2.*cosh(w.*hi).^2+cosh(w.*hi).^2+w.^2.*hi.^2.*mu2.^2+...
        hi.^4.*mu2.^2.*w.^4-2.*w.^2.*hi.^2.*mu2+2.*w.^2.*hi.*mu2.^2.*cosh(w.*hi).^2-w.^2.*hi.^2.*mu2.^2.*cosh(w.*hi).^2+...
        w.^2.*mu2.^2.*hi.^2.*cosh(w.*(-1+hi)).^2+2.*sinh(w.*hi).*mu2.*cosh(w.*(-1+hi)).*sinh(w.*(-1+hi)).*cosh(w.*hi)+mu2.^2.*w.^2+...
        mu2.^2.*cosh(w.*(-1+hi)).^2-mu2.^2.*w.^2.*cosh(w.*hi).^2+hi.^4.*w.^4+hi.^2.*w.^2-2.*hi.^3.*w.^4+hi.^2.*w.^4-2.*cosh(w.*hi).^2.*w.^2.*hi-...
        hi.^2.*w.^2.*cosh(w.*(-1+hi)).^2+cosh(w.*hi).^2.*w.^2-cosh(w.*hi).^2.*cosh(w.*(-1+hi)).^2-2.*w.^2.*hi.*mu2.^2-2.*hi.^2.*w.^4.*mu2-...
        2.*hi.^4.*mu2.*w.^4+hi.^2.*mu2.^2.*w.^4-2.*hi.^3.*mu2.^2.*w.^4+4.*hi.^3.*w.^4.*mu2+2.*mu2.*w.^2.*hi));
    
    
end






Gr_ana_num      =   interp1(lam_anal,Gr_ana,lam_vec);
TotalError      =   0;
for i=1:size(Growthrate,2)
    TotalError      =   TotalError + norm(Gr_ana_num(:)-Growthrate(:,i)) ;
end

AllowedError    =   3.5e-5;

clc
disp(['==============================================='])
if TotalError<AllowedError
    disp(' Two-layer Rayleigh-Taylor benchmark PASSED ' )
else
    disp(' Two-layer Rayleigh-Taylor benchmark FAILED ' )
end
disp(['==============================================='])

if CreatePlot
    % Create plot:
    figure(1),clf
    plot(lam_vec,Growthrate(:,1),'o',lam_vec,Growthrate(:,2),'r+',lam_anal,Gr_ana,'k');
    xlabel('\lambda/H')
    ylabel('q')
    legend('Numerical Q1P0','Numerical Q2Pm1\_global','Analytical','Location','SouthEast')
    title(['RT test, Viscosity contrast =',num2str(R),' H_{interface}=',num2str(Hi),' Error=',num2str(TotalError),' Allowed error=',num2str(AllowedError)])
    
    print('-djpeg',['Benchmark_RayleighTaylor_NoSlipUpperLowerBC_',date])
end



