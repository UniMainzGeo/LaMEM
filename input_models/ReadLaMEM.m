%ReadStokesMG_1
%
%Reads StokesMG output from disk - also in the case that the output was
%   generated and saved on more than one processor
clear
    

addpath('../matlab');  % add I/O routines


plast = 0;

% read first one to learn how many files there should be
fname = './DetachmentFolding3D_Output';     % gravity & positive salt buoyancy

% fname = 'SaltDiapirsFastErosion'
% fname = 'FallingBlock3D'

IntPointProperties = logical(1);

time_step = input('timestep = ');
if isempty(time_step); time_step=0; end

if plast == 1
    visualization = logical(0);
else
    visualization = logical(1);
end

%time_step=3;
dstep=10;
% for time_step=time_step:dstep:2e3
for time_step=time_step:dstep:2000
    
    time_step
    %    SimName = 'Case1_32by128'
    
    Directory = ['Timestep_',num2str(time_step,'%1.6i')]
    
    % Read all the data, and reconstruct it into 3D arrays======================
    [info]   =   PetscBinaryRead([Directory,'/',fname,'.',num2str(0),'.',num2str(time_step+1000000),'.out']);
    num_proc =   info(18);  nnode_x  =   info(1);    nnode_y  =   info(2); nnode_z  =   info(3);
    ElementType = info(10); nel_x    =   info(4);    nel_y    =   info(5); nel_z    =   info(6);
    ngp_vel     = info(11);
    ngp_1D      = ngp_vel^(1/3);
    
    coord_x_3D = zeros(nnode_x,nnode_y,nnode_z);    coord_y_3D = zeros(nnode_x,nnode_y,nnode_z);   coord_z_3D = zeros(nnode_x,nnode_y,nnode_z);
    Vx_3D      = zeros(nnode_x,nnode_y,nnode_z);    Vy_3D      = zeros(nnode_x,nnode_y,nnode_z);   Vz_3D      = zeros(nnode_x,nnode_y,nnode_z);
    T_3D       = zeros(nnode_x,nnode_y,nnode_z);
    CPU_3D      = zeros(nnode_x,nnode_y,nnode_z);                   % CPU's
    
    mu_3D      = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % viscosity
    rho_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % density
    n_3D       = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % powerlaw exponent
    G_3D       = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Elastic shear module
    C_3D       = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Cohesion
    phi_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Friction angle
    k_3D       = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % thermal conductivity
    Cp_3D      = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % heat capacity
    Q_3D       = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % heat production
    alpha_3D   = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % thermal expansivity
    FK_3D      = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Frank-Kamenetskii parameter
    Phases_3D  = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Phase distribution
    Pressure_3D= zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Pressure
    Tau2nd_3D  = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % 2nd invariant deviatoric stress tensor
    E2nd_3D    = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % 2nd invariant deviatoric strainrate tensor
    Txx_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Txx
    Tyy_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Tyy
    Tzz_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Tzz
    Txz_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Txz
    Txy_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Txy
    Tyz_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Tyz
    Exx_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Exx
    Eyy_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Eyy
    Ezz_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Ezz
    Exz_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Exz
    Exy_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Exy
    Eyz_3D     = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Eyz
    NumPart_3D = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % # Particles per integration point
    Strain_3D           = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % Strain
    PlasticStrain_3D    = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);     % PlasticStrain
    
    
    intpx_3D   = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);
    intpy_3D   = zeros(nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);
    intpz_3D   = zeros( nel_x*ngp_1D,nel_y*ngp_1D,nel_z*ngp_1D);
    
    for i=1:num_proc  % read files from all CPU's
        
        if (1==1)
            %newest version
            [info,coord_x,coord_y,coord_z, Vx, Vy, Vz, mu, rho, n,G,...
                C,phi,k,Cp,Q,alpha,FK,intpx, intpy, intpz, charac, temp, Phases,...
                P,Tau_2nd, E2nd, Txx,Tyy,Tzz,Txy,Tyz,Txz, Exx,Eyy,Ezz,Exy,Eyz,Exz, TimeDependentData, NumParticles, ...
                Strain, PlasticStrain]  =   PetscBinaryRead([Directory,'/',fname,'.',num2str(i-1),'.',num2str(time_step+1000000),'.out']);
            
        else
            [info,coord_x,coord_y,coord_z, Vx, Vy, Vz, mu, rho, intpx, intpy, intpz, charac]  =   PetscBinaryRead([Directory,'/',fname,'.',num2str(i-1),'.',num2str(time_step+1000000),'.out']);
            temp=Vx;    n=mu*0; C=mu*0; G=mu*0; phi=mu*0; k=mu*0; Cp=mu*0; Q=mu*0; alpha=mu*0; FK=mu*0;
            Phases=FK;P=FK;Tau_2nd=FK; E2nd=FK; Txx=FK; Tyy=FK; Tzz=FK; Txz=FK; Tyz=FK; Txy=FK;
            
        end
        
        
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
        
        
        
        % (2) Deal with integration-point data
        xs = info(25);    ys = info(26);    zs = info(27);    xm = info(28);    ym = info(29);    zm = info(30);
        num = 1;
        for iel_z = zs+1:zs+zm
            for iel_y = ys+1:ys+ym
                for iel_x = xs+1:xs+xm
                    mu_loc      = mu((num-1)*ngp_vel+1:(num)*ngp_vel);
                    rho_loc     = rho((num-1)*ngp_vel+1:(num)*ngp_vel);
                    n_loc       = n((num-1)*ngp_vel+1:(num)*ngp_vel);
                    G_loc       = G((num-1)*ngp_vel+1:(num)*ngp_vel);
                    C_loc       = C((num-1)*ngp_vel+1:(num)*ngp_vel);
                    phi_loc     = phi((num-1)*ngp_vel+1:(num)*ngp_vel);
                    k_loc       = k((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Cp_loc      = Cp((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Q_loc       = Q((num-1)*ngp_vel+1:(num)*ngp_vel);
                    alpha_loc   = alpha((num-1)*ngp_vel+1:(num)*ngp_vel);
                    FK_loc      = FK((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Phases_loc  = Phases((num-1)*ngp_vel+1:(num)*ngp_vel);
                    NumPart_loc = NumParticles((num-1)*ngp_vel+1:(num)*ngp_vel);
                    
                    P_loc       = P((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Tau2nd_loc  = Tau_2nd((num-1)*ngp_vel+1:(num)*ngp_vel);
                    E2nd_loc    = E2nd((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Txx_loc     = Txx((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Tyy_loc     = Tyy((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Tzz_loc     = Tzz((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Txz_loc     = Txz((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Txy_loc     = Txy((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Tyz_loc     = Tyz((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Exx_loc     = Exx((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Eyy_loc     = Eyy((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Ezz_loc     = Ezz((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Exz_loc     = Exz((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Exy_loc     = Exy((num-1)*ngp_vel+1:(num)*ngp_vel);
                    Eyz_loc     = Eyz((num-1)*ngp_vel+1:(num)*ngp_vel);
                    
                    Strain_loc  = Strain((num-1)*ngp_vel+1:(num)*ngp_vel);
                    PlasticStrain_loc  = PlasticStrain((num-1)*ngp_vel+1:(num)*ngp_vel);
                    
                    
                    
                    intpx_loc   = intpx((num-1)*ngp_vel+1:(num)*ngp_vel);
                    intpy_loc   = intpy((num-1)*ngp_vel+1:(num)*ngp_vel);
                    intpz_loc   = intpz((num-1)*ngp_vel+1:(num)*ngp_vel);
                    if IntPointProperties
                        mu_loc_3D   = permute(reshape(mu_loc,       [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        rho_loc_3D  = permute(reshape(rho_loc,      [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        n_loc_3D    = permute(reshape(n_loc,        [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        G_loc_3D    = permute(reshape(G_loc,        [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        C_loc_3D    = permute(reshape(C_loc,        [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        phi_loc_3D  = permute(reshape(phi_loc,      [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        k_loc_3D    = permute(reshape(k_loc,        [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Cp_loc_3D   = permute(reshape(Cp_loc,       [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Q_loc_3D    = permute(reshape(Q_loc,        [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        alpha_loc_3D= permute(reshape(alpha_loc,    [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        FK_loc_3D   = permute(reshape(FK_loc,       [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Phases_loc_3D=permute(reshape(Phases_loc,   [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        NumPart_loc_3D=permute(reshape(NumPart_loc, [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        
                        Strain_loc_3D           =   permute(reshape(Strain_loc, [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        PlasticStrain_loc_3D    =   permute(reshape(PlasticStrain_loc, [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        
                        P_loc_3D     =permute(reshape(P_loc,        [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Tau2nd_loc_3D=permute(reshape(Tau2nd_loc,   [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        E2nd_loc_3D  =permute(reshape(E2nd_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Txx_loc_3D   =permute(reshape(Txx_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Tyy_loc_3D   =permute(reshape(Tyy_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Tzz_loc_3D   =permute(reshape(Tzz_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Txz_loc_3D   =permute(reshape(Txz_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Tyz_loc_3D   =permute(reshape(Tyz_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Txy_loc_3D   =permute(reshape(Txy_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Exx_loc_3D   =permute(reshape(Exx_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Eyy_loc_3D   =permute(reshape(Eyy_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Ezz_loc_3D   =permute(reshape(Ezz_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Exz_loc_3D   =permute(reshape(Exz_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Eyz_loc_3D   =permute(reshape(Eyz_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        Exy_loc_3D   =permute(reshape(Exy_loc,     [ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        
                        intpx_loc_3D= permute(reshape(intpx_loc,[ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        intpy_loc_3D= permute(reshape(intpy_loc,[ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                        intpz_loc_3D= permute(reshape(intpz_loc,[ngp_1D ngp_1D ngp_1D]),[3 2 1]);
                    end
                    
                    indx        =  (iel_x-1)*ngp_1D+1:(iel_x)*ngp_1D;
                    indy        =  (iel_y-1)*ngp_1D+1:(iel_y)*ngp_1D;
                    indz        =  (iel_z-1)*ngp_1D+1:(iel_z)*ngp_1D;
                    
                    
                    if IntPointProperties
                        
                        mu_3D (   indx, indy, indz )   = mu_loc_3D;
                        rho_3D(   indx, indy, indz )   = rho_loc_3D;
                        n_3D  (   indx, indy, indz )   = n_loc_3D;
                        G_3D  (   indx, indy, indz )   = G_loc_3D;
                        C_3D  (   indx, indy, indz )   = C_loc_3D;
                        phi_3D(   indx, indy, indz )   = phi_loc_3D;
                        k_3D  (   indx, indy, indz )   = k_loc_3D;
                        
                        Cp_3D   (   indx, indy, indz )    = Cp_loc_3D;
                        Q_3D    (   indx, indy, indz )    = Q_loc_3D;
                        alpha_3D(   indx, indy, indz )    = alpha_loc_3D;
                        FK_3D   (   indx, indy, indz )    = FK_loc_3D;
                        Phases_3D(  indx, indy, indz )    = Phases_loc_3D;
                        NumPart_3D(  indx, indy, indz )    = NumPart_loc_3D;
                        
                        
                        Pressure_3D(  indx, indy, indz  ) = P_loc_3D;
                        Tau2nd_3D(    indx, indy, indz  ) = Tau2nd_loc_3D;
                        E2nd_3D(      indx, indy, indz  ) = E2nd_loc_3D;
                        Txx_3D(       indx, indy, indz  ) = Txx_loc_3D;
                        Tyy_3D(       indx, indy, indz  ) = Tyy_loc_3D;
                        Tzz_3D(       indx, indy, indz  ) = Tzz_loc_3D;
                        Txy_3D(       indx, indy, indz  ) = Txy_loc_3D;
                        Tyz_3D(       indx, indy, indz  ) = Tyz_loc_3D;
                        Txz_3D(       indx, indy, indz  ) = Txz_loc_3D;
                        Exx_3D(       indx, indy, indz  ) = Exx_loc_3D;
                        Eyy_3D(       indx, indy, indz  ) = Eyy_loc_3D;
                        Ezz_3D(       indx, indy, indz  ) = Ezz_loc_3D;
                        Exy_3D(       indx, indy, indz  ) = Exy_loc_3D;
                        Eyz_3D(       indx, indy, indz  ) = Eyz_loc_3D;
                        Exz_3D(       indx, indy, indz  ) = Exz_loc_3D;
                        
                        Strain_3D(       indx, indy, indz  ) = Strain_loc_3D;
                        PlasticStrain_3D(       indx, indy, indz  ) = PlasticStrain_loc_3D;
                        
                        
                        
                        intpx_3D((iel_x-1)*ngp_1D+1:(iel_x)*ngp_1D,(iel_y-1)*ngp_1D+1:(iel_y)*ngp_1D, (iel_z-1)*ngp_1D+1:(iel_z)*ngp_1D )   = intpx_loc_3D;
                        intpy_3D((iel_x-1)*ngp_1D+1:(iel_x)*ngp_1D,(iel_y-1)*ngp_1D+1:(iel_y)*ngp_1D, (iel_z-1)*ngp_1D+1:(iel_z)*ngp_1D )   = intpy_loc_3D;
                        intpz_3D((iel_x-1)*ngp_1D+1:(iel_x)*ngp_1D,(iel_y-1)*ngp_1D+1:(iel_y)*ngp_1D, (iel_z-1)*ngp_1D+1:(iel_z)*ngp_1D )   = intpz_loc_3D;
                    end
                    
                    num = num+1;
                end
            end
        end
        
        
    end
    
    %
    if plast == 1
        EXX         = [EXX Exx_3D(7,7,7)];
        TXX         = [TXX Txx_3D(7,7,7)];
    end
    
    Time_vec        =   TimeDependentData(1 , :);
    Vrms_time       =   TimeDependentData(2 , :);
    SlabDepthTime   =   max(coord_z_3D(:))-TimeDependentData(19, :);
    % =========================================================================
    
    characteristic.cmYear				 = charac(1);
    characteristic.Myrs			     = charac(2);
    characteristic.MPa                 = charac(3);
    characteristic.SecYear             = charac(4);
    characteristic.km                  = charac(5);
    characteristic.ThermalExpansivity	 = charac(6);
    characteristic.Strainrate			 = charac(7);
    characteristic.kg					 = charac(8);
    characteristic.Density			 = charac(9);
    characteristic.Viscosity			 = charac(10);
    characteristic.Temperature		 = charac(11);
    characteristic.Velocity			 = charac(12);
    characteristic.Stress				 = charac(13);
    characteristic.Time				 = charac(14);
    characteristic.Length				 = charac(15);
    time                     = info(31);
    dt                       = info(32);
    
    Time_Tip          = Time_vec./characteristic.SecYear./1e6;
    MaxSlabDepthTime  =   max(coord_z_3D(:))-TimeDependentData(19, :);
    
    
    %figure(1), clf
    %slice(intpx_3D,intpy_3D,intpz_3D,mu_3D,mean(intpx_3D(:)),0.9*max(intpy_3D(:)),max(intpz_3D(:))*0.1);
    %axis equal
    %view(0,0)
    
    
    X               = coord_x_3D;
    Y               = coord_y_3D;
    Z               = coord_z_3D;
    Exx             = info(12);   Eyy    = info(13); Ezz = -(Exx+Eyy);
    GravityAngle    = info(33);
    
    if GravityAngle~=0
        % rotate the box (easier for visualization)
        RotMat              = [cos(-GravityAngle/180*pi) -sin(-GravityAngle/180*pi); sin(-GravityAngle/180*pi), cos(-GravityAngle/180*pi)];
        new                 = RotMat*[X(:)'; Z(:)'];
        X(find(X==X))       = new(1,:);
        Z(find(Z==Z))       = new(2,:);
        
        new                 = RotMat*[Vx_3D(:)'; Vz_3D(:)'];
        Vx_3D(find(X==X)) 	= new(1,:);
        Vz_3D(find(Z==Z))   = new(2,:);
        
        
        new                                 = RotMat*[intpx_3D(:)'; intpz_3D(:)'];
        intpx_3D(find(intpx_3D==intpx_3D)) 	= new(1,:);
        intpz_3D(find(intpx_3D==intpx_3D)) 	= new(2,:);
        
    end
    
    
    
    
    Vz_act = Vz_3D + Z*Ezz;
    Vy_act = Vy_3D + Y*Eyy;
    Vx_act = Vx_3D + X*Exx;
    
    
    if visualization
        figure(2),  clf
        hold on,
        
        col = 'b-';
        
        %     X = coord_x_3D./characteristic.km;
        %     Y = coord_y_3D./characteristic.km;
        %     Z = coord_z_3D./characteristic.km;
        
        
        plot3(squeeze(X(:,1,:)),  squeeze(Y(:,1,:)),  squeeze(Z(:,1,:))  ,col)
        plot3(squeeze(X(:,1,:)).',squeeze(Y(:,1,:)).',squeeze(Z(:,1,:)).',col)
        plot3(squeeze(X(1,:,:)),  squeeze(Y(1,:,:)),  squeeze(Z(1,:,:))  ,col)
        plot3(squeeze(X(1,:,:)).',squeeze(Y(1,:,:)).',squeeze(Z(1,:,:)).',col)
        plot3(squeeze(X(:,:,end)),  squeeze(Y(:,:,end)),  squeeze(Z(:,:,end))  ,col)
        plot3(squeeze(X(:,:,end)).',squeeze(Y(:,:,end)).',squeeze(Z(:,:,end)).',col)
        %surf(X(:,:,end),Y(:,:,end),Z(:,:,end), Vz_act(:,:,end)*characteristic.cmYear), colorbar
        surf(X(:,:,end),Y(:,:,end),Z(:,:,end), Vx_3D(:,:,end)*characteristic.cmYear), colorbar
        %surf(X(:,:,end),Y(:,:,end),Z(:,:,end), Vy_3D(:,:,end)*characteristic.cmYear), colorbar
        
        %surf(X(:,:,end),Y(:,:,end),Z(:,:,end), Z(:,:,end)-mean(mean(Z(:,:,end))))
        
        
        view(3), axis equal, axis tight
        
        figure(3), clf
        subplot(221)
         pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,squeeze(rho_3D(:,fix(end/2),:))), shading flat, colorbar
        axis equal, axis tight, title(['Density @ time=',num2str(time/characteristic.SecYear/characteristic.Myrs)])
        step = 2
        hold on
        quiver(squeeze(X(1:step:end,1,1:step:end)),squeeze(Z(1:step:end,1,1:step:end)),squeeze(Vx_3D(1:step:end,1,1:step:end)),squeeze(Vz_3D(1:step:end,1,1:step:end)),'w');
        %quiver(squeeze(X(1:step:end,1,1:step:end))/characteristic.km,squeeze(Z(1:step:end,1,1:step:end))/characteristic.km,squeeze(Vx_act(1:step:end,1,1:step:end)),squeeze(Vz_act(1:step:end,1,1:step:end)),'w');
        
        
        subplot(222)
         pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,log10(squeeze(mu_3D(:,1,:)))), shading flat, colorbar
        %axis equal, axis tight, title(['log10(Viscosity), max vel=',num2str(max(sqrt(Vz_3D(:).^2+Vx_3D(:).^2 + Vy_3D(:).^2))*characteristic.cmYear)])
        %axis equal, axis tight, title(['log10(Viscosity), max vz=',num2str(max(sqrt(Vz_3D(:).^2))*characteristic.cmYear),' max vx=',num2str(max(sqrt(Vx_3D(:).^2))*characteristic.cmYear)])
        axis tight, title(['log10(Viscosity), max vz=',num2str(max(sqrt(Vz_3D(:).^2))*characteristic.cmYear),' max vx=',num2str(max(sqrt(Vx_3D(:).^2))*characteristic.cmYear)])
        
        hold on
        % quiver(squeeze(X(1:step:end,1,1:step:end)),squeeze(Z(1:step:end,1,1:step:end)),squeeze(Vx_3D(1:step:end,1,1:step:end)),squeeze(Vz_3D(1:step:end,1,1:step:end)),'w');
        %axis equal, axis tight
        
        subplot(223)
        pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),(squeeze(T_3D(:,1,:)))), shading flat, colorbar
        %axis equal, axis tight, title(['log10(Viscosity), max vel=',num2str(max(sqrt(Vz_3D(:).^2+Vx_3D(:).^2 + Vy_3D(:).^2))*characteristic.cmYear)])
        axis equal, axis tight, title(['Temperature'])
        hold on
        quiver(squeeze(X(1:step:end,1,1:step:end)),squeeze(Z(1:step:end,1,1:step:end)),squeeze(Vx_3D(1:step:end,1,1:step:end)),squeeze(Vz_3D(1:step:end,1,1:step:end)),'w');
        
        
        subplot(224)
         pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,(squeeze(Phases_3D(:,1,:)))), shading flat, colorbar
        %axis equal, axis tight, title(['log10(Viscosity), max vel=',num2str(max(sqrt(Vz_3D(:).^2+Vx_3D(:).^2 + Vy_3D(:).^2))*characteristic.cmYear)])
        axis equal, axis tight, title(['Phases'])
        hold on
        quiver(squeeze(X(1:step:end,1,1:step:end))/characteristic.km,squeeze(Z(1:step:end,1,1:step:end))/characteristic.km,squeeze(Vx_3D(1:step:end,1,1:step:end)),squeeze(Vz_3D(1:step:end,1,1:step:end)),'w');
        
        Vel_mean=mean(sqrt(Vx_3D(:).^2+Vz_3D(:).^2))
        
        figure(234), clf
        %pcolor(squeeze(X(:,1,:)),squeeze(Z(:,1,:)),squeeze(sqrt(Vz_3D(:,1,:).^2 + Vx_3D(:,1,:).^2))/Vel_mean), colorbar
        pcolor(squeeze(X(:,1,:))/characteristic.km,squeeze(Z(:,1,:))/characteristic.km,squeeze(Vz_act(:,1,:))), colorbar
        %          pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,log10(squeeze(E2nd_3D(:,1,:))*characteristic.Strainrate)), shading flat, colorbar
        
        %pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,(squeeze(Strain_3D(:,1,:)))), shading flat, colorbar
        %    pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,(squeeze(PlasticStrain_3D(:,1,:)))), shading flat, colorbar
        
        
        shading flat, axis equal
        %         title('log10(E2nd)')
        
        title('Vz_{active}')
        
        hold on
        quiver(squeeze(X(1:step:end,1,1:step:end))/characteristic.km,squeeze(Z(1:step:end,1,1:step:end))/characteristic.km,squeeze(Vx_3D(1:step:end,1,1:step:end)),squeeze(Vz_3D(1:step:end,1,1:step:end)),'k');
        
        figure(235), clf
        subplot(121)
%         pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,(squeeze(Phases_3D(:,end/2,:))))
        %shading flat,
        colorbar
        axis equal
        axis tight
        title('Phase distribution')
        
        
        subplot(122)
%         pcolor(squeeze(intpx_3D(:,1,:))/characteristic.km,squeeze(intpz_3D(:,1,:))/characteristic.km,(squeeze(NumPart_3D(:,end/2,:))))
        %shading flat,
        colorbar
        axis equal
        axis tight
        title('# particles per integration point')
        disp(['Total # of particles=',num2str(sum(NumPart_3D(:)))])
        
    end
    
    
    %============================================================
    if (1==0)
        Part.x   = []; Part.y     = []; Part.z   = [];  Part.num = []; Part.phase = []; Part.cpu = [];
        Part.ix  = []; Part.iy    = []; Part.iz   = []; Part.T   = [];
        
        [info,Partic]    =   PetscBinaryRead(['Particles.',num2str(i-1),'.',num2str(time_step+1000000),'.out']);
        num_particle_prop   =   info(3);
        
        Part.x     =   [Part.x;     Partic( 1:num_particle_prop:end)];
        Part.y     =   [Part.y;     Partic( 2:num_particle_prop:end)];
        Part.z     =   [Part.z;     Partic( 3:num_particle_prop:end)];
        Part.num   =   [Part.num;   Partic( 4:num_particle_prop:end)];
        Part.phase =   [Part.phase; Partic( 5:num_particle_prop:end)];
        Part.cpu   =   [Part.cpu;   Partic( 6:num_particle_prop:end)];
        Part.ix    =   [Part.ix;    Partic( 7:num_particle_prop:end)];
        Part.iy    =   [Part.iy;    Partic( 8:num_particle_prop:end)];
        Part.iz    =   [Part.iz;    Partic( 9:num_particle_prop:end)];
        Part.T     =   [Part.T;     Partic(13:num_particle_prop:end)];
        
        subplot(211)
        Part.x = Part.x*characteristic.Length/characteristic.km;
        Part.y = Part.y*characteristic.Length/characteristic.km;
        Part.z = Part.z*characteristic.Length/characteristic.km;
        
        
        ind = find(Part.phase==0);      %deactivated particles
        plot(Part.x(ind),Part.z(ind),'r+')
        
        ind = find(Part.phase==1);      %deactivated particles
        plot(Part.x(ind),Part.z(ind),'k+')
        
        ind = find(Part.phase==2);      %deactivated particles
        plot(Part.x(ind),Part.z(ind),'y+')
        
    end
    
    
    if (1==0)
        load HeatFlowAndCrustalThicknessData.mat
        
        % Plot Western US data
        X_surf = X(:,:,end);
        Y_surf = Y(:,:,end);
        Z_surf = Z(:,:,end)*1e3;
        
        % Modify modeled data to fit height @ lower left corner
        diff_z = Z_surf(1)-topol_z(1);
        Z_surf=Z_surf-diff_z;
        
        
        dlat         =  (247-235)/(nnode_x-1);
        dlong        =  (42 -29 )/(nnode_y-1);
        [ylong,xlat] = meshgrid([29:dlong:42],[235:dlat:247]);
        figure(333)
        subplot(121)
        pcolor(qx,qy,topol_z); shading flat; caxis([-4500 4000]), colorbar
        title('Observed topography')
        states = shaperead('usastatehi', 'UseGeoCoords', true);
        Ca     = states(5);
        Ca.Lon = Ca.Lon+360;
        hold on
        plot(Ca.Lon,Ca.Lat,'k')
        
        subplot(122)
        pcolor(xlat,ylong,Z_surf); shading flat; caxis([-4500 4000])
        colorbar
        states = shaperead('usastatehi', 'UseGeoCoords', true);
        Ca     = states(5);
        Ca.Lon = Ca.Lon+360;
        hold on
        plot(Ca.Lon,Ca.Lat,'k')
        title('Modeled topography')
        
        
        % Interpolate Ca coordinates to km scale
        CaLonKm = griddata(qx,qy,qx_km,Ca.Lon,Ca.Lat);
        CaLatKm = griddata(qx,qy,qy_km,Ca.Lon,Ca.Lat);
        CaTopKm = griddata(X(:,:,end),Y(:,:,end),Z_surf,CaLonKm, CaLatKm);
        
        
        
        
        
        % Transform into Lat-lon coordinates
        X1          = [(X + 500)./1000]*12 + 235;
        Y1          = [(Y + 750)./1500]*13 + 29;
        Z1          = Z;
        Z1(:,:,end) = Z_surf/1e3/100;
        
        col = 'k';
        figure(4),clf, hold on
        plot3(squeeze(X(:,1,:)),  squeeze(Y(:,1,:)),  squeeze(Z(:,1,:))  ,col)
        plot3(squeeze(X(:,1,:)).',squeeze(Y(:,1,:)).',squeeze(Z(:,1,:)).',col)
        plot3(squeeze(X(1,:,:)),  squeeze(Y(1,:,:)),  squeeze(Z(1,:,:))  ,col)
        plot3(squeeze(X(1,:,:)).',squeeze(Y(1,:,:)).',squeeze(Z(1,:,:)).',col)
        plot3(squeeze(X(:,:,end)),  squeeze(Y(:,:,end)),  squeeze(Z(:,:,end))  ,col)
        plot3(squeeze(X(:,:,end)).',squeeze(Y(:,:,end)).',squeeze(Z(:,:,end)).',col)
        surf(X(:,:,end),Y(:,:,end),Z(:,:,end), Z_surf)
        
        caxis([-4500 4000])
        plot3(CaLonKm,CaLatKm,CaTopKm/1e3+10,'w','linewidth',2)
        
        
        %colorbar
        view(3), axis equal, axis tight
        
        % velocity
        stepx = 4; stepz=4;
        X_sl = squeeze(X(1:stepx:end,1,1:stepz:end));
        Y_sl = squeeze(Y(1:stepx:end,1,1:stepz:end));
        Z_sl = squeeze(Z(1:stepx:end,1,1:stepz:end));
        
        Vx_sl  = squeeze(Vx_3D(1:stepx:end,1,1:stepz:end));
        Vy_sl  = squeeze(Vy_3D(1:stepx:end,1,1:stepz:end));
        Vz_sl  = squeeze(Vz_3D(1:stepx:end,1,1:stepz:end));
        h = quiver3(X_sl(:),Y_sl(:),Z_sl(:),Vx_sl(:),Vy_sl(:),Vz_sl(:),'m')
        set(h,'AutoScaleFactor',1,'linewidth',2);
        
        
        stepx = 4; stepz=4;
        X_sl = squeeze(X(1,1:stepx:end,1:stepz:end));
        Y_sl = squeeze(Y(1,1:stepx:end,1:stepz:end));
        Z_sl = squeeze(Z(1,1:stepx:end,1:stepz:end));
        
        Vx_sl  = squeeze(Vx_3D(1,1:stepx:end,1:stepz:end));
        Vy_sl  = squeeze(Vy_3D(1,1:stepx:end,1:stepz:end));
        Vz_sl  = squeeze(Vz_3D(1,1:stepx:end,1:stepz:end));
        h = quiver3(X_sl(:),Y_sl(:),Z_sl(:),Vx_sl(:),Vy_sl(:),Vz_sl(:),'m')
        set(h,'AutoScaleFactor',1,'linewidth',2);
        
        
        title(['Time = ',num2str(round(time/characteristic.SecYear/1e3)*1e3),' years, max vertical velocity =',num2str(max(Vz_3D(:))*characteristic.cmYear,'%2.2f'),' cm/year'])
        
        xlabel('Width [km]')
        ylabel('Length [km]')
        zlabel('Depth [km]')
        
    end
    
    
    if (1==1)
        Matlab2VTK_intp;         % write VTK files
        Matlab2VTK;
    end
    
    
    if (1==0)
        figure(234), clf
        % plot subduction zone
        pcolor(squeeze(intpx_3D(1,:,:))/characteristic.km,squeeze(intpz_3D(1,:,:))/characteristic.km,log10(squeeze(mu_3D(1,:,:)))),
        %pcolor(squeeze(intpx_3D(1,:,:))/characteristic.km,squeeze(intpz_3D(1,:,:))/characteristic.km,squeeze(Phases_3D(1,:,:)))
        
        hold on
        
        shading flat,
        h=colorbar
        title(h,['log_{10}(\mu)'])
        
        %title(['Time =',num2str(time/characteristic.SecYear/characteristic.Myrs,'%1.1f'),' Myrs, 32x128x2 nodes'])
        title(['Time =',num2str(time/characteristic.SecYear/characteristic.Myrs,'%1.1f'),' Myrs, 64x256x2 nodes'])
        
        axis equal, axis tight
        
        figure(234)
        set(gcf,'PaperPositionMode','auto')
        set(gcf, 'Renderer', 'painter')
        set(get(gcf,'JavaFrame'),'Maximized',1);
        print('-djpeg90','-r300',[SimName,'_',num2str(time/characteristic.SecYear/1e6,'%3.0f'),'Myrs_phases.jpg'])
        
        
    end
    
    if visualization
        figure(12223)
        surf(X(:,:,end)/1e3, Y(:,:,end)/1e3, Z(:,:,end))
        
        title(['Width of domain = ',num2str(round((max(X(:,1,end))-min(X(:,1,end)))/1e3)),'km'])
        
    end
    
    
    

%     pause
    %            error('stop')
end

if plast == 1
    [EXX, ix] = sort(EXX);
    TXX = TXX(ix);
    figure(1)
    plot(abs(EXX),abs(TXX),'x')
end


