%ReadParticles
%
%Reads particles from disc - also in the case that those particles were
%generated and saved on more than one processor
clear, close all


addpath('./LaMEM_MATLAB');  % add I/O routines

for time_step = 0:1:10000;
    
    
    
    % read first one to learn how many files there should be
    [info]   =   PetscBinaryRead(['Particles.',num2str(0),'.',num2str(time_step+1000000),'.out']);
    num_proc =   info(2);
    
    Part.x   = []; Part.y     = []; Part.z   = [];  Part.num = []; Part.phase = []; Part.cpu = [];
    Part.ix  = []; Part.iy    = []; Part.iz   = []; Part.T   = [];
    
    for i=1:num_proc  % read all files
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
        
        
    end
    
    fname = 'Diapir_output';
    %fname = 'Subduction_output';
    
    %fname = 'Convection_output';
%     fname = 'Subduction3D_Output';
    fname = 'FallingBlock3D';
    
    %[info,coord_x,coord_y,coord_z, Vx, Vy, Vz, mu, rho, intpx, intpy, intpz,characteristic]  =   PetscBinaryRead([fname,'.',num2str(i-1),'.',num2str(time_step+1000000),'.out']);
    
    [info,coord_x,coord_y,coord_z, Vx, Vy, Vz, mu, rho, n,G,...
        C,phi,k,Cp,Q,alpha,FK,intpx, intpy, intpz, characteristic, temp]  =   PetscBinaryRead([fname,'.',num2str(i-1),'.',num2str(time_step+1000000),'.out']);
    
    
    xs = info(19);    ys = info(20);    zs = info(21);    xm = info(22);    ym = info(23);    zm = info(24);
    X  = reshape(coord_x,[xm ym zm]);       Y = reshape(coord_y,[xm ym zm]);        Z = reshape(coord_z,  [xm ym zm]);
    
    
    char.cmYear				 = characteristic(1);
    char.Myrs			     = characteristic(2);
    char.MPa                 = characteristic(3);
    char.SecYear             = characteristic(4);
    char.km                  = characteristic(5);
    char.ThermalExpansivity	 = characteristic(6);
    char.Strainrate			 = characteristic(7);
    char.kg					 = characteristic(8);
    char.Density			 = characteristic(9);
    char.Viscosity			 = characteristic(10);
    char.Temperature		 = characteristic(11);
    char.Velocity			 = characteristic(12);
    char.Stress				 = characteristic(13);
    char.Time				 = characteristic(14);
    char.Length				 = characteristic(15);
    
    Part.x = Part.x*char.Length;
    Part.y = Part.y*char.Length;
    Part.z = Part.z*char.Length;
    
    
    miny=0.0; maxy=10.1;
    
    
    
    figure(1),clf, hold on
    
    ind = find(Part.phase<0);
    
    Part.x(ind)     =[];
    Part.y(ind)     =[];
    Part.z(ind)     =[];
    Part.cpu(ind)   =[];
    Part.phase(ind) =[];
    Part.ix(ind)     =[];
    Part.iy(ind)     =[];
    Part.iz(ind)     =[];
    
    
    
    ind = find(Part.phase==0 );      % element found on cpu 0
    plot3(Part.x(ind)/1e3,Part.y(ind)/1e3,Part.z(ind)/1e3,'k.')
    
    ind = find(Part.phase==1 );      % element found on cpu1
    plot3(Part.x(ind)/1e3,Part.y(ind)/1e3,Part.z(ind)/1e3,'r.')
   
       ind = find(Part.phase==2 );      % element found on cpu1
     plot3(Part.x(ind)/1e3,Part.y(ind)/1e3,Part.z(ind)/1e3,'y.')
   
     axis equal, axis tight
    view(180,0)
    %
    %  ind = find(Part.cpu==3 );      % element found on cpu1
    %  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'r.')
    
    
    title(time_step)
    pause
end

%
% ind = find(Part.phase==0);      %deactivated particles
%  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'r+')

%  ind = find(Part.phase==1 );      %deactivated particles
%  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'k+')
%
%  ind = find(Part.phase==2 & Part.x<-4e5);      %deactivated particles
%  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'y+')
% %
%  ind = find(Part.num==2671);      %deactivated particles
%  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'bo','markersize',10)
%
%


if 1==1
    
    % ind = find(Part.ix==0 & Part.iy==0 & Part.iz==6);      %advanced search and not found
    % plot3(Part.x(ind),Part.y(ind),Part.z(ind),'r*')
    %
    plot3(squeeze(X(:,1,:)),  squeeze(Y(:,1,:)),  squeeze(Z(:,1,:))  ,'b-')
    plot3(squeeze(X(:,1,:)).',squeeze(Y(:,1,:)).',squeeze(Z(:,1,:)).','b-')
    plot3(squeeze(X(1,:,:)),  squeeze(Y(1,:,:)),  squeeze(Z(1,:,:))  ,'b-')
    plot3(squeeze(X(1,:,:)).',squeeze(Y(1,:,:)).',squeeze(Z(1,:,:)).','b-')
    
    plot3(squeeze(X(:,:,1)),  squeeze(Y(:,:,1)),  squeeze(Z(:,:,1))  ,'b-')
    plot3(squeeze(X(:,:,1)).',squeeze(Y(:,:,1)).',squeeze(Z(:,:,1)).','b-')
    
    plot3(squeeze(X(:,:,end)),  squeeze(Y(:,:,end)),  squeeze(Z(:,:,end))  ,'b-')
    plot3(squeeze(X(:,:,end)).',squeeze(Y(:,:,end)).',squeeze(Z(:,:,end)).','b-')
    
    
    view(0,0)
    
    
    
    %===================================================================
    NumBgPx = 100;
    NumBgPy = 100;
    NumBgPz = 1;
    dx      = 1000/(NumBgPx-1);
    dy      = 2500/(NumBgPy-1);
    dz      =  300/(NumBgPz-1);
    
    x_left  = -500;
    y_front = -750;
    z_bot   = -300;
    
    
    [x3d,y3d,z3d] = meshgrid([x_left+0*dx/2:dx:(NumBgPx-1)*dx+x_left],[y_front+0*dy/2:dy:(NumBgPy-1)*dy+y_front],[z_bot+0*dz/2:dz:(NumBgPz-1)*dz+z_bot]);
    
    T3d     =  z3d;T3d=T3d*0;
    Phase3D =  z3d; Phase3D=Phase3D*0;
    
    
    load HeatFlowAndCrustalThicknessData
    Crust_interpolated  = griddata(crust2x_km,crust2y_km, crust2, squeeze(x3d(:,:,end)),squeeze(y3d(:,:,end)),'nearest');
    q_surf_interpolated = griddata(qx_km,     qy_km,      q, squeeze(x3d(:,:,end)),squeeze(y3d(:,:,end)),'nearest');
    Topo_interpolated   = griddata(qx_km,     qy_km,      topol_z, squeeze(x3d(:,:,end)),squeeze(y3d(:,:,end)),'nearest');
    
    %Crust_interpolated = griddata(crust2x_km,crust2y_km,crust2, coord_x_3D(:,:,end)/1e3,coord_y_3D(:,:,end)/1e3,'nearest');
    
    Tmantle     = 1300;
    k           = 3;
    kappa       = 1e-6;
    
    
    % Set the heatflow profile (pointwise) and set phases
    for ix=1:size(z3d,1)
        for iy=1:size(z3d,2)
            
            % Compute heatflow and thermal structure-----------------
            SecYear      =  3600*24*365.25;
            q_surf       =  q_surf_interpolated(ix,iy)*1e-3;                                      % in W/m2
            ThermalAge   =  (k*Tmantle/q_surf)^2/pi/kappa;                      % thermal age [s]
            
            z_depth      =  z3d(ix,iy,:)*1e3;                                   % in  m
            T_halfspace  =  Tmantle*erf(-z_depth./(2*sqrt(kappa*ThermalAge)));  % halfspace cooling model
            
            T3d(ix,iy,:) =  T_halfspace;
            
            % Compute phases ----------------------------------------
            crustal_thickness   = -Crust_interpolated(ix,iy);               %  if we use crust 2
            %  crustal_thickness   = -crustEars_new(ix,iy);        %  for the crustEARS model
            
            
            ind                 = find(z3d(ix,iy,:)>crustal_thickness);
            Topography          = Topo_interpolated(ix,iy);
            if (Topography>-2000)
                Phase3D(ix,iy,ind)  = 1;        % Continental crust
            else
                Phase3D(ix,iy,ind)  = 2;        % Oceanic crust
            end
        end
        ix
    end
    
    %
    % ind = find(Phase3D==2);
    % x3d=x3d(ind);
    % y3d=y3d(ind);
    % z3d=z3d(ind);
    % Phase3D=Phase3D(ind);
    
    
    %  ind = find(Part.phase==1 & Part.x<470e3);      %deactivated particles
    %  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'y+')
    %
    %
    %   ind = find(Part.phase==1 & Part.y<720e3);      %deactivated particles
    %  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'y+')
    
    
    
    
    %
    %  ind = find(Part.phase==2 );      %deactivated particles
    %  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'y+')
    
    %
    % ind = find(Part.ix==0 & Part.iy==0  & Part.iz==24 & Part.phase==2);      %deactivated particles
    % plot3(Part.x(ind),Part.y(ind),Part.z(ind),'r+')
    %
    %
    %  ind = find(Part.phase==0 );      %deactivated particles
    %  plot3(Part.x(ind),Part.y(ind),Part.z(ind),'g+')
    %
    
    ind = find(Part.ix==0 & Part.iy==20  & Part.iz==24 & Part.phase==2);      %deactivated particles
    plot3(Part.x(ind),Part.y(ind),Part.z(ind),'y+')
    
    ind = find(Part.ix==0 & Part.iy==20  & Part.iz==24 & Part.phase==0);      %deactivated particles
    plot3(Part.x(ind),Part.y(ind),Part.z(ind),'r+')
    
    
    % ind = find(x3d*1e3<-4.4531e+005 & x3d*1e3>-4.9219e+005 & y3d*1e3<2.2852e+005 & y3d*1e3> 1.9336e+005 & Phase3D==2);
    % plot3(x3d(ind)*1e3,y3d(ind)*1e3,z3d(ind)*1e3,'y.')
    %
    % ind = find(x3d*1e3<-4.4531e+005 & x3d*1e3>-4.9219e+005 & y3d*1e3<2.2852e+005 & y3d*1e3> 1.9336e+005 & Phase3D==0);
    % plot3(x3d(ind)*1e3,y3d(ind)*1e3,z3d(ind)*1e3,'r.')
    
    view(-2*90,0)
    %===================================================================
    
    coord_x_3D=X;
    coord_y_3D=Y;
    coord_z_3D=Z;
    
end
