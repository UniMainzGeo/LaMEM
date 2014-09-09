%ModifyInitialParticles_LithosphericChannels
%
% (1) Reads in the initial particles (for each processor separately)
% (2) Modifies the particle phases
% (3) Writes the particle information back to file (again for each processor)
%
%

clear, clf

addpath('../matlab/')

Directory = 'InitialParticles';

% Read in the initial particles on rank 0 (to find the total # of procs)
[num_proc] = ReadInitialParticles(0,Directory);


Width_Channel = 20;     % in km
RadiusDiapir  = 20;     % in km
ThickLithos   = 100;    % in km 


l_c = 100e3;

figure(1),clf
style = {'r.','b.','k.','y.','g.','m.','c.','yo','go'}
for rank=0:num_proc-1     % loop over all processors
    
    % Read particles
    [num_proc, Part, Partic, info] = ReadInitialParticles(rank,Directory);
    
    % Modify the particles (based on [x,y,z] coordinates of the particles)
    Part.x = Part.x*l_c/1e3; % in km
    Part.y = Part.y*l_c/1e3; % km
    Part.z = Part.z*l_c/1e3; % km
    
    ind = find(Part.z>-ThickLithos & Part.z<-40); Part.phase(ind) = 1;   % mantle lithosphere
    
    % add layering to upper & lower crust
    Thick_layer = 5;
    for depth=-200:2*Thick_layer:0
        
        % lower crust 1
        z_start = depth;
        z_end   = z_start + Thick_layer;
        ind     = find(Part.z>z_start & Part.z<z_end & Part.z>-40 & Part.z<-20);    % lower crust 1
        Part.phase(ind) = 2;
        
        
        % lower crust 2
        z_start = depth+ Thick_layer;
        z_end   = z_start + 2*Thick_layer;
        ind     = find(Part.z>z_start & Part.z<z_end & Part.z>-40 & Part.z<-20);    % lower crust 2
        Part.phase(ind) = 3;
        
          % upper crust 1
        z_start = depth;
        z_end   = z_start + Thick_layer;
        ind     = find(Part.z>z_start & Part.z<z_end & Part.z>-20 & Part.z<-0);    % upper crust 1 
        Part.phase(ind) = 4;
        
        
        % upper crust 2
        z_start = depth+ Thick_layer;
        z_end   = z_start + 2*Thick_layer;
        ind     = find(Part.z>z_start & Part.z<z_end & Part.z>-20 & Part.z<0);    % upper crust 2
        Part.phase(ind) = 5;
        
    end
    
    % sticky air
    ind = find(Part.z>0 ); Part.phase(ind) = 6;   %
    
    % Add channel through mantle lithosphere 
    ind = find(Part.phase==1 & abs(Part.x)<Width_Channel/2 ); Part.phase(ind) = 7;   % channel
    
    % add diapir
    
     ind = find( ((Part.x).^2 + (Part.z+(ThickLithos+RadiusDiapir-1)).^2)< RadiusDiapir.^2  );
      Part.phase(ind) = 8; 
     
    
    
    
    % Plot the particles
    num_phases = max(Part.phase);
    for iphase=0:num_phases
        ind= find(Part.phase==iphase);
        style_number = iphase+1;
        
        plot3(Part.x(ind), Part.y(ind), Part.z(ind), style{style_number});
        hold on
        
    end
    axis equal
    
    
    Part.x = Part.x*1e3/l_c; % in km
    Part.y = Part.y*1e3/l_c; % km
    Part.z = Part.z*1e3/l_c; % km
    
    % Save the info of that rank into a LaMEM-readable file
     WriteInitialParticles(rank,Part, Partic, info);
    
end



xlabel('Width [km]')
zlabel('Depth [km]')
view(0,0);