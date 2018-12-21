function [Location] = TracingParticles_LaMEM(Velocity, Location, time, time_end_Myrs, TimeIntegrationSteps, Mode);
% This routine uses LaMEM matlab files with velocities to trace a particle forwards/backwards
% in time. The velocity field is assumed to be fixed with time.
% We assume GEO-units, with lengths in [km] and velocities in [cm/yr].
%
% A 4th order Runga-Kutta method is used for advecting.
%
% Usage:
%
%  [Location] = TracingParticles(Velocity, StartingLocation, time_end_Myrs, TimeIntegrationSteps, Mode);;
%
%       Input:
%               Directory               -   Directory which has the MVEP2 output files
%               starting_filenumber     -   Starting file #
%               StartingLocation        -   Structure with the following fields:
%                                           StartingLocation.x - in [m]
%                                           StartingLocation.z - in [m]
%               Step                    -   Spacing between the Outputs
%                                           (Attention the higher the
%                                           spacing the less exact the
%                                           tracing)
%               Mode                    -   'Forward' or 'Backwards' in time 
%
%

% TODO:
% [Note: a future version should take time-varying velocity fields into
% account]


% % Directory       = 'CrustalConvectionMultiplePulses_HiRes';
% % starting_filenumber = 10;
% Directory = pwd;
% starting_filenumber = 13;

% % Indicate the starting location in m's - you can indicate several starting points
% StartingLocation_x  = 1900e3;
% StartingLocation_z  = -100e3;

% start_BackwardTracing;

switch Mode
    case 'Forward'
        TimeStepDirection =  1;
   
    case 'Backwards'
        TimeStepDirection = -1;
       
    otherwise
        error('Choose whether you advect Forward or Backward in time') 
end
            

numTimesteps 	= length(Velocity);

% Error checking - if loaded directly from LaMEM, the dimensions of the
% arrays are likely inconsistent with MESHGRID, and the interpolation
% routines will give trouble
x_vec           = Velocity{1}.x(1,:,1);
if (max(x_vec)-min(x_vec))<1e-3
    Shift = true;
    warning('provided arrays inconsistent with MESHGRID; shifting them')
    
    for i=1:numTimesteps
        % check based on
        names           = fieldnames(Velocity{i});
        
        for ii=1:length(names)
            data        =   getfield(Velocity{i},names{ii});
            data        =   shiftdim(data,1);
            Velocity{i} =   setfield(Velocity{i},names{ii}, data);
        end
        
    end
else
    Shift = false;
end

for i=1:numTimesteps
    Velocity{i}.x_vec = squeeze(Velocity{i}.x(1,:,1));
    Velocity{i}.y_vec = squeeze(Velocity{i}.y(:,1,1));
    Velocity{i}.z_vec = squeeze(Velocity{i}.z(1,1,:));
end


% Filename of current timestep
BoundingBoxSize     =   100;          % in km


% Create row vectors out of it
Location.x          =   Location.x(:);     % dimensionless
Location.y          =   Location.y(:);     % dimensionless
Location.z          =   Location.z(:);     % dimensionless
Location.time       =   time;


% Only check
BoundingBoxSize     =   BoundingBoxSize;

min_x               =   min(Location.x(:,end)) - BoundingBoxSize;
max_x               =   max(Location.x(:,end)) + BoundingBoxSize;
min_y               =   min(Location.y(:,end)) - BoundingBoxSize;
max_y               =   max(Location.y(:,end)) + BoundingBoxSize;
min_z               =   min(Location.z(:,end)) - BoundingBoxSize;
max_z               =   max(Location.z(:,end)) + BoundingBoxSize;

% coordinate grids
X                   =   Velocity{1}.x;
Y                   =   Velocity{1}.y;
Z                   =   Velocity{1}.z;

ind_nearby_nodes    =   find( (X    > min_x)  & (X < max_x ) & (Y    > min_y)  & (Y < max_y ) & (Z > min_z)  & (Z < max_z ) );
while length(ind_nearby_nodes)<10
    display('Bounding Box is too small .. Made it larger to find nearby nodes and integration points!')
    BoundingBoxSize     =   2*BoundingBoxSize;
    min_x               =   min(Location.x(:,end))-BoundingBoxSize;
    max_x               =   max(Location.x(:,end))+BoundingBoxSize;
    min_y               =   min(Location.y(:,end)) - BoundingBoxSize;
    max_y               =   max(Location.y(:,end)) + BoundingBoxSize;
    min_z               =   min(Location.z(:,end))-BoundingBoxSize;
    max_z               =   max(Location.z(:,end)+BoundingBoxSize);
    
    ind_nearby_nodes    =   find( (X    > min_x)  & (X < max_x ) & (Y    > min_y)  & (Y < max_y ) & (Z > min_z)  & (Z < max_z ) );
end

% Determine if we have one timestep as onput (=steady state) or multiple
% ones
if numTimesteps==1
    SteadyState = true;
else
    SteadyState = false;
end


% Initialize
Location.Vx_cmYr    =   [];
Location.Vy_cmYr    =   [];
Location.Vz_cmYr    =   [];
Location.time_Myrs  =   0;

%% Compute Particle Location and Parameters for every timestep without the first one

% Define number of steps we want to take
% IntegrationSteps = 10;       % default value
dt               = (abs(time_end_Myrs)-time)/TimeIntegrationSteps*1e6;       % timestep in yrs


for num =1:TimeIntegrationSteps
    
    
    if ~SteadyState
        % in case we have multiple timesteps we need to find out where in
        % time we are
        error('to be done')
        id_current = 1;
        id_old     = 1;
        
    else
        id_current = 1;
        id_old     = 1;
    end
    
    % Retrieve data of current step
    
    
    % Compute relevant points
    min_x               =   min(Location.x(:,end))-BoundingBoxSize;
    max_x               =   max(Location.x(:,end))+BoundingBoxSize;
    min_y               =   min(Location.y(:,end)) - BoundingBoxSize;
    max_y               =   max(Location.y(:,end)) + BoundingBoxSize;
    min_z               =   min(Location.z(:,end))-BoundingBoxSize;
    max_z               =   max(Location.z(:,end)+BoundingBoxSize);
    
    % to be modified to take 3 index arrays into account
    ix                  =   find( (Velocity{id_current}.x_vec>min_x) &  (Velocity{id_current}.x_vec<max_x));
    iy                  =   find( (Velocity{id_current}.y_vec>min_y) &  (Velocity{id_current}.y_vec<max_y));
    iz                  =   find( (Velocity{id_current}.z_vec>min_z) &  (Velocity{id_current}.z_vec<max_z));
    
%     ind_nearby_nodes    =   find( (X    > min_x)  & (X < max_x ) & (Y    > min_y)  & (Y < max_y ) & (Z > min_z)  & (Z < max_z ) );
    
    % Compute grids and velocity at dt/2 (use linear interpolation in time)
    x_half              =   ( Velocity{id_old}.x  + Velocity{id_current}.x )/2;
    y_half              =   ( Velocity{id_old}.y  + Velocity{id_current}.y )/2;
    z_half              =   ( Velocity{id_old}.z  + Velocity{id_current}.z )/2;
    
    Vx_half             =   ( Velocity{id_old}.Vx + Velocity{id_current}.Vx )/2;
    Vy_half             =   ( Velocity{id_old}.Vy + Velocity{id_current}.Vy )/2;
    Vz_half             =   ( Velocity{id_old}.Vz + Velocity{id_current}.Vz )/2;
    
    % We follow the wikipedia nomenclature here for RK4, with the
    % difference that h=-dt  (see
    % http://en.wikipedia.org/wiki/Runge?Kutta_methods) plus that we do
    % it in 2D:
    %
    % y_n+1 = y_n + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    % t_n+1 = t_n - dt;
    %
    % where
    %     k1 = -dt*f(t_n       , y_n         )
    %     k2 = -dt*f(t_n-1/2*dt, y_n + 1/2*k1)
    %     k3 = -dt*f(t_n-1/2*dt, y_n + 1/2*k2)
    %     k4 = -dt*f(t_n-    dt, y_n +     k3)
    
         
    Vx1     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vx_half(iy,ix,iz),Location.x(:,end),Location.y(:,end),Location.z(:,end));
    Vy1     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vy_half(iy,ix,iz),Location.x(:,end),Location.y(:,end),Location.z(:,end));
    Vz1     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vz_half(iy,ix,iz),Location.x(:,end),Location.y(:,end),Location.z(:,end));
    
    k1_x    = TimeStepDirection*dt*Vx1/1e5;     % transform from cm/yr to km/yr
    k1_y    = TimeStepDirection*dt*Vy1/1e5;
    k1_z    = TimeStepDirection*dt*Vz1/1e5;
    
    Vx2     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vx_half(iy,ix,iz),Location.x(:,end)  + 0.5*k1_x,Location.y(:,end)  + 0.5*k1_y,Location.z(:,end)  + 0.5*k1_z);
    Vy2     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vy_half(iy,ix,iz),Location.x(:,end)  + 0.5*k1_x,Location.y(:,end)  + 0.5*k1_y,Location.z(:,end)  + 0.5*k1_z);
    Vz2     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vz_half(iy,ix,iz),Location.x(:,end)  + 0.5*k1_x,Location.y(:,end)  + 0.5*k1_y,Location.z(:,end)  + 0.5*k1_z);
    k2_x    = TimeStepDirection*dt*Vx2/1e5;
    k2_y    = TimeStepDirection*dt*Vy2/1e5;
    k2_z    = TimeStepDirection*dt*Vz2/1e5;
    
   
    Vx3     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vx_half(iy,ix,iz),Location.x(:,end)  + 0.5*k2_x,Location.y(:,end)  + 0.5*k2_y,Location.z(:,end)  + 0.5*k2_z);
    Vy3     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vy_half(iy,ix,iz),Location.x(:,end)  + 0.5*k2_x,Location.y(:,end)  + 0.5*k2_y,Location.z(:,end)  + 0.5*k2_z);
    Vz3     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vz_half(iy,ix,iz),Location.x(:,end)  + 0.5*k2_x,Location.y(:,end)  + 0.5*k2_y,Location.z(:,end)  + 0.5*k2_z);
    k3_x    = TimeStepDirection*dt*Vx3/1e5;
    k3_y    = TimeStepDirection*dt*Vy3/1e5;
    k3_z    = TimeStepDirection*dt*Vz3/1e5;
    
    
    Vx4     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vx_half(iy,ix,iz),Location.x(:,end)  + 0.5*k3_x,Location.y(:,end)  + 0.5*k3_y,Location.z(:,end)  + 0.5*k3_z);
    Vy4     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vy_half(iy,ix,iz),Location.x(:,end)  + 0.5*k3_x,Location.y(:,end)  + 0.5*k3_y,Location.z(:,end)  + 0.5*k3_z);
    Vz4     = interp3(x_half(iy,ix,iz),y_half(iy,ix,iz),z_half(iy,ix,iz),Vz_half(iy,ix,iz),Location.x(:,end)  + 0.5*k3_x,Location.y(:,end)  + 0.5*k3_y,Location.z(:,end)  + 0.5*k3_z);
    
    k4_x    = TimeStepDirection*dt*Vx4/1e5;
    k4_y    = TimeStepDirection*dt*Vy4/1e5;
    k4_z    = TimeStepDirection*dt*Vz4/1e5;
    
    % Compute location of the particle @ OldStep
    Location.x(:,end+1) = Location.x(:,end) + 1/6*(k1_x + 2*k2_x + 2*k3_x + k4_x);
    Location.y(:,end+1) = Location.y(:,end) + 1/6*(k1_y + 2*k2_y + 2*k3_y + k4_y);
    Location.z(:,end+1) = Location.z(:,end) + 1/6*(k1_z + 2*k2_z + 2*k3_z + k4_z);
    
    % Store velocity on the particles
    Location.Vx_cmYr(:,end+1) = (Location.x(:,end)-Location.x(:,end-1))/dt*1e5;      % average velocity
    Location.Vy_cmYr(:,end+1) = (Location.y(:,end)-Location.y(:,end-1))/dt*1e5;      % average velocity
    Location.Vz_cmYr(:,end+1) = (Location.z(:,end)-Location.z(:,end-1))/dt*1e5;      % average velocity in cm/yr
    
    % Compute properties at the particle (temperature, etc)
    time = time + TimeStepDirection*dt;
    Location.time_Myrs(end+1)= time/1e6;
    
   
end






% time=CurrentStep.time;

% Transform into dimensional units before outputting
% SecYear = 3600*24*365.25;
% Location.time_years             =   Location.time*CHAR.Time/SecYear;                    % in  years
% Location.x                      =   Location.x*CHAR.Length;                             % in m
% Location.z                      =   Location.z*CHAR.Length;                             % in m
% Location.Temperature_Celcius    =   Location.Temperature*CHAR.Temperature-273.15;       % in C
% 
% Location = rmfield(Location,'Temperature');     % to avoid confusion in which units this is
% Location = rmfield(Location,'time_years');      % to avoid confusion in which units this is
% 










