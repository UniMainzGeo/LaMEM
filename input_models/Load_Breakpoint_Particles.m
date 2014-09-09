%Load_Breakpoint_particles
%
%Reads particles from disc that are stored in a breakpoint file
clear, close all


addpath('./LaMEM_MATLAB');  % add I/O routines

Breakpoint_number = 0;
num_proc          = 16;
num_particle_prop = 28;



[coord_x_3D,coord_y_3D,coord_z_3D,info] = ReadInitialMesh('InitialMesh','MeshAfterRemeshing');


% read first one to learn how many files there should be
Part.x   = []; Part.y     = []; Part.z   = [];  Part.num = []; Part.phase = []; Part.cpu = [];
Part.ix  = []; Part.iy    = []; Part.iz  = []; Part.T   = [];
Part.eta = []; Part.zeta  = []; Part.phi = []; 

for i=0:num_proc-1  % read all Particle files
    [Partic]   =   PetscBinaryRead(['BreakPoint_Particles_',num2str(Breakpoint_number),'.',num2str(i),'.breakpoint']);
    
    Part.x     =   [Part.x;     Partic( 1:num_particle_prop:end)];
    Part.y     =   [Part.y;     Partic( 2:num_particle_prop:end)];
    Part.z     =   [Part.z;     Partic( 3:num_particle_prop:end)];
    Part.num   =   [Part.num;   Partic( 4:num_particle_prop:end)];
    Part.phase =   [Part.phase; Partic( 5:num_particle_prop:end)];
    Part.cpu   =   [Part.cpu;   Partic( 6:num_particle_prop:end)];
    Part.ix    =   [Part.ix;    Partic( 7:num_particle_prop:end)];
    Part.iy    =   [Part.iy;    Partic( 8:num_particle_prop:end)];
    Part.iz    =   [Part.iz;    Partic( 9:num_particle_prop:end)];
    
    Part.eta    =   [Part.eta;    Partic( 10:num_particle_prop:end)];
    Part.zeta    =   [Part.zeta;    Partic( 11:num_particle_prop:end)];
    Part.phi    =   [Part.phi;    Partic( 12:num_particle_prop:end)];
    
    
    Part.T     =   [Part.T;     Partic(13:num_particle_prop:end)];
end


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


% Plot the particles
style = {'r.','b.','k.','y.','g.','m.','c.','bo'}
num_phases = max(Part.phase);
for iphase=0:num_phases
    ind= find(Part.phase==iphase);
    style_number = iphase+1;
    
    plot3(Part.x(ind), Part.y(ind), Part.z(ind), style{style_number});
    hold on
    
end
axis equal

axis equal, axis tight
view(45,45)





