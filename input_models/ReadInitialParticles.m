function [num_proc, Part, Partic, info] = ReadInitialParticles(rank,Directory)
% ReadInitialParticles
%
% This routine reads the initial particle information that's created by a LaMEM simulation
% 
% Note: The directory in which the initial particles are located is
% specified, as well as the rank
%
%
% $Id$ Boris Kaus

if nargin==0
    Directory   =   'InitialParticles';
    rank        =   0;
elseif nargin==1
    Directory   =   'InitialParticles';
end

addpath('./LaMEM_MATLAB');  % add I/O routines


[info]   =   PetscBinaryRead(['./',Directory,'/Initial_Particles.0.out']);
num_proc =   info(1);

Part.x   = []; Part.y     = []; Part.z   = [];  Part.num = []; Part.phase = []; 

[info,Partic]    =   PetscBinaryRead(['./',Directory,'/Initial_Particles.',num2str(rank),'.out']);
num_particle_prop   =   info(3);
    
%
Part.x     =   [Part.x;     Partic( 1:num_particle_prop:end)];
Part.y     =   [Part.y;     Partic( 2:num_particle_prop:end)];
Part.z     =   [Part.z;     Partic( 3:num_particle_prop:end)];
Part.num   =   [Part.num;   Partic( 4:num_particle_prop:end)];
Part.phase =   [Part.phase; Partic( 5:num_particle_prop:end)];

