function [ ] = WriteInitialParticles(rank, Part, Partic, info)
% ReadInitialParticles
%
% This routine reads the initial particle information that's created by a LaMEM simulation
% 
% Note: The directory in which the initial particles are located is
% specified, as well as the rank
%
%
% $Id$ Boris Kaus

Directory   =   'InitialParticles';
mkdir(Directory)
addpath('./LaMEM_MATLAB');  % add I/O routines

%
num_particle_prop                   =   info(3);
Partic( 1:num_particle_prop:end)    =   Part.x;
Partic( 2:num_particle_prop:end)    =   Part.y;
Partic( 3:num_particle_prop:end)    =   Part.z;
Partic( 4:num_particle_prop:end)    =   Part.num;
Partic( 5:num_particle_prop:end)    =   Part.phase;

fname = ['./',Directory,'/Initial_Particles.',num2str(rank),'.out'];

% Save mesh ---------------------------------------------------------------
PetscBinaryWrite(fname,info, Partic );
%--------------------------------------------------------------------------

