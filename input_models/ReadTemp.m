%ReadTemp
%
clear

time_step = 0;

% Read all the data, and reconstruct it into 3D arrays======================
[A, rhs]   =   PetscBinaryRead(['Temperature.dat']);
