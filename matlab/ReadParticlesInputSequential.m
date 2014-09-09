function A = ReadParticlesInputSequential(filename,options)
%
% options.flippedxy     = 0;
% options.W             = 1120e3;
% options.L             = 885e3;
% options.H             = 420e3;
% options.CL            = 1e6;
% options.npart_x       = 3;
% options.npart_y       = 3;
% options.npart_z       = 3;
% 
% A = ReadParticlesInputSequential('ParticlesInput3D.dat',options);
%
% Tobias Baumann, Mainz University, 2013
%
% =========================================================================
A         = options;
PhaseVec  = PetscBinaryRead(filename);

if(options.flippedxy==1)
    A.nump_x        = PhaseVec(2);
    A.nump_y        = PhaseVec(3);
else
    A.nump_x        = PhaseVec(3);
    A.nump_y        = PhaseVec(2);
end
A.nump_z = PhaseVec(1);
x_left   = 0;
y_front  = 0;
z_bot    = 0;
A.Phase  = PhaseVec(4:length(PhaseVec(4:end))/2+3);
A.Temp   = PhaseVec(length(PhaseVec(4:end))/2+4:end);
A.x      = x_left :A.W/(A.nump_x-1):x_left+A.W;
A.y      = y_front:A.L/(A.nump_y-1):y_front+A.L;
A.z      = z_bot  :A.H/(A.nump_z-1):z_bot+A.H;
A.Phase  = reshape(A.Phase,length(A.x),length(A.y),length(A.z));
if(options.flippedxy==1) 
    amat = reshape(A.Phase,length(A.y),length(A.x),length(A.z));
    num = 1;
    for k=1:A.nump_z
        for j=1:A.nump_y
            for i=A.nump_x:-1:1
                A.Phase(num)=amat(j,i,k);
                num=num+1;
            end 
        end
    end
end

end
