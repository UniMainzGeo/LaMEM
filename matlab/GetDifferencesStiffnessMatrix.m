function [diff] = GetDifferencesStiffnessMatrix(folder1,folder2)
%==========================================================================
%
%   compare stiffness matrixes  
%   
%   This tool was created in order to observe differences in the stiffness 
%   matrix of a certain timestep, for the initial run and a 2nd run restarted
%   from a breakpoint.
%
%   just specify the names of the two folders (full paths) where you've 
%   saved the data
%
%   expl. usage:
%   diff = GetDifferencesStiffnessMatrix('/Users/tobibaumann/01_LaMEM/bin/E01_Breakpointcheck/debug_00','/Users/tobibaumann/01_LaMEM/bin/E01_Breakpointcheck/debug_01_after')
%
%   (created:  August 2012, Tobias Baumann)
%
%==========================================================================

workingdir=pwd;



cd(folder1);
disp(['load' folder1]);
[VV, VP, PV, PP, approx_S, InfoVec]   =  PetscBinaryRead('stiffness.dat');
[f,g, sol_vel, sol_P]                 =  PetscBinaryRead('rhs_vector.dat');

% struct for timestep before saving breakpoint
before.InfoVec=InfoVec;
%before.InfoVec_T=InfoVec_T
before.PP=PP;
before.PV=PV;
before.VP=VP;
before.VV=VV;
before.approx_S=approx_S;
before.sol_P=sol_P;
before.sol_vel=sol_vel;
before.f=f;
before.g=g;

cd(folder2);
disp(['load' folder2]);
[VV, VP, PV, PP, approx_S, InfoVec]   =  PetscBinaryRead('stiffness.dat');
[f,g, sol_vel, sol_P]                 =  PetscBinaryRead('rhs_vector.dat');
% struct for timestep after saving the breakpoint
after.InfoVec=InfoVec;
%after.InfoVec_T=InfoVec_T
after.PP=PP;
after.PV=PV;
after.VP=VP;
after.VV=VV;
after.approx_S=approx_S;
after.sol_P=sol_P;
after.sol_vel=sol_vel;
after.f=f;
after.g=g;


% get the differnce
diff.InfoVec    =  after.InfoVec    - before.InfoVec ;
diff.PP         =  after.PP         - before.PP;
diff.PV         =  after.PV         - before.PV;
diff.VP         =  after.VP         - before.VP;
diff.VV         =  after.VV         - before.VV;
diff.approx_S   =  after.approx_S   - before.approx_S;
diff.sol_P      =  after.sol_P      - before.sol_P;
diff.sol_vel    =  after.sol_vel    - before.sol_vel;
diff.f          =  after.f          - before.f;
diff.g          =  after.g          - before.g;

diff.sol_vel_sum=sum(diff.sol_vel);
diff.sol_P_sum  =sum(diff.sol_P);
diff.f_sum      =sum(diff.f);
diff.g_sum      =sum(diff.g);


diff.PP_min_approxS = min(min(diff.approx_S));
diff.PP_max_approxS = max(max(diff.approx_S));

diff.PP_min     = min(min(diff.PP));
diff.PP_max     = max(max(diff.PP));

diff.PV_min     = min(min(diff.PV));
diff.PV_max     = max(max(diff.PV));

diff.VP_min     = min(min(diff.VP));
diff.VP_max     = max(max(diff.VP));

diff.VV_min     = min(min(diff.VV));
diff.VV_max     = max(max(diff.VV));

if((abs(diff.PP_max)+ ...
    abs(diff.VV_max)+ ...
    abs(diff.VP_max)+ ...
    abs(diff.PV_max)+ ...
    abs(diff.PP_max_approxS)+ ...
    abs(diff.sol_P_sum)+ ...
    abs(diff.sol_vel_sum)+ ...
    abs(diff.f_sum)+ ...
    abs(diff.g_sum)) == 0)
    disp('no differences observed');

else
    
    disp('differences observed!');
    disp(['Results save here' workingdir '/differences.mat']);
    
    save differences
    % plotting
    close all

    figure('name','VV');
    imagesc(diff.VV); colorbar; axis equal; title('diff VV');

    figure('name','PP');
    imagesc(diff.PP); colorbar; axis equal; title('diff PP');

    figure('name','PV');
    imagesc(diff.PV); colorbar; axis equal; title('diff PV');

    figure('name','VP');
    imagesc(diff.VP); colorbar; axis equal; title('diff VP');

%     savepdf(1,'VV_diff');
%     savepdf(2,'PP_diff');
%     savepdf(3,'PV_diff');
%     savepdf(4,'VP_diff');
end

cd(workingdir);
end