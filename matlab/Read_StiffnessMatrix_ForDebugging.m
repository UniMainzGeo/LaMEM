%Read_StiffnessMatrix_ForDebugging
%
% reads the stiffness matrix from disk and solves it in matlab

clear, clf

addpath('./LaMEM_MATLAB')

% Read the files from disk:
disp(['Reading files'])
tic

% Read stiffness matrixes & rhs
[VV, VP, PV, PP, approx_S, InfoVec]   =  PetscBinaryRead('stiffness.dat');
[f,g, sol_vel, sol_P]                 =  PetscBinaryRead('rhs_vector.dat');

[T_MAT, InfoVec_T]                      =  PetscBinaryRead('stiffness_T.dat');
[rhs_T, T_sol]                          =  PetscBinaryRead('rhs_vector_T.dat');

%[f]                       =  PetscBinaryRead('rhs_vector.dat');

% g= zeros(length(PP),1);
disp(['Done in ',num2str(toc),'s'])

error('stop here')
%VP = PV';

incomp          = PV*sol_vel;
ForceBalance    = VP*sol_P + VV*sol_vel;    % force balance equations


error_incomp = max(abs( incomp ))
error_Force1 = max(abs( ForceBalance ))


%PV=VP';


% construct a 3D numbering matrix
nel_x                       =   InfoVec(1);
nel_y                       =   InfoVec(2);
nel_z                       =   InfoVec(3);
numVx = zeros(nel_x+1,nel_y+1,nel_z+1);
numVy = zeros(nel_x+1,nel_y+1,nel_z+1);
numVz = zeros(nel_x+1,nel_y+1,nel_z+1);
num=0;
for k=1:nel_z+1
    for j=1:nel_y+1
        for i=1:nel_x+1
            numVx(i,j,k) = num;
            num=num+1;
            
            numVy(i,j,k) = num;
            num=num+1;
            
            numVz(i,j,k) = num;
            num=num+1;
            
        end
    end
end



if 1==0
    % Solve with a direct solver
    A           	=   [VV VP; PV PP];
    Rhs             =   [f; g];
    Sol             =   A\Rhs;                  % Solve
    Sol_Vel        =   Sol(1:size(VV,1));
    Sol_P          =   Sol(size(VV,1)+1:end);
    disp(['solved with a direct method'])
    
end

if 1==0
    % solve with gmres
    tol = 1e-12;
    maxit = 4;
    restart = 30;
    
    if 1==1
        disp(['explictly computing Schur; better have sufficient memory!'])
        Vi = inv(VV);
        S = PV*Vi*VP;        % schur
        save SchurINPUT Vi VV VP PV PP approx_S S f g  Rhs
    end
    
    
    B = [VV VP;  PV*0 -S];          % preconditioning matrix
    [x,flag,relres,iter, resvec1] = gmres(A,Rhs,restart,tol,maxit,B);
    
    B = [VV VP;  PV*0 -approx_S];   % preconditioning matrix
    [x,flag,relres,iter, resvec2] = gmres(A,Rhs,restart,tol,maxit,B);
    
    B = [VV VP;  PV*0 -diag(diag(S))]; % preconditioning matrix
    [x,flag,relres,iter, resvec3] = gmres(A,Rhs,restart,tol,maxit,B);
    
    
    B = [VV VP;  PV*0 -PV*inv(diag(diag(VV)))*VP]; % preconditioning matrix
    [x,flag,relres,iter, resvec4] = gmres(A,Rhs,restart,tol,maxit,B);
    
    
    loglog(1:length(resvec1),resvec1,'-',...
        1:length(resvec2),resvec2,'-',...
        1:length(resvec3),resvec3,'-',...
        1:length(resvec4),resvec4,'-')
    
    
    
    
    legend('S=B*F^{-1}*B^T','1/eta*M_p','diag(S)','B*D^{-1}*B^T')
    xlabel('iteration')
    ylabel('residual')
    title('gmres for FC Stokes system; aspect ratio of elements ~1 (lines) ')
    
    
    
    
    
    
    %     Sol             =   x;                  % Solve
    %     Sol_Vel        =   Sol(1:size(VV,1));
    %     Sol_P          =   Sol(size(VV,1)+1:end);
end



if 1==0
    
    % Solve with a direct solver
    A           	=   [VV VP; PV PP*0];
    Rhs             =   [f; g];
    Sol             =   A\Rhs;                  % Solve
    Sol_Vel1        =   Sol(1:size(VV,1));
    Sol_P1          =   Sol(size(VV,1)+1:end);
    
    % Solve with powell hestenes iterations
    kappa           =   1e3;
    invPP           =   inv(PP);
    %invPP           =   spdiags(1./diag(PP),[0],size(PP,1),size(PP,1));
    VV_new          =   VV  +   kappa*VP*invPP*PV;
    max_div         =   realmax;
    P               =   zeros(size(g));
    while max_div>1e-6
        f_new       =   f-VP*P;
        Vel         =   VV_new\f_new;
        dP          =   kappa*invPP*PV*Vel;
        P           =   P + dP;
        
        Div         =   PV*Vel;
        max_div     =   max(abs(Div))
    end
    Sol_Vel         =   Vel;
    Sol_P           =   P;
    
    Sol_P(1:4:end)  =   Sol_P(1:4:end)- (Sol_P(1) - Sol_P1(1));     % make null-space equal for comparison
    
    
    Sol_Vel     = Sol_Vel1;
    
    error_vel = norm(Sol_Vel-Sol_Vel1)
    %     error_P   = norm(Sol_P  -Sol_P1  )
    
    
end

if 1==1
    %     Sol_Vel     = Sol_Vel1;
    
    
    
    %--------------------------------------------------------------------------
    %
    % Form element numbering etc.
    nel_x                       =   InfoVec(1);
    nel_y                       =   InfoVec(2);
    nel_z                       =   InfoVec(3);
    nnode_x                     =   InfoVec(4);
    nnode_y                     =   InfoVec(5);
    nnode_z                     =   InfoVec(6);
    
    node        =   1;
    dof_number  =   1;
    for iz=1:nnode_z;
                 for iy=1:nnode_y;
                     for ix=1:nnode_x;
          
     
                NodeNumber(ix,iy,iz)   =   node;
                
                for dof=1:3
                    DOF_Number(ix,iy,iz,dof)    =    dof_number;
                    dof_number                  =    dof_number+1;
                end
                
                node                   =   node+1;
            end
        end
    end
    
    % Set boundary conditions
    Vx_EQ           = squeeze(DOF_Number(:,:,:,1));
    Vy_EQ           = squeeze(DOF_Number(:,:,:,2));
    Vz_EQ           = squeeze(DOF_Number(:,:,:,3));
    %
    %--------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    Vx3D            =   zeros(size(Vx_EQ));
    Vy3D            =   zeros(size(Vx_EQ));
    Vz3D            =   zeros(size(Vx_EQ));
    
    Rho            =   zeros(size(Vx_EQ));
    for iz=1:size(Vx_EQ,3)
        Vx3D(:,:,iz) =   Sol_Vel(Vx_EQ(:,:,iz));
        Vy3D(:,:,iz) =   Sol_Vel(Vy_EQ(:,:,iz));
        Vz3D(:,:,iz) =   Sol_Vel(Vz_EQ(:,:,iz));
        
        Rho(:,:,iz) =   f(Vz_EQ(:,:,iz));
    end
    
    
    Vx2d    = squeeze(Vx3D(fix(end/2),:,:))';
    Vz2d    = squeeze(Vz3D(fix(end/2),:,:))';
    
    figure(1),clf
    contourf(Vx2d)
    hold on
    quiver(Vx2d,Vz2d)
    xlabel('x')
    ylabel('z')
    h=colorbar; title(h,'Vx')
    
    
    
    Vx2d    = squeeze(Vx3D(:,fix(end/2),:))';
    Vz2d    = squeeze(Vz3D(:,fix(end/2),:))';
    
    figure(2),clf
    contourf(Vx2d)
    hold on
    quiver(Vx2d,Vz2d)
    xlabel('x')
    ylabel('z')
    h=colorbar; title(h,'Vx')
    
    
end
