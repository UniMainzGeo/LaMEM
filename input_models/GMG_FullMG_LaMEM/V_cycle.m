function u = V_cycle(MG_Info,u, Rhs_vec,level,nu1, nu2, omega_1, omega_2, Method);    % Performs a V cycle; recursive.
% Performs a single V-cycle and employs a direct solve @ the coarsest level
%
% Boris Kaus, 2010


% Smooth
u     = jacobi(u, Rhs_vec, MG_Info, level, omega_1, omega_2, nu1);

if (level==1)                     % Coarsest grid
    % Retrieve data
%     C                       =   MG_Info.Matrixes{level}.C;
%     F                       =   MG_Info.Matrixes{level}.F;
%     Bt                      =   MG_Info.Matrixes{level}.Bt;
%     B                       =   MG_Info.Matrixes{level}.B;
%     A                       =   [F Bt; B 1e-10*eye(size(C))];
%     
    numV                    =   MG_Info.Numbering{1}.numV;
    numP                    =   MG_Info.Numbering{1}.numP;
    
    A                       =   MG_Info.Matrixes{level}.A;
   
    % We use a direct solver on the coarsest grid to solve the system:
    %  | VV VP | |v|   |f|
    %  | PV PP | |p| = |g|
    % 
    %  the PP block is zero, but for numerical reasons we put a small value
    %  on the diagonal:
    A(numV+1:end,numV+1:end) = 1e-10*eye(numP,numP);
    
    DirichletNodes          =   MG_Info.BC{level}.DirichletNodes;
    Free                    =   1:size(A,1);
    Free(DirichletNodes)    =   [];
    
    % Perform a direct solve
    u(Free)                 =   A(Free,Free)\Rhs_vec(Free);     
    
    return
end

% Compute residual
r   =   Residual(MG_Info.Matrixes{level}.A, u, Rhs_vec, MG_Info, level);


% Compute residual
fc      = Restrict(r,level, MG_Info, Method);                                   % Right-hand side of coarser grid
uc      = zeros(size(fc));                                                      % Zero initial guess on coarse grid
uc      = V_cycle(MG_Info,uc, fc,level-1,nu1, nu2, omega_1, omega_2, Method);   % Recursive call
u_add   = Prolongate(uc,level-1, MG_Info, Method);                              % Add correction from coarser grid to 



u       = u + u_add;                    % Interpolate and add correction

                  
u       = jacobi(u, Rhs_vec, MG_Info, level, omega_1, omega_2, nu2);  % Perform nu2 relaxation sweeps

return