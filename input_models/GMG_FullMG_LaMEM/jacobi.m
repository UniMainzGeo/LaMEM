function Sol_vec = jacobi(Sol_vec, Rhs_vec, MG_Info, level,omega_1,omega_2,nu1);
% Perform a jacobi smoother step
%

% Retrieve data
DirichletNodes          =   MG_Info.BC{level}.DirichletNodes;
Free                    =   MG_Info.BC{level}.Free;
A                       =   MG_Info.Matrixes{level}.A;
A_precondition          =   MG_Info.Matrixes{level}.A_precondition;

InvDiag                 =   1./diag(A_precondition);       % this particular preconditioner uses the diagonal only
InvDiagfree             =   InvDiag(Free);
numV                    =   MG_Info.Numbering{level}.numV;
for itime=1:nu1

    % Compute the residual
    [Res, Res_free]                   =   Residual(A, Sol_vec, Rhs_vec, MG_Info, level);

%  For debugging: transfer residuals to 3D grids
%
%     Error_Vx =  Res(MG_Info.Numbering{level}.Number_Vx);
%     Error_Vy =  Res(MG_Info.Numbering{level}.Number_Vy);
%     Error_Vz =  Res(MG_Info.Numbering{level}.Number_Vz);
%     Error_P  =  Res(MG_Info.Numbering{level}.Number_P);
%     
%     Vx       =  Sol_vec(MG_Info.Numbering{level}.Number_Vx);
%     Vy       =  Sol_vec(MG_Info.Numbering{level}.Number_Vy);
%     Vz       =  Sol_vec(MG_Info.Numbering{level}.Number_Vz);
%     P        =  Sol_vec(MG_Info.Numbering{level}.Number_P);
%     
%     
% %     figure(1), plot(log10(abs(Res)),'o'), drawnow
% %     disp(['   Iter ',num2str(itime),' error=',num2str(norm(Res))])
     
    % Compute correction, which simply uses 
    dR                      =   zeros(size(Rhs_vec));
    dR(Free)                =   InvDiagfree.*Res_free;

   
    dR(DirichletNodes)      =   0;
    dR                      =   dR*omega_1;               	% As applied to velocity equations
    dR(numV+1:end)          =   dR(numV+1:end)*omega_2;     % As applied to pressure equations
    
    
    % Add correction
    Sol_vec                 =   Sol_vec + dR;
    
end
