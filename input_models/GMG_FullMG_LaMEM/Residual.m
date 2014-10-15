function [Res, Res_free] = Residual(A, Sol_vec, Rhs_vec, MG_Info, level);
% Computes the residual
%

% Retrieve data
Free                    =   MG_Info.BC{level}.Free;
DirichletNodes          =   MG_Info.BC{level}.DirichletNodes;

Res_free                =   Rhs_vec(Free) - A(Free,Free)*Sol_vec(Free);                         % residual of stokes equations


Res                     =   zeros(size(Sol_vec));
Res(Free)               =   Res_free;