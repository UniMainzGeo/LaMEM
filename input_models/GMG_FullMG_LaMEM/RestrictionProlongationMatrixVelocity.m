% function [Rp,Pp] = RestrictionProlongationMatrixVelocity(level_fine, MG_Info);
%
% It interpolates the Vx values, defined at Vx-points, from a
% finer to a coarser mesh, using a technique recommended by a number of
% publications that is some kind of weighted linear interpolation.
%
% In addition, we need to take care that we have additional 'dummy' planes
% in the velocity matrixes

% debugging
level_fine      = 2;


level_coarse    = level_fine-1;


NumVx_coarse     =  MG_Info.Numbering{level_coarse}.Number_Vx;
NumVx_fine       =  MG_Info.Numbering{level_fine}.Number_Vx;

Rp              = sparse(max(NumVx_coarse(:)),max(NumVx_fine(:)));

for iy=2:size(NumVx_coarse,1)
    for ix=2:size(NumVx_coarse,2)
        for iz=1:size(NumVx_coarse,3)
            
            % Vx  is averaged from 12 surrounding cells
            ind_fine    = [];               ind_val = [];
            ind_fine    = [ind_fine, ];     ind_val = [ind_val, ];
            
            
            
%             ind_coarse  = NumP_coarse(iy,ix,iz);
            
%             Rp(ind_coarse, ind_fine(:)) = 1/8;
            
        end
    end
end

% Prolongation matrix is simply the transpose of the restriction one multiplied by a factor:
% Pp                      =   Rp'*8;



% % Testing the restriction operatator
% P_fine                  = (MG_Info.Grids{level_fine}.ZP);        % a P-value at a fine grid
% 
% % Set P-solution in solution vector
% Sol_vec                 =   zeros(max(NumP_fine(:)),1);
% Sol_vec(NumP_fine(:))   =   P_fine(:);
% 
% 
% % Restrict
% sol_vec_coarse          =   Rp*Sol_vec;
% 
% 
% % extract 3D matrix
% P_coarse                =   sol_vec_coarse(MG_Info.Numbering{level_coarse}.Number_P);
% 
% 
% 
% % Prolongate 
% sol_vec_coarse          =   ones(size(sol_vec_coarse));

% 
% sol_fine_prolongated    =   Pp'*sol_vec_coarse;
% P_fine_prolongated      =   sol_fine_prolongated(MG_Info.Numbering{level_fine}.Number_P);
% 



