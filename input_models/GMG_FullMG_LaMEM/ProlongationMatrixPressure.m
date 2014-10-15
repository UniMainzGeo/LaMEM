function [P]  	=  ProlongationMatrixPressure(level_fine, MG_Info);

level_coarse = level_fine-1;
%% Pressure
% We use simple injection for this
NumP_coarse     =  MG_Info.Numbering{level_coarse}.Number_P;
NumP_fine       =  MG_Info.Numbering{level_fine}.Number_P;

P               = sparse(max(NumP_coarse(:)),max(NumP_fine(:)));

for iy=1:size(NumP_coarse,1)
    for ix=1:size(NumP_coarse,2)
        for iz=1:size(NumP_coarse,3)
            
            % Pressure is averaged from the 8 surrounding cells
            ind_fine    = NumP_fine(2*iy-1:2*iy,2*ix-1:2*ix,2*iz-1:2*iz);   
            ind_coarse  = NumP_coarse(iy,ix,iz);
            
            P(ind_coarse, ind_fine(:)) = 1;
            
        end
    end
end

P  = P';


