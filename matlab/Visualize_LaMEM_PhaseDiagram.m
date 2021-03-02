function [] = Visualize_LaMEM_PhaseDiagram(PD_name)
% This file plots a LaMEM phase diagram.
%
% Usage:
%
%   1) Make sure that the LaMEM/matlab is on the path
%       e.g. addpath('../../matlab')
%   2) Visualize a LaMEM phase diagram with
%       Visualize_LaMEM_PhaseDiagram('TestPD.in')
%
%       That requires you to be in the same directory as 'TestPD.in'



% Indicate name
% PD_name = 'TestPD.in';

sym         =   0;
H           =	dlmread(PD_name, ' ', [49 0 54 0]);
Tmin_LaMEM  =   H(1);
dT_LaMEM    =   H(2);
nt          =   H(3);
Pmin_LaMEM  =   H(4);
dP_LaMEM    =   H(5);
np          =   H(6);

M           =   dlmread(PD_name, ' ', 55, 0);

rho         =   reshape(M(:,3),np,nt);
rho_fluid   =   reshape(M(:,1),np,nt);
melt        =   reshape(M(:,2),np,nt);
T_C         =   reshape(M(:,4),np,nt)-273.15;       % (LaMEM: degree K), Here: Celcius
P           =   reshape(M(:,5),np,nt)./1e3;         % kbar  - LaMEM input must be bar!


% Clear strange values in the PD and give warnings
NaN_in_Diagram = false;
ind = find(rho == 0);
if ~isempty(ind)
    warning('Rho had some weird values that may result in NaN and crash the LaMEM calculation, consider changing your PD')
    rho(ind) = 3000;
end
ind = find(rho_fluid == 0);
if ~isempty(ind)
    warning('Rho_fluid had some 0 values, consider changing your PD')
    rho_fluid(ind) = 3000;
end
ind = find(T_C <= 0);
if ~isempty(ind)
    warning('T had some weird values, PD is symmetrized')
    sym = 1;
end
ind = find(P <= 0);
if ~isempty(ind)
    warning('P had some weird values, PD is symmetrized')
    sym = 1;
end

ind = find(isnan(P));
if ~isempty(ind)
    warning('P has NaN values, suggesting a problem with the Perple_X calculation at those points')
    
    NaN_in_Diagram = true;
end
ind = find(isnan(T_C));
if ~isempty(ind)
    warning('T has NaN values, suggesting a problem with the Perple_X calculation at those points')
    
    NaN_in_Diagram = true;
end

% if sym == 1
%     % symmetrize
%     rho         = rho(1:end-1,1:end-1);
%     rho_fluid   = rho_fluid(1:end-1,1:end-1);
%     melt        = melt(1:end-1,1:end-1);
%     T_C         = T_C(1:end-1,1:end-1);
%     P           = P(1:end-1,1:end-1);
% end

% Check that Pmin/dP, Tmin/dT are the same as indicated for the LaMEM
% diagram
dP      =   (P(2,2)      -   P(1,1))*1e3;
dT      =   (T_C(2,2)    -   T_C(1,1));
Pmin    =   P(1,1)*1e3;
Tmin    =   T_C(1,1)+273.15;

if abs(dP-dP_LaMEM)>1e-10
    warning(['difference in dP:  dP_LaMEM=',num2str(dP_LaMEM),'  dP=',num2str(dP)])
end
if abs(dT-dT_LaMEM)>1e-10
    warning(['difference in dT:  dT_LaMEM=',num2str(dT_LaMEM),'  dT=',num2str(dT)])
end

if abs(Pmin-Pmin_LaMEM)>1e-10
    warning(['difference in Pmin: Pmin_LaMEM=',num2str(Pmin_LaMEM),'  Pmin=',num2str(Pmin)])
end
if abs(Tmin-Tmin_LaMEM)>1e-10
    warning(['difference in Tmin: Tmin_LaMEM=',num2str(Tmin_LaMEM),'  Tmin=',num2str(Tmin)])
end





figure(1),clf
subplot(221)
pcolor(T_C,P,rho)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('density rock w/out melt!')

subplot(222)
pcolor(T_C,P,melt)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('Melt content')

subplot(223)
pcolor(T_C,P,rho_fluid)
ylabel('P [kbar]')
xlabel('T')
title('density of fluid')
shading interp
colorbar

subplot(224)
rho_average = (1-melt).*rho + melt.*rho_fluid;
pcolor(T_C,P,rho_average)
ylabel('P [kbar]')
xlabel('T')
title('Combined density of fluid + solid (should NOT be passed to LaMEM!)')
shading interp
colorbar





% Correct the diagram if there are NaN's in P or T
if NaN_in_Diagram 
    
    ind                 = [find(isnan(T_C)); find(isnan(P)); find(isnan(rho))]; % indices of values to be corrected
    ind_correct         = find(T_C==T_C);
    ind_correct(ind)    = [];
    [ix,iz]             = ind2sub(size(T_C), ind);
    
    for i=1:length(ind) % do a few times, in case the neighbor is NaN too
        % Correct T and P values
        for i=1:length(iz)
            if ix(i)<size(P,1)
                P(ix(i),iz(i)) = P(ix(i)+1,iz(i));
            else
                P(ix(i),iz(i)) = P(ix(i)-1,iz(i));
            end
            
            if iz(i)<size(P,2)
                T_C(ix(i),iz(i)) = T_C(ix(i),iz(i)+1);
            else
                T_C(ix(i),iz(i)) = T_C(ix(i),iz(i)-1);
            end
        end
        
    end
    
    % at this stage P & T_C are correct
    rho(ind)       = griddata(P(ind_correct), T_C(ind_correct), rho(ind_correct),          P(ind), T_C(ind));
    rho_fluid(ind) = griddata(P(ind_correct), T_C(ind_correct), rho_fluid(ind_correct),    P(ind), T_C(ind));
    melt(ind)      = griddata(P(ind_correct), T_C(ind_correct), melt(ind_correct),         P(ind), T_C(ind));
    
    M1      =   [rho_fluid(:), melt(:), rho(:), T_C(:)+273.15 P(:)*1e3];        % output matrix    
   
    
    % generate a new LaMEM phase diagram that is corrected 
    PD_name_new = [PD_name(1:end-3),'_corrected.in'];
    fid         =   fopen(PD_name,'r');
    fid_new     =   fopen(PD_name_new,'w');
    for i=1:49
        line = fgetl(fid);
        fprintf(fid_new,'%s \n',line);
    end
    fprintf(fid_new,'%.6f \n',min(T_C(:))+273.15 );
    fprintf(fid_new,'%.6f \n',dT);
    fprintf(fid_new,'%i \n',nt);
    fprintf(fid_new,'%.6f \n',min(P(:))*1e3);
    fprintf(fid_new,'%.6f \n',dP);
    fprintf(fid_new,'%i \n',np);
    
    for i=1:size(M1,1)
         fprintf(fid_new,'%.5f %.5f %.5f %.5f %.5f \n',M1(i,:));
    end
    
    
    fclose(fid);
    fclose(fid_new);
    
    
    disp(['I made an attempt to autocorrect the phase diagram and saved it under ',PD_name_new])
end
    
