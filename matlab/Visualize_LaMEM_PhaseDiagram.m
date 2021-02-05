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
np          =   H(6);
nt          =   H(3);
Tmin_LaMEM  =   H(1);
dT_LaMEM    =   H(2);
Pmin_LaMEM  =   H(4);
dP_LaMEM    =   H(5);

M           =   dlmread(PD_name, ' ', 55, 0);

rho         =   reshape(M(:,3),np,nt);
rho_fluid   =   reshape(M(:,1),np,nt);
melt        =   reshape(M(:,2),np,nt);
T_C         =   reshape(M(:,4),np,nt)-273.15;       % (LaMEM: degree K), Here: Celcius
P           =   reshape(M(:,5),np,nt)./1e3;         % kbar  - LaMEM input must be bar!


% Clear strange values in the PD and give warnings
ind = find(rho == 0);
if ~isempty(ind)
    warning('Rho had some 0 values, consider changing your PD')
    rho(ind) = 3000;
end
ind = find(rho_fluid == 0);
if ~isempty(ind)
    warning('Rho_fluid had some 0 values, consider changing your PD')
    rho_fluid(ind) = 3000;
end
ind = find(T_C == 0);
if ~isempty(ind)
    warning('T had some weird values, PD is symmetrized')
    sym = 1;
end
ind = find(P == 0);
if ~isempty(ind)
    warning('P had some weird values, PD is symmetrized')
    sym = 1;
end

if sym == 1
    % symmetrize
    rho         = rho(1:end-1,1:end-1);
    rho_fluid   = rho_fluid(1:end-1,1:end-1);
    melt        = melt(1:end-1,1:end-1);
    T_C         = T_C(1:end-1,1:end-1);
    P           = P(1:end-1,1:end-1);
end

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




















