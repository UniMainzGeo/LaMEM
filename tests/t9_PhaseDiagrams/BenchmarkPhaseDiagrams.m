%% Test phase diagrams LaMEM
% This file plots the phase diagram
clear
close

% Indicate name
PD_name = 'TestPD.in';

sym=0;

H = dlmread(PD_name, ' ', [49 0 54 0]);
np = H(6);
nt = H(3);
M = dlmread(PD_name, ' ', 55, 0);

rho         =   reshape(M(:,3),np,nt);
rho_fluid   =   reshape(M(:,1),np,nt);
melt        =   reshape(M(:,2),np,nt);
T           =   reshape(M(:,4),np,nt);       % degree C
P           =   reshape(M(:,5),np,nt)./1e3;  % kbar  - LaMEM input must be bar!

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
ind = find(T == 0);
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
    T           = T(1:end-1,1:end-1);
    P           = P(1:end-1,1:end-1);
end

figure(1),clf
subplot(221)
pcolor(T,P,rho)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('density rock w/out melt!')

subplot(222)
pcolor(T,P,melt)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('Melt content')

subplot(223)
pcolor(T,P,rho_fluid)
ylabel('P [kbar]')
xlabel('T')
title('density of fluid')
shading interp
colorbar

subplot(224)
rho_average = (1-melt).*rho + melt.*rho_fluid;
pcolor(T,P,rho_average)
ylabel('P [kbar]')
xlabel('T')
title('Combined density of fluid + solid (should NOT be passed to LaMEM!)')
shading interp
colorbar




















