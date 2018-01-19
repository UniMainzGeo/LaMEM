%% Test phase diagrams LaMEM
clear
close

PD_name = 'TestPD.in';


H = dlmread(PD_name, ' ', [50 0 55 1]);
np = H(5);
nt = H(2);
M = dlmread(PD_name, ' ', 56, 0);

rho = reshape(M(:,3),np,nt);
rho_fluid = reshape(M(:,1),np,nt);
melt = reshape(M(:,2),np,nt);
T = reshape(M(:,4),np,nt);
P = reshape(M(:,5),np,nt)./1e3;

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
    rho = rho(1:end-1,1:end-1);
    rho_fluid = rho_fluid(1:end-1,1:end-1);
    melt = melt(1:end-1,1:end-1);
    T = T(1:end-1,1:end-1);
    P = P(1:end-1,1:end-1);
end

figure(10),clf
pcolor(P,T,rho)
shading interp
colorbar

figure(11),clf
pcolor(P,T,melt)
shading interp
colorbar

figure(12),clf
pcolor(P,T,rho_fluid)
shading interp
colorbar

%% Interpolate
n = 100;
P_int = linspace(0,2500,n)./1e3; %kbar
T_int = linspace(0,300,n); % C
x = linspace(0,10,n);

% find hotter block
ind = find(x<7.5 & x > 2.5);
T_int(ind) = 800;

T_int = T_int+273.15; % C -> K

hold on
figure(10)
line(P_int,T_int,'LineWidth',4,'Color','r')
figure(11)
line(P_int,T_int,'LineWidth',4,'Color','r')
figure(12)
line(P_int,T_int,'LineWidth',4,'Color','r')
hold off

% interpolate along line
rho_int = interp2(P,T,rho,P_int,T_int);
melt_int = interp2(P,T,melt,P_int,T_int);
rho_fluid_int = interp2(P,T,rho_fluid,P_int,T_int);

rho_int_final = (melt_int.*rho_fluid_int) + ((1-melt_int).*rho_int);

figure(20),clf
subplot(1,2,1)
plot(x,rho_int_final);
subplot(1,2,2)
plot(x,melt_int)


abc = 1;



























