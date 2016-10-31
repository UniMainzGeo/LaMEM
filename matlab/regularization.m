
clear 

eta   = 1e25; 
tau_y = 1e8;
emin  = 1e-18;
emax  = 1e-16;
n     = 1000;

cf_eta_min = 10;
n_pw       = 100;

% set strain rate range
eps   = emin:((emax-emin)/n):emax;

% get visco-plastic transition strain rate
%eps_y = tau_y/eta/2;


% plastic viscosity
plast = tau_y./eps./2;

% quasi-harmonic averaging
hmean = 1./(1./eta + 1./plast);

% power-law approximation
cf = 1/n_pw - 1;
etapw = tau_y./2.*eps.^cf;


% visco-plastic approximation
etast = eta/cf_eta_min;
etavp = (etast + plast)/(1 + etast/eta);

figure(1), clf
hold on;

xp = [emin, emax];
yp = [eta, eta];

plot(xp,  yp,    '-k');
plot(eps, plast, '-r'); 
plot(eps, hmean, '-b'); 
plot(eps, etapw, '-g'); 
plot(eps, etavp, '-m'); 

set(gca, 'xscale', 'log');

xlabel('Strain Rate')
ylabel('Efective viscosity')
legend('Linear','Plastic','Harmonic','Power-Law','Visco-Plastic')


title('Plasticity Regularization')
box on
grid on
ylim([0 3*eta])



% figure(2), clf
% hold on;
% 
% xp = [emin, emax];
% 
% plot(eps, H, '-k'); 
% 
% set(gca, 'xscale', 'log');
% 
% xlabel('Strain Rate')
% ylabel('Heavyside function')
% legend('Heavyside')
% title('Plasticity Regularization')
% box on
% grid on
% 
% 
% figure(3), clf
% hold on;
% 
% xp = [emin, emax];
% yp = [eta, eta];
% 
% plot(xp,  yp,    '-k');
% plot(eps, plast, '-r'); 
% plot(eps, hmean, '-b'); 
% 
% set(gca, 'xscale', 'log');
% 
% xlabel('Strain Rate')
% ylabel('Efective viscosity')
% legend('Linear','Plastic','Quasi-Harmonic')
% title('Quasi-Harmonic viscosity averaging')
% box on
% grid on
% ylim([0 3*eta])

